// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
 Parses a gmsh file into our native mesh structure.

 See readme for list of mesh restrictions.

 NOTE:
  - File is not optimized, nor polished.
  - Not all mesh invariants are enforced.
*/
package serialization


import "base:intrinsics"
import "core:bufio"
import "core:io"
import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:os/os2"
import "core:slice"
import "core:strconv"
import "core:strings"

import fem "../fe_core"


GMSH_EXPECTED_MSH_VERSION :: "2.2 1 8"

Gmsh_Validation_Error :: enum {
	Parse_Error,
	Unsupported_Element,
	Incorrect_Dimension,
}

Gmsh_Error :: union {
	io.Error,
	os2.Error,
	Gmsh_Validation_Error,
}

Gmsh_Warning :: enum {
	Skipped_Element,
}

@(private = "file")
ALIGNMENT_TOLERANCE :: 1e-12

// numbers align with gmsh element type integers.
@(private = "file")
Gmsh_Element_Type :: enum {
	MSH_LINE_2     = 1,
	MSH_LINE_3     = 8,
	MSH_TRIANGLE_3 = 2,
	MSH_TRIANGLE_6 = 9,
	MSH_QUAD_4     = 3,
	MSH_QUAD_9     = 10,
	MSH_TETRA_4    = 4,
	MSH_TETRA_10   = 11,
	MSH_HEX_8      = 5,
	MSH_HEX_27     = 12,
	MSH_POINT      = 15,
}

@(private = "file")
GMSH_ELEMENT_NUM_NODES := #sparse[Gmsh_Element_Type]int {
	.MSH_LINE_2     = 2,
	.MSH_LINE_3     = 3,
	.MSH_TRIANGLE_3 = 3,
	.MSH_TRIANGLE_6 = 6,
	.MSH_QUAD_4     = 4,
	.MSH_QUAD_9     = 9,
	.MSH_TETRA_4    = 4,
	.MSH_TETRA_10   = 10,
	.MSH_HEX_8      = 8,
	.MSH_HEX_27     = 27,
	.MSH_POINT      = 1,
}


@(private = "file")
SUPPORTED_ELEMENTS :: bit_set[Gmsh_Element_Type] {
	.MSH_LINE_2,
	.MSH_LINE_3,
	.MSH_TRIANGLE_3,
	.MSH_TRIANGLE_6,
	.MSH_QUAD_4,
	.MSH_QUAD_9,
	.MSH_TETRA_4,
	.MSH_TETRA_10,
	.MSH_HEX_8,
	.MSH_HEX_27,
	.MSH_POINT,
}


@(private = "file")
Raw_Element :: struct {
	type:         Gmsh_Element_Type,
	tags:         []i32,
	node_indices: []int,
}

@(private = "file")
gmsh_type_to_element_info :: proc(type: Gmsh_Element_Type) -> (fem.Element_Type, fem.Order) {
	switch type {
	case .MSH_LINE_2:
		return .Line, .Linear
	case .MSH_LINE_3:
		return .Line, .Quadratic
	case .MSH_TRIANGLE_3:
		return .Triangle, .Linear
	case .MSH_TRIANGLE_6:
		return .Triangle, .Quadratic
	case .MSH_QUAD_4:
		return .Quadrilateral, .Linear
	case .MSH_QUAD_9:
		return .Quadrilateral, .Quadratic
	case .MSH_TETRA_4:
		return .Tetrahedron, .Linear
	case .MSH_TETRA_10:
		return .Tetrahedron, .Quadratic
	case .MSH_HEX_8:
		return .Hexahedron, .Linear
	case .MSH_HEX_27:
		return .Hexahedron, .Quadratic
	case .MSH_POINT:
		return .Point, .Linear
	case:
		unreachable()
	}
}

gmsh_parse :: proc(
	path: string,
	allocator := context.allocator,
) -> (
	mesh: fem.Mesh,
	warn: bit_set[Gmsh_Warning],
	err: Gmsh_Error,
) {
	scratch: virtual.Arena

	if err := virtual.arena_init_growing(&scratch); err != nil {log.panic("Allocation Failure: ", err)}
	defer virtual.arena_destroy(&scratch)

	context.temp_allocator = virtual.arena_allocator(&scratch)
	context.allocator = allocator

	file := os2.open(path) or_return
	defer os2.close(file)

	reader: bufio.Reader
	bufio.reader_init(&reader, os2.to_stream(file))
	defer bufio.reader_destroy(&reader)

	// format
	swap_bytes: bool
	{
		expect_ascii_line(&reader, "$MeshFormat") or_return
		expect_ascii_line(&reader, GMSH_EXPECTED_MSH_VERSION) or_return
		endianess := read_number(i32, &reader, false) or_return
		swap_bytes = endianess != 1
		bufio.reader_read_byte(&reader) or_return // newline between the above and $EndMeshFormat
		expect_ascii_line(&reader, "$EndMeshFormat") or_return
	}

	// physical names
	mesh.boundary_names = make(map[string]fem.Boundary_ID)
	mesh.section_names = make(map[string]fem.Section_ID)
	names := make(map[string]struct {
			id:  int,
			dim: fem.Dimension,
		}, context.temp_allocator)
	{
		expect_ascii_line(&reader, "$PhysicalNames") or_return
		num_groups := parse_int(read_ascii_line(&reader) or_return) or_return

		max_dimension: int
		for _ in 0 ..< num_groups {
			line := read_ascii_line(&reader) or_return

			components := strings.split(line, " ", context.temp_allocator)
			if len(components) < 3 {
				return {}, warn, .Parse_Error
			}

			group_dim := parse_int(components[0]) or_return
			if group_dim > max_dimension {max_dimension = group_dim}

			name := strings.clone(strings.trim(components[2], "\""))
			names[name] = {parse_int(components[1]) or_return, fem.Dimension(group_dim)}
		}
		expect_ascii_line(&reader, "$EndPhysicalNames") or_return

		mesh.dim = fem.Dimension(max_dimension)

		for key, value in names {
			if value.dim < mesh.dim {
				mesh.boundary_names[key] = fem.Boundary_ID(value.id)
			} else {
				mesh.section_names[key] = fem.Section_ID(value.id)
			}
		}
	}

	//nodes
	gmsh_id_to_node_id := make(map[i32]int, context.temp_allocator)
	nodes := make([dynamic]fem.Vec3, context.temp_allocator)
	{
		expect_ascii_line(&reader, "$Nodes") or_return
		num_nodes := parse_int(read_ascii_line(&reader) or_return) or_return

		for i in 0 ..< num_nodes {
			id := read_number(i32, &reader, swap_bytes) or_return
			x := read_number(f64, &reader, swap_bytes) or_return
			y := read_number(f64, &reader, swap_bytes) or_return
			z := read_number(f64, &reader, swap_bytes) or_return

			gmsh_id_to_node_id[id] = i
			append(&nodes, fem.Vec3{x, y, z})
		}
		bufio.reader_read_byte(&reader) // newline before $End
		expect_ascii_line(&reader, "$EndNodes") or_return
	}

	// TODO: tons of repeats in here, should abstract out.
	// TODO: edge/face orientation

	// raw elements
	raw_primary_elements := make([dynamic]Raw_Element, context.temp_allocator)
	raw_boundary_elements := make([dynamic]Raw_Element, context.temp_allocator)
	{
		expect_ascii_line(&reader, "$Elements") or_return
		num_elements := parse_int(read_ascii_line(&reader) or_return) or_return

		for elements_read: i32 = 0; elements_read < i32(num_elements); {
			type_integer := read_number(i32, &reader, swap_bytes) or_return
			num_in_group := read_number(i32, &reader, swap_bytes) or_return
			num_tags := read_number(i32, &reader, swap_bytes) or_return

			type := Gmsh_Element_Type(type_integer)
			if type not_in SUPPORTED_ELEMENTS {return {}, warn, .Unsupported_Element}

			// we skip elements with no tags, no physical group.
			if num_tags == 0 {
				warn += {.Skipped_Element}
				block_size := (GMSH_ELEMENT_NUM_NODES[type] * size_of(i32)) + (int(num_tags) * size_of(i32)) + size_of(i32)
				io.seek(io.to_seeker(bufio.reader_to_stream(&reader)), i64(block_size), .Current) or_return
				elements_read += num_in_group
				continue
			}

			for _ in 0 ..< num_in_group {
				elements_read += 1
				element: Raw_Element
				element.type = type
				_ = read_number(i32, &reader, swap_bytes) or_return //id

				element.tags = make([]i32, num_tags, context.temp_allocator)
				for j in 0 ..< len(element.tags) {
					element.tags[j] = read_number(i32, &reader, swap_bytes) or_return
				}

				node_count := GMSH_ELEMENT_NUM_NODES[type]
				element.node_indices = make([]int, node_count)

				for j in 0 ..< node_count {
					node_idx := read_number(i32, &reader, swap_bytes) or_return
					node_id := gmsh_id_to_node_id[node_idx]
					element.node_indices[j] = node_id
				}
				fem_type, _ := gmsh_type_to_element_info(type)
				#partial switch fem.element_dim(fem_type) {
				case fem.Dimension(int(mesh.dim) - 1):
					append(&raw_boundary_elements, element)
				case mesh.dim:
					append(&raw_primary_elements, element)
				case:
					return {}, warn, .Incorrect_Dimension
				}
			}
		}
	}

	MAX_SUB_ENTITY_VERTICES :: 4 // quadrilateral faces worst case so at most 4 vertices.

	// build mesh from raw elements

	Vertex_Key :: fem.Vec3
	Entity_Key :: [MAX_SUB_ENTITY_VERTICES]fem.Entity_ID // sorted array of vertex id's

	create_key :: proc(verts: []fem.Entity_ID) -> (key: Entity_Key) {
		copy(key[:], verts[:])
		slice.sort(key[:len(verts)])
		return key
	}

	element_is_affine :: proc(el: fem.Mesh_Element) -> bool {
		if el.order != .Linear {return false}
		#partial switch el.type {
		case .Quadrilateral:
			c1 := linalg.length(el.nodes[0] - el.nodes[1] - (el.nodes[2] - el.nodes[3])) < ALIGNMENT_TOLERANCE
			c2 := linalg.length(el.nodes[2] - el.nodes[1] - (el.nodes[3] - el.nodes[0])) < ALIGNMENT_TOLERANCE
			return c1 && c2
		case .Hexahedron:
			return false // TODO: perform the check for an affine hex
		}
		return true
	}

	found_vertices := make(map[Vertex_Key]fem.Entity_ID, context.temp_allocator)
	found_edges := make(map[Entity_Key]fem.Entity_ID, context.temp_allocator)
	found_faces := make(map[Entity_Key]fem.Entity_ID, context.temp_allocator)

	//maps a facet, to the first entity that found it.
	facet_entity_adjacency := make(map[Entity_Key]struct {
			entity:      fem.Entity_ID,
			facet_index: int,
		}, context.temp_allocator)

	next_new_vert: fem.Entity_ID
	next_new_edge: fem.Entity_ID
	next_new_face: fem.Entity_ID

	element_builder := make([dynamic]fem.Mesh_Element)

	for element, i in raw_primary_elements {
		mesh_element := fem.Mesh_Element{}
		defer append(&element_builder, mesh_element)

		mesh_element.type, mesh_element.order = gmsh_type_to_element_info(element.type)
		mesh_element.section = fem.Section_ID(element.tags[0])
		mesh_element.boundaries = make([]Maybe(fem.Boundary_ID), fem.element_num_facets(mesh_element.type))
		mesh_element.id = fem.Entity_ID(i)

		node_builder := make([dynamic]fem.Vec3)
		for node_index in element.node_indices {append(&node_builder, nodes[node_index])}
		mesh_element.nodes = node_builder[:]

		mesh_element.affine = element_is_affine(mesh_element)

		vert_connectivity_builder := make([dynamic]fem.Entity_ID)
		for vert_nodes, vert_index in fem.element_sub_entity_nodes(mesh_element.type, mesh_element.order, .D0) {
			vert_node := vert_nodes[0] // all vertices just have one node
			if val, exists := found_vertices[mesh_element.nodes[vert_node]]; exists {
				append(&vert_connectivity_builder, val)
				continue
			}
			found_vertices[mesh_element.nodes[vert_node]] = next_new_vert
			append(&vert_connectivity_builder, next_new_vert)
			if fem.element_dim(mesh_element.type) == .D1 {
				facet_entity_adjacency[create_key({next_new_vert})] = {fem.Entity_ID(i), vert_index}
			}
			next_new_vert += 1
		}

		mesh_element.downward_connectivity[.D0] = vert_connectivity_builder[:]

		if fem.element_dim(mesh_element.type) < .D2 {continue}

		// Build edges
		edge_connectivity_builder := make([dynamic]fem.Entity_ID)
		edge_orientations := make([dynamic]bool)
		for edge_nodes, edge_index in fem.element_sub_entity_nodes(mesh_element.type, mesh_element.order, .D1) {
			v0_local := edge_nodes[0]
			v1_local := edge_nodes[1]
			v0_global := found_vertices[mesh_element.nodes[v0_local]]
			v1_global := found_vertices[mesh_element.nodes[v1_local]]

			append(&edge_orientations, v0_global < v1_global)

			edge_key := create_key({v0_global, v1_global})

			if edge_key in found_edges {
				append(&edge_connectivity_builder, found_edges[edge_key])
			} else {
				found_edges[edge_key] = next_new_edge
				append(&edge_connectivity_builder, next_new_edge)
				if fem.element_dim(mesh_element.type) == .D2 {facet_entity_adjacency[edge_key] = {fem.Entity_ID(i), edge_index}}
				next_new_edge += 1
			}


		}
		mesh_element.downward_connectivity[.D1] = edge_connectivity_builder[:]
		mesh_element.edge_orientation = edge_orientations[:]

		if fem.element_dim(mesh_element.type) < .D3 {continue}

		face_connectivity_builder := make([dynamic]fem.Entity_ID)
		face_orientations := make([dynamic]u8)
		for face_nodes, face_index in fem.element_sub_entity_nodes(mesh_element.type, mesh_element.order, .D2) {
			num_face_verts := len(
				fem.element_sub_entity_nodes(fem.element_facet_type(mesh_element.type, face_index), mesh_element.order, .D0),
			)

			face_verts := make([dynamic]fem.Entity_ID, context.temp_allocator)
			for v in 0 ..< num_face_verts {
				v_local := face_nodes[v]
				v_global := found_vertices[mesh_element.nodes[v_local]]
				append(&face_verts, v_global)
			}

			face_key := create_key(face_verts[:])

			permutation := [4]int{}
			for local, local_idx in face_verts {
				for sorted, sorted_idx in face_key[:len(face_verts)] {
					if local == sorted {
						permutation[local_idx] = sorted_idx
						break
					}
				}
			}

			//facet_type := fem.element_facet_type(mesh_element.type, face_index)
			// for reference_perm, orientation_idx in fem.FACE_ORIENTATIONS[facet_type] {
			// 	if slice.equal(permutation[:], reference_perm) {
			// 		mesh_element.face_orientation[face_index] = u8(orientation_idx)
			// 		break
			// 	}
			// }


			if face_key in found_faces {
				append(&face_connectivity_builder, found_faces[face_key])
			} else {
				found_faces[face_key] = next_new_face
				append(&face_connectivity_builder, next_new_face)
				facet_entity_adjacency[face_key] = {fem.Entity_ID(i), face_index}
				next_new_face += 1
			}
		}

		mesh_element.downward_connectivity[.D2] = face_connectivity_builder[:]
	}

	for element in raw_boundary_elements {
		type, order := gmsh_type_to_element_info(element.type)
		key: Entity_Key
		switch fem.element_dim(type) {
		case .D0:
			key = {found_vertices[nodes[element.node_indices[0]]], 0, 0, 0}
		case .D1:
			key = create_key({found_vertices[nodes[element.node_indices[0]]], found_vertices[nodes[element.node_indices[1]]]})
		case .D2:
			num_face_verts := len(fem.element_sub_entity_nodes(type, order, .D0))
			face_verts := make([dynamic]fem.Entity_ID, context.temp_allocator)
			for i in 0 ..< num_face_verts {
				append(&face_verts, found_vertices[nodes[element.node_indices[i]]])
			}
			key = create_key(face_verts[:])
		case .D3:
			unreachable()
		case:
			unreachable()
		}
		assert(key in facet_entity_adjacency)
		finder := facet_entity_adjacency[key]
		element_builder[finder.entity].boundaries[finder.facet_index] = fem.Boundary_ID(element.tags[0])
	}
	mesh.num_vertices = int(next_new_vert)
	mesh.num_edges = int(next_new_edge)
	mesh.num_faces = int(next_new_face)
	mesh.num_cells = 0 if mesh.dim < .D3 else len(mesh.elements)
	mesh.elements = element_builder[:]

	return mesh, warn, nil
}

@(private = "file")
parse_int :: proc(s: string) -> (int, Gmsh_Error) {
	val, ok := strconv.parse_int(s)
	if !ok {return 0, .Parse_Error}
	return val, nil
}

@(private = "file")
expect_ascii_line :: proc(b: ^bufio.Reader, expected: string) -> Gmsh_Error {
	line := read_ascii_line(b) or_return
	if line != expected {return .Parse_Error}
	return nil
}

@(private = "file")
read_ascii_line :: proc(b: ^bufio.Reader) -> (line: string, err: Gmsh_Error) {
	line_bytes := bufio.reader_read_slice(b, '\n') or_return
	line = strings.trim_space(string(line_bytes))
	return line, nil
}

@(private = "file")
read_number :: proc(
	$T: typeid,
	b: ^bufio.Reader,
	swap_bytes: bool,
) -> (
	v: T,
	err: io.Error,
) where intrinsics.type_is_numeric(T) {
	bytes: [size_of(T)]u8
	for i in 0 ..< size_of(T) {bytes[i] = bufio.reader_read_byte(b) or_return}

	t := (cast(^T)raw_data(bytes[:]))^
	if swap_bytes {t = intrinsics.byte_swap(t)}

	return t, err
}
