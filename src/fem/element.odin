// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Static reference element topology and geometric mappings from reference to physical.
 Everything here is element-local, connectivity concerns are part of the mesh.

 ASSUMPTIONS:
  - All geometry is assumed to be embedded in 3D space.
    This simplifies geometry code while allowing for surface FEM too.
  - Standard FEM node ordering is used.

 NOTES:
  - The term "facet" is used to describe the face of an element in any dimension.
    A line on a 2d mesh, or a point on a 1d mesh.

 SUPPORT:
  - Common element types at linear and quadratic order.
*/
package fem

import "core:log"
import "core:math/linalg"
import "core:slice"

Vec3 :: linalg.Vector3f64
Mat3 :: linalg.Matrix3x3f64

AMBIENT_DIM :: 3

Element_Type :: enum u8 {
	Point,
	Line,
	Triangle,
	Quadrilateral,
	Tetrahedron,
}

Order :: enum u8 {
	Linear,
	Quadratic,
}

Dimension :: enum {
	D0,
	D1,
	D2,
	D3,
}

// Instantiation of a particular reference element (type, order) at a particular physical location (nodes)
Element :: struct {
	type:   Element_Type,
	order:  Order,
	affine: bool, // quads and hexes may be linear order, but non-affine, hence the seperate flag.
	nodes:  []Vec3,
}

Reference_Element :: struct {
	dimension:           Dimension,
	sub_entity_nodes:    [Order][Dimension][][]int,
	facet_element_types: []Element_Type,
	facet_tangents:      [][]Vec3,
}


@(rodata)
REFERENCE_ELEMENTS := [Element_Type]Reference_Element {
	.Point         = REFERENCE_POINT,
	.Line          = REFERENCE_LINE,
	.Triangle      = REFERENCE_TRIANGLE,
	.Quadrilateral = REFERENCE_QUADRILATERAL,
	.Tetrahedron   = REFERENCE_TETRAHEDRON,
}

// Topology (Decided to wrap in access functions in case of table changes or safety checks)

MAX_FACETS :: 6 //enough for hexahedrons.

element_dim :: #force_inline proc(type: Element_Type) -> Dimension {
	return REFERENCE_ELEMENTS[type].dimension
}

element_sub_entity_nodes :: #force_inline proc(type: Element_Type, order: Order, entity_dimension: Dimension) -> [][]int {
	assert(type != .Point)
	assert(entity_dimension < element_dim(type))
	return REFERENCE_ELEMENTS[type].sub_entity_nodes[order][entity_dimension]
}

element_facet_type :: #force_inline proc(parent_type: Element_Type, facet: int) -> Element_Type {
	assert(parent_type != .Point)
	return REFERENCE_ELEMENTS[parent_type].facet_element_types[facet]
}

element_num_facets :: #force_inline proc(parent_type: Element_Type) -> int {
	assert(parent_type != .Point)
	return len(REFERENCE_ELEMENTS[parent_type].facet_element_types)
}

element_facet_dim :: #force_inline proc(parent_type: Element_Type) -> Dimension {
	assert(parent_type != .Point)
	return Dimension(int(element_dim(parent_type)) - 1)
}

element_facet_tangents :: #force_inline proc(parent_type: Element_Type, facet: int) -> []Vec3 {
	assert(parent_type != .Point)
	return REFERENCE_ELEMENTS[parent_type].facet_tangents[facet]
}

// used by visualization code to get rid of some issues when splitting higher-order cells
// into many smaller linear sub cells.
element_reduce_to_linear :: proc(element: Element) -> (reduced: Element) {
	reduced = element
	reduced.order = .Linear
	reduced.nodes = element.nodes[:len(element_sub_entity_nodes(element.type, .Linear, .D0))]
	return reduced
}

// This query currently does not have enough table information to be a quick lookup.
// It is expected to be used outside of hot paths, for sparsity pattern generation or similar.
// If perf becomes an issue, this logic can be fully, or partially moved into reference element tables.
element_entity_on_facet :: proc(parent_type: Element_Type, facet: int, entity_dim: Dimension, entity_index: int) -> bool {
	assert(parent_type != .Point)
	assert(entity_dim <= element_facet_dim(parent_type))

	facet_nodes := element_sub_entity_nodes(parent_type, .Linear, element_facet_dim(parent_type))[facet]
	entity_nodes := element_sub_entity_nodes(parent_type, .Linear, entity_dim)[entity_index]

	for node in entity_nodes {
		if !slice.contains(facet_nodes, node) {return false}
	}

	return true
}

// Geometry

// Maps a reference point to a physical location.
compute_physical_point :: proc(element: Element, ref_point: Vec3) -> (p: Vec3) {
	for node, idx in element.nodes {
		p += node * raw_shape_value(element.type, element.order, idx, ref_point)
	}
	return p
}

// Same as above, prefer this one if sample is known.
compute_physical_point_sample :: proc(element: Element, s: Sample) -> (p: Vec3) {
	for node, idx in element.nodes {
		val := LAGRANGE_REFERENCE_TABLES[element.type][element.order][s.group].values[s.index][idx]
		p += val * node
	}
	return p
}

// Computes the jacobian which is the "derivative" of coordinate transform at some reference point
compute_jacobian :: proc(element: Element, ref_point: Vec3) -> (j: Mat3) {
	for node, idx in element.nodes {
		j += linalg.outer_product(node, raw_shape_gradient(element.type, element.order, idx, ref_point))
	}
	return j
}

// Same as above, prefer if sample is known.
compute_jacobian_sample :: proc(element: Element, s: Sample) -> (j: Mat3) {
	for node, idx in element.nodes {
		grad := LAGRANGE_REFERENCE_TABLES[element.type][element.order][s.group].gradients[s.index][idx]
		j += linalg.outer_product(node, grad)
	}
	return j
}

// Computes the jacobian of an individual facet of the element.
compute_facet_jacobian :: proc(element: Element, parent_jacobian: Mat3, facet: int) -> (j: Mat3) {
	for tangent, i in element_facet_tangents(element.type, facet) {j[i] = parent_jacobian * tangent}
	return j
}

// Returns an outward normal of a given facet, requires both the facet and parent jacobians.
compute_facet_normal :: proc(facet_element: Element_Type, facet_jacobian: Mat3, parent_jacobian: Mat3) -> Vec3 {
	assert(element_dim(facet_element) != .D3)
	switch element_dim(facet_element) {
	case .D0:
		return linalg.normalize(parent_jacobian[0])
	case .D1:
		edge_tangent := linalg.normalize(facet_jacobian[0])
		surface_normal := linalg.normalize(linalg.cross(parent_jacobian[0], parent_jacobian[1]))
		return linalg.normalize(linalg.cross(edge_tangent, surface_normal))
	case .D2:
		normal := linalg.cross(facet_jacobian[0], facet_jacobian[1])
		return linalg.normalize(normal)
	case .D3:
		unreachable()
	case:
		unreachable()
	}
}

// works for 1d elements only.
element_tangent :: proc(element: Element_Type, jacobian: Mat3) -> Vec3 {
	assert(element_dim(element) == .D1)
	return linalg.normalize(jacobian[0])
}

// Returns the determinant, or determinant-like quantity, from the given jacobian.
// The jacobian must have been constructed from an element with dimension == `dim`.
measure_from_jacobian :: proc(dim: Dimension, jacobian: Mat3) -> f64 {
	switch dim {
	case .D0:
		return 1
	case .D1:
		return linalg.length(jacobian[0])
	case .D2:
		return linalg.length(linalg.cross(jacobian[0], jacobian[1]))
	case .D3:
		return linalg.determinant(jacobian)
	case:
		unreachable()
	}
}

// Returns the inverse, or inverse-like quantity from the given jacobian.
// The jacobian must have been constructed from an element with dimension == `dim`.
pullback_from_jacobian :: proc(dim: Dimension, jacobian: Mat3) -> (transform: Mat3) {
	switch dim {
	case .D0:
		return {}
	case .D1:
		len_sq := linalg.length2(jacobian[0])
		transform[0] = jacobian[0] / len_sq
		return linalg.transpose(transform)
	case .D2:
		J: matrix[3, 2]f64
		J[0] = jacobian[0]
		J[1] = jacobian[1]
		G_inv := linalg.inverse(linalg.transpose(J) * J)
		mat2x3 := linalg.transpose(J * G_inv)
		transform = matrix[3, 3]f64{
			mat2x3[0, 0], mat2x3[0, 1], mat2x3[0, 2],
			mat2x3[1, 0], mat2x3[1, 1], mat2x3[1, 2],
			0, 0, 0,
		}
		return linalg.transpose(transform)
	case .D3:
		return linalg.transpose(linalg.inverse(jacobian))
	case:
		unreachable()
	}
}


// Reference element definitions

@(private = "file")
REFERENCE_POINT: Reference_Element : {dimension = .D0}

@(private = "file")
REFERENCE_LINE: Reference_Element : {
	dimension = .D1,
	facet_element_types = {.Point, .Point},
	sub_entity_nodes = {.Linear = #partial{.D0 = {{0}, {1}}}, .Quadratic = #partial{.D0 = {{0}, {1}}}},
}

@(private = "file")
REFERENCE_TRIANGLE: Reference_Element : {
	dimension = .D2,
	facet_element_types = {.Line, .Line, .Line},
	facet_tangents = {{{1, 0, 0}}, {{-1, 1, 0}}, {{0, -1, 0}}},
	sub_entity_nodes = {
		.Linear = #partial{.D0 = {{0}, {1}, {2}}, .D1 = {{0, 1}, {1, 2}, {2, 0}}},
		.Quadratic = #partial{.D0 = {{0}, {1}, {2}}, .D1 = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}}},
	},
}

@(private = "file")
REFERENCE_QUADRILATERAL: Reference_Element : {
	dimension = .D2,
	facet_element_types = {.Line, .Line, .Line, .Line},
	facet_tangents = {{{1, 0, 0}}, {{0, 1, 0}}, {{-1, 0, 0}}, {{0, -1, 0}}},
	sub_entity_nodes = {
		.Linear = #partial{.D0 = {{0}, {1}, {2}, {3}}, .D1 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}}},
		.Quadratic = #partial{.D0 = {{0}, {1}, {2}, {3}}, .D1 = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}}},
	},
}

@(private = "file")
REFERENCE_TETRAHEDRON: Reference_Element : {
	dimension = .D3,
	facet_element_types = {.Triangle, .Triangle, .Triangle, .Triangle},
	facet_tangents = {{{0, 1, 0}, {1, 0, 0}}, {{1, 0, 0}, {0, 0, 1}}, {{-1, 1, 0}, {-1, 0, 1}}, {{0, 0, 1}, {0, 1, 0}}},
	sub_entity_nodes = {
		.Linear = #partial{
			.D0 = {{0}, {1}, {2}, {3}},
			.D1 = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}},
			.D2 = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}},
		},
		.Quadratic = #partial{
			.D0 = {{0}, {1}, {2}, {3}},
			.D1 = {{0, 1, 4}, {1, 2, 5}, {2, 0, 6}, {0, 3, 7}, {1, 3, 8}, {2, 3, 9}},
			.D2 = {{0, 1, 2, 4, 5, 6}, {0, 1, 3, 4, 8, 7}, {1, 2, 3, 5, 9, 8}, {2, 0, 3, 6, 7, 9}},
		},
	},
}
