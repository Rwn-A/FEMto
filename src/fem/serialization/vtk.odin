package serialization

import "core:bufio"
import sa "core:container/small_array"
import "core:encoding/base64"
import "core:fmt"
import "core:io"
import "core:mem"
import "core:os"
import "core:strconv"

import fem "../"

MAX_OUTPUT_FIELD_COMPONENTS :: 9

Output_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: [][MAX_OUTPUT_FIELD_COMPONENTS]f64)

Output_Field :: struct {
	friendly_name:  string,
	components:     int,
	data:           rawptr,
	value_provider: Output_Proc,
}


VTK_Raw_Mesh :: struct {
	vertices:     []fem.Vec3,
	connectivity: []i32,
	offsets:      []i32,
	types:        []u8,
	output_order: fem.Basis_Order,
}

VTK_Element_Type :: enum u8 {
	Point         = 1,
	Line          = 3,
	Triangle      = 5,
	Quadrilateral = 9,
	Tetrahedron   = 10,
	Hexahedron    = 12,
}

@(rodata)
VTK_ELEMENT_TYPE_FROM_NATIVE := [fem.Element_Type]VTK_Element_Type {
	.Point         = .Point,
	.Line          = .Line,
	.Triangle      = .Triangle,
	.Quadrilateral = .Quadrilateral,
	.Tetrahedron   = .Tetrahedron,
	.Hexahedron    = .Hexahedron,
}

Output_Variable_Data :: struct {
	system: fem.System,
	var:    fem.Var_Handle,
	data:   fem.Vector,
}

output_field_from_system_variable :: proc(
	$SPACE_TYPE: typeid,
	data: ^Output_Variable_Data,
	name: string,
) -> (
	of: Output_Field,
) {
	of.friendly_name = name
	of.components = 1 when SPACE_TYPE == fem.Grad_Space(.Scalar) else 3
	of.data = data

	of.value_provider = proc(
		mapped: fem.Mapped_Element,
		time: f64,
		data: rawptr,
		out: [][MAX_OUTPUT_FIELD_COMPONENTS]f64,
	) {
		od := cast(^Output_Variable_Data)data
		bd := fem.system_var_bd(od.system, od.var)

		when SPACE_TYPE == fem.Grad_Space(.Scalar) {
			space := fem.basis_grad_space(mapped, bd, .Scalar)
			for p in 0 ..< fem.space_points(space) {
				for dof in 0 ..< fem.space_arity(space) {
					coeff := od.data[fem.system_global_dof(od.system, od.var, mapped.element.id, fem.Local_DOF(dof))]
					out[p][0] += coeff * fem.space_value(space, p, dof)
				}
			}
		} else when SPACE_TYPE == fem.Grad_Space(.Vector) {
			space := fem.basis_grad_space(mapped, bd, .Vector)
			for p in 0 ..< fem.space_points(space) {
				for dof in 0 ..< fem.space_arity(space) {
					coeff := od.data[fem.system_global_dof(od.system, od.var, mapped.element.id, fem.Local_DOF(dof))]
					v := coeff * fem.space_value(space, p, dof)
					out[p][0] += v.x
					out[p][1] += v.y
					out[p][2] += v.z
				}
			}
		} else {
			#assert(false, "unimplemented space in output proc")
		}
	}

	return of

}


vtk_create_visualization_mesh :: proc(
	mesh: fem.Mesh,
	order: fem.Basis_Order,
	allocator := context.allocator,
) -> VTK_Raw_Mesh {
	context.allocator = allocator

	vertices := make([dynamic]fem.Vec3)
	connectivity := make([dynamic]i32)
	offsets := make([dynamic]i32)
	types := make([dynamic]u8)

	for el, i in mesh.elements {
		base_vertex_idx := i32(len(vertices))

		element: fem.Mesh_Element = el
		element.el = fem.element_reduce_to_linear(element)
		subcell := fem.map_element(element, &fem.SUBCELL_POINT_RULES[element.type][order])
		defer fem.mapped_destroy(&subcell)

		for i in 0 ..< len(subcell.rule.points) {
			append(&vertices, subcell.im.point[i])
		}


		for subcell in fem.REFERENCE_ELEMENTS[element.type].subcell[order].connectivity {
			for node_index in subcell {append(&connectivity, base_vertex_idx + i32(node_index))}
			append(&offsets, i32(len(connectivity)))
			append(&types, u8(VTK_ELEMENT_TYPE_FROM_NATIVE[element.type]))
		}
	}
	return {
		vertices = vertices[:],
		connectivity = connectivity[:],
		offsets = offsets[:],
		types = types[:],
		output_order = order,
	}
}

vtk_destroy_visualization_mesh :: proc(mesh: VTK_Raw_Mesh, allocator := context.allocator) {
	context.allocator = allocator
	delete(mesh.vertices)
	delete(mesh.connectivity)
	delete(mesh.offsets)
	delete(mesh.types)
}

// xml-based vtk file specific

@(private = "file")
XML_Tag :: enum {
	VTKFile,
	UnstructuredGrid,
	Piece,
	Points,
	DataArray,
	Cells,
	CellData,
	PointData,
	Collection,
	DataSet,
	AppendedData,
}

@(private = "file")
MAX_NESTED_TAGS :: 16

@(private = "file")
XML_Writer :: struct {
	open: sa.Small_Array(MAX_NESTED_TAGS, XML_Tag),
	w:    io.Stream,
}

@(thread_local, private = "file")
xml_writer_buffer: [2048]u8

@(private = "file")
XML_Attribute :: struct {
	key:   string,
	value: string,
}


@(private = "file")
xml_open_tag :: proc(xml_w: ^XML_Writer, tag: XML_Tag, attrib: []XML_Attribute = {}) {
	sa.push_back(&xml_w.open, tag)
	fmt.wprintf(xml_w.w, "<%s", tag)
	defer fmt.wprintf(xml_w.w, ">")
	for a in attrib {fmt.wprintf(xml_w.w, " %s = \"%s\"", a.key, a.value)}
}

@(private = "file")
xml_close_tag :: proc(xml_w: ^XML_Writer) {
	fmt.wprintfln(xml_w.w, "</%s>", sa.pop_back(&xml_w.open))
}

// vtu specifc


write_vtu :: proc(
	path: string,
	mesh: fem.Mesh,
	viz_mesh: VTK_Raw_Mesh,
	time: f64,
	fields: []Output_Field,
	allocator := context.allocator,
) -> (
	err: os.Error,
) {
	context.allocator = allocator

	xml_w: XML_Writer
	fd := os.open(path, {.Read, .Write, .Create, .Trunc}) or_return
	defer os.close(fd)

	xml_w.w = io.to_writer(os.to_writer(fd))

	fmt.wprintln(xml_w.w, "<?xml version=\"1.0\"?>")

	xml_open_tag(
		&xml_w,
		.VTKFile,
		{{"type", "UnstructuredGrid"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}, {"header_type", "UInt64"}},
	)
	defer xml_close_tag(&xml_w)

	b0, b1: [32]u8
	num_points := strconv.write_int(b0[:], cast(i64)len(viz_mesh.vertices), 10)
	num_cells := strconv.write_int(b1[:], cast(i64)len(viz_mesh.offsets), 10)

	xml_open_tag(&xml_w, .UnstructuredGrid)
	xml_open_tag(&xml_w, .Piece, {{"NumberOfPoints", num_points}, {"NumberOfCells", num_cells}})

	appended := make([dynamic]u8)
	defer delete(appended)

	appended_write :: proc(appended: ^[dynamic]u8, data: []u8) -> (offset: string) {
		@(static) b: [32]u8
		offset = strconv.write_int(b[:], cast(i64)len(appended), 10)
		len_bytes := transmute([8]u8)uint(len(data))
		append(appended, ..(len_bytes[:]))
		append(appended, ..data)
		return offset
	}


	{
		xml_open_tag(&xml_w, .Points)
		defer xml_close_tag(&xml_w)

		offset := appended_write(&appended, mem.slice_to_bytes(viz_mesh.vertices))

		xml_open_tag(
			&xml_w,
			.DataArray,
			{{"type", "Float64"}, {"NumberOfComponents", "3"}, {"format", "appended"}, {"offset", offset}},
		)
		defer xml_close_tag(&xml_w)
	}

	{
		xml_open_tag(&xml_w, .Cells)
		defer xml_close_tag(&xml_w)

		offset_connectivity := appended_write(&appended, mem.slice_to_bytes(viz_mesh.connectivity))

		xml_open_tag(
			&xml_w,
			.DataArray,
			{{"type", "Int32"}, {"Name", "connectivity"}, {"format", "appended"}, {"offset", offset_connectivity}},
		)
		xml_close_tag(&xml_w)

		offset_offsets := appended_write(&appended, mem.slice_to_bytes(viz_mesh.offsets))
		xml_open_tag(
			&xml_w,
			.DataArray,
			{{"type", "Int32"}, {"Name", "offsets"}, {"format", "appended"}, {"offset", offset_offsets}},
		)
		xml_close_tag(&xml_w)

		offset_types := appended_write(&appended, mem.slice_to_bytes(viz_mesh.types))
		xml_open_tag(
			&xml_w,
			.DataArray,
			{{"type", "UInt8"}, {"Name", "types"}, {"format", "appended"}, {"offset", offset_types}},
		)
		xml_close_tag(&xml_w)
	}

	// actual field data
	{
		xml_open_tag(&xml_w, .PointData)
		defer xml_close_tag(&xml_w)

		point_data := make([][dynamic]f64, len(fields))
		for i in 0 ..< len(fields) {point_data[i] = make([dynamic]f64)}
		defer {
			delete(point_data)
			for i in 0 ..< len(fields) {delete(point_data[i])}
		}

		for el, element_id in mesh.elements {
			element: fem.Mesh_Element = el
			element.el = fem.element_reduce_to_linear(element)
			subcell := fem.map_element(element, &fem.SUBCELL_POINT_RULES[element.type][viz_mesh.output_order])
			defer fem.mapped_destroy(&subcell)


			out := make([][MAX_OUTPUT_FIELD_COMPONENTS]f64, len(subcell.rule.points))
			defer delete(out)

			for field, i in fields {
				mem.zero_slice(out)
				field.value_provider(subcell, time, field.data, out)
				for &point in out {
					append(&point_data[i], ..(point[:field.components]))
				}
			}
		}

		for field, i in fields {
			data_offset := appended_write(&appended, mem.slice_to_bytes(point_data[i][:]))
			str_components := strconv.write_int(b0[:], cast(i64)field.components, 10)
			xml_open_tag(
				&xml_w,
				.DataArray,
				{
					{"type", "Float64"},
					{"Name", field.friendly_name},
					{"format", "appended"},
					{"NumberOfComponents", str_components},
					{"offset", data_offset},
				},
			)
			xml_close_tag(&xml_w)
		}

	}

	xml_close_tag(&xml_w) // close piece
	xml_close_tag(&xml_w) // close unstructered grid

	xml_open_tag(&xml_w, .AppendedData, {{"encoding", "raw"}})
	fmt.wprintf(xml_w.w, "_")
	io.write(xml_w.w, appended[:])
	xml_close_tag(&xml_w)

	return nil
}

write_pvd :: proc(path: string, vtk_files: []string, times: []f64) -> os.Error {
	xml_w: XML_Writer
	fd := os.open(path, {.Read, .Write, .Create, .Trunc}) or_return
	defer os.close(fd)

	xml_w.w = io.to_writer(os.to_writer(fd))

	fmt.wprintln(xml_w.w, "<?xml version=\"1.0\"?>")

	xml_open_tag(&xml_w, .VTKFile, {{"type", "Collection"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}})
	defer xml_close_tag(&xml_w)

	xml_open_tag(&xml_w, .Collection)
	defer xml_close_tag(&xml_w)

	for entry in soa_zip(vtk_file = vtk_files, time = times) {
		b: [32]u8
		time_str := strconv.write_float(b[:], entry.time, 'f', 8, 64)
		xml_open_tag(
			&xml_w,
			.DataSet,
			{{"timestep", time_str}, {"group", ""}, {"part", "0"}, {"file", entry.vtk_file}},
		)
		xml_close_tag(&xml_w)
	}

	return nil

}
