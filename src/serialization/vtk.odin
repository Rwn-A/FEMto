// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package serialization

import "core:bufio"
import sa "core:container/small_array"
import "core:encoding/base64"
import "core:fmt"
import "core:io"
import "core:mem"
import "core:os/os2"
import "core:strconv"

import fem "../fe_core"

VTK_Raw_Mesh :: struct {
	vertices:     []fem.Vec3,
	connectivity: []i32,
	offsets:      []i32,
	types:        []u8,
	rule:         fem.Subcell_Rule,
}

VTK_Element_Type :: enum u8 {
	Point         = 1,
	Line          = 3,
	Triangle      = 5,
	Quadrilateral = 9,
	Tetrahedron   = 10,
}

@(rodata)
VTK_ELEMENT_TYPE_FROM_NATIVE := [fem.Element_Type]VTK_Element_Type {
	.Point         = .Point,
	.Line          = .Line,
	.Triangle      = .Triangle,
	.Quadrilateral = .Quadrilateral,
	.Tetrahedron   = .Tetrahedron,
}

vtk_create_visualization_mesh :: proc(mesh: fem.Mesh, rule: fem.Subcell_Rule, allocator := context.allocator) -> VTK_Raw_Mesh {
	context.allocator = allocator

	vertices := make([dynamic]fem.Vec3)
	connectivity := make([dynamic]i32)
	offsets := make([dynamic]i32)
	types := make([dynamic]u8)

	for &element, i in mesh.elements {
		base_vertex_idx := i32(len(vertices))

		subcell, count := fem.subcell_for(element, rule, {.Physical_Point})
		subcell.element.el = fem.element_reduce_to_linear(element)
		for i in 0 ..< count {
			append(&vertices, fem.compute_physical_point_context(subcell, i))
		}

		for subcell in fem.SUBCELL_CONNECTIVITY[element.type][rule] {
			for node_index in subcell {append(&connectivity, base_vertex_idx + i32(node_index))}
			append(&offsets, i32(len(connectivity)))
			append(&types, u8(VTK_ELEMENT_TYPE_FROM_NATIVE[element.type]))
		}
	}
	return {vertices = vertices[:], connectivity = connectivity[:], offsets = offsets[:], types = types[:], rule = rule}
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
	fields: []Output_Field,
	allocator := context.allocator,
) -> (
	err: os2.Error,
) {
	context.allocator = allocator

	xml_w: XML_Writer
	fd := os2.open(path, {.Read, .Write, .Create, .Trunc}) or_return
	defer os2.close(fd)

	xml_w.w = io.to_writer(os2.to_writer(fd))

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

		for element, element_id in mesh.elements {
			geometry_required := fem.Geometry_Options{}
			for field in fields {geometry_required += field.geometry_required}

			subcell, count := fem.subcell_for(element, viz_mesh.rule, geometry_required)

			for point in 0 ..< count {
				for field, i in fields {
					data := field.value_provider(subcell, point, field.data)
					append(&point_data[i], ..(data[:field.components]))
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

write_pvd :: proc(path: string, vtk_files: []string, times: []f64) -> os2.Error {
	xml_w: XML_Writer
	fd := os2.open(path, {.Read, .Write, .Create, .Trunc}) or_return
	defer os2.close(fd)

	xml_w.w = io.to_writer(os2.to_writer(fd))

	fmt.wprintln(xml_w.w, "<?xml version=\"1.0\"?>")

	xml_open_tag(
		&xml_w,
		.VTKFile,
		{{"type", "Collection"}, {"version", "0.1"}, {"byte_order", "LittleEndian"}},
	)
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