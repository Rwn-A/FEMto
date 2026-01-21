// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Connectivity between elements, and their sub-entities.
*/
package fem

Boundary_ID :: distinct uint
Entity_ID :: distinct u32 // just to cut down on peak memory in mesh storage.
Section_ID :: distinct u32

Mesh :: struct {
	elements:                                      []Mesh_Element, // primary element in the mesh, lines in 1d, quads/tris in 2d etc.
	num_vertices, num_edges, num_faces, num_cells: int,
	boundary_names:                                map[string]Boundary_ID,
	section_names:  							   map[string]Section_ID,
	dim:                                           Dimension,
}

Mesh_Element :: struct {
	using el:              Element,
	downward_connectivity: [Dimension][]Entity_ID,
	adjacency:             []union {
		Boundary_ID,
		Entity_ID,
	},
	section: Section_ID,
	edge_orientation:      []bool,
	face_orientation:      []u8,
}

FACE_ORIENTATIONS := #partial [Element_Type][][]int {
	.Triangle      = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}},
	.Quadrilateral = {
		{0, 1, 2, 3},
		{1, 2, 3, 0},
		{2, 3, 0, 1},
		{3, 0, 1, 2},
		{0, 3, 2, 1},
		{3, 2, 1, 0},
		{2, 1, 0, 3},
		{1, 0, 3, 2},
	},
}
