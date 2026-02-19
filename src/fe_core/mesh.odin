// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Connectivity between elements, and their sub-entities.
*/
package fem

// u32 is just to cut down on memory
Boundary_ID :: distinct u32
Entity_ID :: distinct u32
Section_ID :: distinct u32

Mesh :: struct {
	elements:                                      []Mesh_Element, // primary element in the mesh, lines in 1d, quads/tris in 2d etc.
	num_vertices, num_edges, num_faces, num_cells: int,
	boundary_names:                                map[string]Boundary_ID,
	section_names:                                 map[string]Section_ID,
	dim:                                           Dimension,
}

Mesh_Element :: struct {
	using el:              Element,
	downward_connectivity: [Dimension][]Entity_ID,
	boundaries:            []Maybe(Boundary_ID), // for each facet, maybe its a boundary facet.
	section:               Section_ID,
	edge_orientation:      []bool,
	id:                    Entity_ID, // same as its index in mesh.elements
	//TODO: face_orientation
}

boundary_set :: proc(element: Mesh_Element) -> (r: bit_set[0..<MAX_FACETS]) {
	for bound, index in element.boundaries {
		if _, is := bound.?; is {r += {index}}
	}
	return r
}