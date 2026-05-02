// Simple mesh structure, assumes conforming meshes.
package fem

import "core:math/linalg"
import "core:slice"

Boundary_ID :: distinct u32
Entity_ID :: distinct u32
Section_ID :: distinct u32

Mesh :: struct {
	elements:                                      []Mesh_Element,
	encountered_elements:                          bit_set[Element_Type],
	num_vertices, num_edges, num_faces, num_cells: int,
	boundary_names:                                map[string]Boundary_ID,
	section_names:                                 map[string]Section_ID,
	dim:                                           Dimension,
}

Mesh_Element :: struct {
	using el:              Element,
	downward_connectivity: [Dimension][]Entity_ID,
	boundary_facets:       bit_set[0 ..< MAX_FACETS],
	boundary_ids:          [MAX_FACETS]Boundary_ID, // only indices in the above set valid.
	section:               Section_ID,
	edge_orientation:      []bool,
	id:                    Entity_ID,
	// TODO: face orientation
}

// TEMPORARY, minorily improves threaded assembly performance.
spatially_order_mesh :: proc(mesh: ^Mesh) {
	if len(mesh.elements) == 0 {return}

	bb_min := mesh.elements[0].nodes[0]
	bb_max := bb_min
	for el in mesh.elements {
		c := element_centroid(el)
		bb_min = linalg.min(bb_min, c)
		bb_max = linalg.max(bb_max, c)
	}

	extent := bb_max - bb_min
	inv := Vec3 {
		extent.x > 0 ? 1.0 / extent.x : 0,
		extent.y > 0 ? 1.0 / extent.y : 0,
		extent.z > 0 ? 1.0 / extent.z : 0,
	}

	Key_Elem :: struct {
		key: u64,
		idx: int,
	}
	keys := make([]Key_Elem, len(mesh.elements))
	defer delete(keys)

	for el, i in mesh.elements {
		c := element_centroid(el)
		nx := u32((c.x - bb_min.x) * inv.x * 0x1fffff)
		ny := u32((c.y - bb_min.y) * inv.y * 0x1fffff)
		nz := u32((c.z - bb_min.z) * inv.z * 0x1fffff)
		keys[i] = {morton3d(nx, ny, nz), i}
	}

	slice.sort_by(keys, proc(a, b: Key_Elem) -> bool {return a.key < b.key})

	visited := make([]bool, len(mesh.elements))
	defer delete(visited)

	for start in 0 ..< len(mesh.elements) {
		if visited[start] || keys[start].idx == start {continue}
		tmp := mesh.elements[start]
		j := start
		for {
			next := keys[j].idx
			visited[j] = true
			if next == start {mesh.elements[j] = tmp; break}
			mesh.elements[j] = mesh.elements[next]
			j = next
		}
	}

	morton3d :: proc(x, y, z: u32) -> u64 {
		expand :: proc(v: u32) -> u64 {
			n := u64(v) & 0x1fffff
			n = (n | n << 32) & 0x1f00000000ffff
			n = (n | n << 16) & 0x1f0000ff0000ff
			n = (n | n << 8) & 0x100f00f00f00f00f
			n = (n | n << 4) & 0x10c30c30c30c30c3
			n = (n | n << 2) & 0x1249249249249249
			return n
		}
		return expand(x) | (expand(y) << 1) | (expand(z) << 2)
	}

	element_centroid :: proc(el: Element) -> (c: Vec3) {
		for n in el.nodes {c += n}
		return c * (1.0 / f64(len(el.nodes)))
	}
}
