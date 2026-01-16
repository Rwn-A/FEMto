// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
Global continuity helpers for FE fields.

Does not deal with local fields like those used in HDG.
*/
package fem

import "../la"

import "core:log"

DOF_Layout :: struct {
	coeffs:  []f64,
	mapping: [][]int,
}

Continuity_Method :: enum {
	CG,
	DG,
	HDG,
	EDG,
	IEDG,
}

create_layout :: proc(
	mesh: Mesh,
	method: Continuity_Method,
	$basis_type: typeid,
	order: Order,
	allocator := context.allocator,
	temp_allocator := context.temp_allocator,
) -> (
	layout: DOF_Layout,
) {
	assert(method == .CG, "Unimplemented: non-continous galerkin dof layouts.")

	context.allocator = allocator
	context.temp_allocator = temp_allocator

	Key :: struct {
		entity_id:   Entity_ID,
		basis_index: int,
	}

	vertex_map := make(map[Key]int, context.temp_allocator)
	edge_map := make(map[Key]int, context.temp_allocator)
	face_map := make(map[Key]int, context.temp_allocator)
	layout.mapping = make([][]int, len(mesh.elements), context.allocator)

	next_global: int
	for element, element_idx in mesh.elements {
		b := basis_create(basis_type, element, order)
		layout.mapping[element_idx] = make([]int, b.arity)
		for basis_index in 0 ..< b.arity {
			support := basis_support(b, basis_index)

			//interior dofs
			if support.entity_dim == element_dim(element.type) {
				layout.mapping[element_idx][basis_index] = next_global
				next_global += 1
				continue
			}

			switch support.entity_dim {
			case .D0:
				vertex_key := Key{element.downward_connectivity[.D0][support.entity_index], support.entity_basis_index}
				if vertex_key in vertex_map {
					layout.mapping[element_idx][basis_index] = vertex_map[vertex_key]
				} else {
					layout.mapping[element_idx][basis_index] = next_global
					vertex_map[vertex_key] = next_global
					next_global += 1
				}
			case .D1:
				local_global_idx := orient_edge_dof(element, support)
				edge_key := Key{element.downward_connectivity[.D1][support.entity_index], local_global_idx}
				if edge_key in edge_map {
					layout.mapping[element_idx][basis_index] = edge_map[edge_key]
				} else {
					layout.mapping[element_idx][basis_index] = next_global
					edge_map[edge_key] = next_global
					next_global += 1
				}
			case .D2:
				if support.entity_basis_arity > 0 {assert(false, "cannot use basis, requires face orientation")}
				face_key := Key{element.downward_connectivity[.D2][support.entity_index], support.entity_basis_index}
				if face_key in face_map {
					layout.mapping[element_idx][basis_index] = face_map[face_key]
				} else {
					layout.mapping[element_idx][basis_index] = next_global
					face_map[face_key] = next_global
					next_global += 1
				}
			case .D3:
				unreachable() // has to be interior
			case:
				unreachable()
			}
		}
	}

	layout.coeffs = make([]f64, next_global)

	return layout
}

create_sparsity_from_layout :: proc(mesh: Mesh, test, trial: DOF_Layout, allocator := context.allocator) -> (sp: la.Sparsity) {
	context.allocator = allocator

	rows := len(test.coeffs)
	sp.row_ptrs = make([]int, rows + 1, allocator)

	row_cols := make([]map[int]bool, rows, allocator)

	for i in 0 ..< rows {row_cols[i] = make(map[int]bool, allocator)}

	for element_id in 0 ..< len(mesh.elements) {
		test_dofs := test.mapping[element_id]
		trial_dofs := trial.mapping[element_id]
		for row in test_dofs {
			for col in trial_dofs {
				row_cols[row][col] = true
			}
		}
	}

	sp.row_ptrs[0] = 0
	for i in 0 ..< rows {sp.row_ptrs[i + 1] = sp.row_ptrs[i] + len(row_cols[i])}

	total_nnz := sp.row_ptrs[rows]
	sp.columns = make([]int, total_nnz, allocator)

	pos := 0
	for row in 0 ..< rows {
		for col in row_cols[row] {
			sp.columns[pos] = col
			pos += 1
		}
	}

	return sp
}


// returns the cannonical element-sub-entity-local index of a basis.
// so if an edge has 2 basis, say basis 3, and basis 4. basis 3 is the edges 0th basis and 4 is its 1st.
// if the edge is oriented wrong, when called with the support for basis 3 you will get a 1, rather than a 0.
orient_edge_dof :: proc(element: Mesh_Element, support: Basis_Support) -> int {
	assert(support.entity_dim == .D1)
	if element.edge_orientation[support.entity_index] {return support.entity_basis_index}
	return support.entity_basis_arity - 1 - support.entity_basis_index
}
