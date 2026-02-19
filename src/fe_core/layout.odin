// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package fem

import "core:mem"
import "core:slice"

import "../la"

MAX_FIELDS :: 16

Field_Handle :: distinct uint

Layout :: struct {
    fds: [MAX_FIELDS]Field_Descriptor,
    dof_layouts: [MAX_FIELDS]DOF_Layout,
    bound_maps: [MAX_FIELDS]Boundary_DOF_Map,
    field_sizes: [MAX_FIELDS]int,
    coupling: [MAX_FIELDS]bit_set[0..<MAX_FIELDS],
    // only valid for fields that were coupled.
    sparsities: [MAX_FIELDS][MAX_FIELDS]la.Sparsity,
    allocator: mem.Allocator,
    num_fields: int,
}

Field_Descriptor :: struct {
	family:             Basis_Family,
	order:              Order,
}

DOF_Layout :: struct {
	mapping:     [][]int,
	num_coeffs: int,
}

Boundary_DOF_Map :: map[Boundary_ID][]Boundary_DOF // for each bound id, which global dofs are "on it".

Boundary_DOF :: struct {
    local: int,
    element: Entity_ID,
}

// Constraint mask would be a parallel array beside the block vector for the problem.
Constraint_Mask :: struct {
    is_constrained: []bool,
    field_offsets: []int,
}

layout_create :: proc(allocator := context.allocator) -> Layout { return {allocator = allocator} }

// does not destroy any allocated vectors or matrices, but does invalidate the sparsity patterns used for the matrices.
layout_destroy :: proc(layout: ^Layout) {
    for field in 0..<layout.num_fields {
        for &local_map in layout.dof_layouts[field].mapping {delete(local_map)}
        delete(layout.dof_layouts[field].mapping)

        for coupled in layout.coupling[field] {
            delete(layout.sparsities[field][coupled].row_ptrs)
            delete(layout.sparsities[field][coupled].columns)
        }
    }
}

layout_add_field :: proc(layout: ^Layout, mesh: Mesh, fd: Field_Descriptor) -> Field_Handle {
    context.allocator = layout.allocator

    defer layout.num_fields += 1

    layout.dof_layouts[layout.num_fields] = create_dof_layout(mesh, fd)
    layout.fds[layout.num_fields] = fd
    layout.field_sizes[layout.num_fields] = layout.dof_layouts[layout.num_fields].num_coeffs

    layout.bound_maps[layout.num_fields] = create_boundary_dof_map(mesh, layout^, Field_Handle(layout.num_fields))

    return Field_Handle(layout.num_fields)
}

// coupling is always both-ways, no need to layout_couple(b, a)
layout_couple :: proc(layout: ^Layout, mesh: Mesh, a, b: Field_Handle) {
    context.allocator = layout.allocator

    layout.coupling[a] += {int(b)}
    layout.coupling[b] += {int(a)}

    layout.sparsities[a][b] = create_sparsity_from_layout(mesh, layout.dof_layouts[a], layout.dof_layouts[b])
    layout.sparsities[b][a] = create_sparsity_from_layout(mesh, layout.dof_layouts[b], layout.dof_layouts[a])
}


layout_make_vec :: proc(layout: ^Layout, allocator := context.allocator) -> la.Block_Vector {
    return la.block_vector_create(layout.field_sizes[:layout.num_fields], allocator)
}

layout_make_empty_constraint_mask :: proc(layout: ^Layout, allocator := context.allocator) -> Constraint_Mask {
    return constraint_mask_create(layout.field_sizes[:layout.num_fields], allocator)
}

// caller is responsible for deleting all but the sparsity patterns used.
layout_make_matrix :: proc(layout: ^Layout, allocator := context.allocator) -> la.Block_Sparse_Matrix {
    bs := la.block_sparse_create(layout.num_fields, layout.num_fields, allocator)
    for row in 0..<layout.num_fields {
        for col in layout.coupling[row] {
            sm := new(la.Sparse_Matrix, allocator)
            sm^ = la.sparse_matrix_from_sparsity(layout.sparsities[row][col], allocator)
            la.block_sparse_set_block(bs, row, col, sm)
        }
    }
    return bs
}

layout_local_problem :: proc(
    layout: ^Layout,
    element: Mesh_Element,
    allocator := context.allocator
) -> (la.Block_Vector, la.Block_Dense_Matrix) {

    infos: [MAX_FIELDS]Basis_Info
    for fd, i in layout.fds[:layout.num_fields] { infos[i] = basis_info_from_fd(fd, element)}
    return local_problem_from_basis(infos[:layout.num_fields], allocator)
}

// block matrix is in order of basis info's passed.
local_problem_from_basis :: proc(basis_infos: []Basis_Info, allocator := context.allocator) -> (la.Block_Vector, la.Block_Dense_Matrix) {
    sizes: [MAX_FIELDS]int

    assert(len(basis_infos) <= MAX_FIELDS, "Local problem is limited to MAX_FIELDS fields.")

    for info, i in basis_infos {sizes[i] = info.arity}

    m := la.block_dense_create(sizes[:len(basis_infos)], sizes[:len(basis_infos)], allocator)
    v := la.block_vector_create(sizes[:len(basis_infos)], allocator)

    return v, m
}

layout_scatter_local :: proc(
    layout: Layout,
    element_id: Entity_ID,
    A_g: la.Block_Sparse_Matrix,
    b_g: la.Block_Vector,
    A_l: la.Block_Dense_Matrix,
    b_l: la.Block_Vector,
    constraint_mask: Constraint_Mask,
) {
	for test_fh in 0..<layout.num_fields {
		local_r := la.block_vector_view(b_l, test_fh)
		global_r := la.block_vector_view(b_g, test_fh)
		test_constrained := constraint_mask_view(constraint_mask, test_fh)

		for trial_fh in 0..<layout.num_fields {
		    local_j := la.block_dense_view(A_l, test_fh, trial_fh)
		    global_j := la.block_sparse_view(A_g, test_fh, trial_fh)
		    trial_constrained := constraint_mask_view(constraint_mask, trial_fh)

			for local_test in 0..<local_j.rows {
				global_test := layout.dof_layouts[test_fh].mapping[element_id][local_test]

				if test_constrained[global_test] {continue}
				global_r[global_test] += local_r[local_test]

				if global_j == nil {continue}

				for local_trial in 0..<local_j.columns {
					J_ij := local_j.values[la.idx(local_j, local_test, local_trial)]
					global_trial := layout.dof_layouts[trial_fh].mapping[element_id][local_trial]
					if !trial_constrained[global_trial] {
						global_j.values[la.idx(global_j, global_test, global_trial)] += J_ij
					}
				}
			}
		}
	}
}

layout_finalize_constraints :: proc(layout: Layout, A_g: la.Block_Sparse_Matrix, b_g: la.Block_Vector, mask: Constraint_Mask) {
	for fh in 0..<layout.num_fields {
		for is_constrained, i in constraint_mask_view(mask, fh) {
			if !is_constrained{continue}

			b := la.block_vector_view(b_g, fh)

			for fh_b in 0..<layout.num_fields {
				if block := la.block_sparse_view(A_g, fh, fh_b); block != nil {
					la.sparse_mat_zero_row(block^, i)
				}
			}
			if diag := la.block_sparse_view(A_g, fh, fh); diag != nil {
				diag^.values[la.idx(diag^, i, i)] = 1
			}

			b[i] = 0
		}
	}
}



// not generic on basis type, cannot be used to get actual basis values or gradients just metadata.
basis_info_from_fd :: proc(fd: Field_Descriptor, element: Element) -> Basis_Info {
	switch fd.family {
	case .LS:
		return basis_create(Basis_LS, element, fd.order).info
	case .LV:
		return basis_create(Basis_LV, element, fd.order).info
	case:
		unreachable()
	}
}


create_dof_layout :: proc(
    mesh: Mesh,
    fd: Field_Descriptor,
    allocator := context.allocator,
) -> (layout: DOF_Layout) {
    orient_edge_dof :: proc(element: Mesh_Element, support: Basis_Support) -> int {
        assert(support.entity_dim == .D1)
        if element.edge_orientation[support.entity_index] { return support.entity_basis_index }
        return support.entity_basis_arity - 1 - support.entity_basis_index
    }

    assign_shared_dof :: proc(m: ^map[Key]int, key: Key, next_global: ^int) -> int {
        if existing, ok := m[key]; ok { return existing }
        idx := next_global^
        m[key] = idx
        next_global^ += 1
        return idx
    }

    context.allocator = allocator

    Key :: struct { entity_id: Entity_ID, basis_index: int }
    vertex_map := make(map[Key]int)
    edge_map   := make(map[Key]int)
    face_map   := make(map[Key]int)

    defer {
        delete(vertex_map)
        delete(edge_map)
        delete(face_map)
    }

    layout.mapping = make([][]int, len(mesh.elements))
    next_global: int

    for element, element_idx in mesh.elements {
        b := basis_info_from_fd(fd, element)
        layout.mapping[element_idx] = make([]int, b.arity)

        for basis_index in 0 ..< b.arity {
            support := basis_support(b, basis_index)
            conn  := element.downward_connectivity

            global_dof: ^int = &layout.mapping[element_idx][basis_index]
            switch support.entity_dim {
            case element_dim(element.type):
                global_dof^ = next_global
                next_global += 1
            case .D0:
                global_dof^ = assign_shared_dof(&vertex_map, {conn[.D0][support.entity_index], support.entity_basis_index}, &next_global)
            case .D1:
                global_dof^ = assign_shared_dof(&edge_map, {conn[.D1][support.entity_index], orient_edge_dof(element, support)}, &next_global)
            case .D2:
                assert(support.entity_basis_arity <= 1, "Basis with multiple face supported dofs are currently unimplemented.")
                global_dof^ = assign_shared_dof(&face_map, {conn[.D2][support.entity_index], support.entity_basis_index}, &next_global)
            case .D3: unreachable()
            case: unreachable()
            }
        }
    }

    layout.num_coeffs = next_global
    return
}

create_boundary_dof_map :: proc(mesh: Mesh, layout: Layout, fh: Field_Handle, allocator := context.allocator) -> (bdm: Boundary_DOF_Map) {
    context.allocator = allocator

    bdm = make(map[Boundary_ID][]Boundary_DOF)

    tmp := make(map[Boundary_ID][dynamic]Boundary_DOF)
    defer { for _, &arr in tmp { delete(arr) }; delete(tmp) }

    fd := layout.fds[fh]
    for element in mesh.elements {
        basis := basis_info_from_fd(fd, element)
        for bnd_facet, facet_index in element.boundaries {
            id, has := bnd_facet.?
            if !has { continue }
            if id not_in tmp { tmp[id] = make([dynamic]Boundary_DOF) }
            for dof in 0..<basis.arity {
                support := basis_support(basis, dof)
                if !element_entity_on_facet(element.type, facet_index, support.entity_dim, support.entity_index) { continue }
                bnd_dof: Boundary_DOF = {
                    local = dof,
                    element = element.id
                }
                append(&tmp[id], bnd_dof)
            }
        }
    }

    for id, &arr in tmp { bdm[id] = slice.clone(arr[:]) }

    return bdm

}

create_sparsity_from_layout :: proc(
    mesh: Mesh,
    test, trial: DOF_Layout,
    allocator := context.allocator,
) -> (sp: la.Sparsity) {
    context.allocator = allocator
    rows := test.num_coeffs
    sp.row_ptrs = make([]int, rows + 1)

    row_cols := make([]map[int]bool, rows)
    for &m in row_cols { m = make(map[int]bool) }
    defer { for &m in row_cols {delete(m)}; delete(row_cols) }

    for element_id in 0 ..< len(mesh.elements) {
        for row in test.mapping[element_id] {
            for col in trial.mapping[element_id] { row_cols[row][col] = true }
        }
    }

    for i in 0 ..< rows { sp.row_ptrs[i + 1] = sp.row_ptrs[i] + len(row_cols[i]) }

    sp.columns = make([]int, sp.row_ptrs[rows])
    pos := 0
    for row in 0 ..< rows {
        for col in row_cols[row] {
            sp.columns[pos] = col
            pos += 1
        }
    }

    return sp
}

// finds the absolute index in a system vector a given local dof, of a given field
layout_global_pos :: proc(layout: Layout, field: Field_Handle, element_id: Entity_ID, local_dof: int) -> int {
    return layout.dof_layouts[field].mapping[element_id][local_dof]
}

constraint_mask_view :: proc(cm: Constraint_Mask, field: int) -> []bool {
     return cm.is_constrained[cm.field_offsets[field]:cm.field_offsets[field + 1]]
}

constraint_mask_create :: proc(field_sizes: []int, allocator := context.allocator) -> Constraint_Mask {
     n := len(field_sizes)

    offsets := make([]int, n + 1, allocator)
    for i in 0..<n { offsets[i + 1] = offsets[i] + field_sizes[i] }

    return {
        make([]bool, offsets[n], allocator),
        offsets,
    }
}