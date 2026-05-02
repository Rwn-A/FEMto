package fem

import "core:mem"
import "core:slice"

import "base:intrinsics"

import "./infra"

MAX_VARS :: 8

// Handles are contiguous u8s starting from 0; field_offsets has num_vars+1 slots.
Var_Handle :: distinct u8

// Flat index into a global solution vector or constraint mask.
Global_DOF :: distinct int

// Index into the per-element local DOF array for one variable.
Local_DOF :: distinct int

// Index into Sparse_Matrix.values.
Sparse_Slot :: distinct int

System_Description :: struct {
	variables: [MAX_VARS]Variable_Description,
	num_vars:  Var_Handle,
}

System :: struct {
	using desc:    System_Description,
	dof_layouts:   [MAX_VARS]DOF_Layout,
	sparsity:      Sparsity,
	field_offsets: []int, // length num_vars+1
	allocator:     mem.Allocator,
}

Variable_Description :: struct {
	bd:         Basis_Descriptor,
	components: int, // COMPONENTS PER BASIS DOF, NOT NESSECARILY VALUE COMPONENTS (RT BASIS VECTOR FIELD BUT 1 COMPONENT HERE.)
	coupled_to: bit_set[0 ..< MAX_VARS],
}

DOF_Layout :: struct {
	local_to_global: [][]int,
	boundary_dofs:   map[Boundary_ID][]Boundary_DOF,
	num_coeffs:      int,
}

Boundary_DOF :: struct {
	element_id: Entity_ID,
	local_dofs: []Local_DOF,
}

Vector_Orphan :: struct {
	slot:  Global_DOF,
	value: f64,
}

Matrix_Orphan :: struct {
	slot:  Sparse_Slot,
	value: f64,
}

Thread_Partition :: struct {
	owned: [MAX_VARS]infra.Range,
	lhs:   [dynamic]Matrix_Orphan,
	rhs:   [dynamic]Vector_Orphan,
}

Thread_Partitions :: struct {
	partitions: []^Thread_Partition,
	allocator:  mem.Allocator,
}

// ----------------------
// Local assembly context
// ----------------------

Local_System :: struct {
	mat:     Dense_Matrix,
	rhs:     Vector,
	offsets: []int, // length num_vars+1, same layout as System.field_offsets jsut for local vars
}

local_system_destroy :: proc(ls: ^Local_System, allocator := context.allocator) {
	dense_matrix_destroy(ls.mat, allocator)
	delete(ls.rhs, allocator)
	delete(ls.offsets, allocator)
}

@(private)
local_system_offset :: #force_inline proc(ls: Local_System, var: Var_Handle, dof: Local_DOF) -> int {
	return ls.offsets[int(var)] + int(dof)
}

local_system_mat_add :: #force_inline proc(
	ls: Local_System,
	test_var: Var_Handle,
	#any_int test_dof: Local_DOF,
	trial_var: Var_Handle,
	#any_int trial_dof: Local_DOF,
	value: f64,
) {
	r := local_system_offset(ls, test_var, test_dof)
	c := local_system_offset(ls, trial_var, trial_dof)
	ls.mat.values[idx(ls.mat, r, c)] += value
}

local_system_rhs_add :: #force_inline proc(
	ls: Local_System,
	test_var: Var_Handle,
	#any_int test_dof: Local_DOF,
	value: f64,
) {
	ls.rhs[local_system_offset(ls, test_var, test_dof)] += value
}

local_system_mat_get :: #force_inline proc(
	ls: Local_System,
	test_var: Var_Handle,
	test_dof: Local_DOF,
	trial_var: Var_Handle,
	trial_dof: Local_DOF,
) -> f64 {
	r := local_system_offset(ls, test_var, test_dof)
	c := local_system_offset(ls, trial_var, trial_dof)
	return ls.mat.values[idx(ls.mat, r, c)]
}

local_system_rhs_get :: #force_inline proc(ls: Local_System, test_var: Var_Handle, test_dof: Local_DOF) -> f64 {
	return ls.rhs[local_system_offset(ls, test_var, test_dof)]
}

// ---------------------------
// System field-offset helpers
// ---------------------------

// Global DOF range [begin, end) for a variable in the flat solution vector.
system_var_range :: #force_inline proc(sys: System, var: Var_Handle) -> (begin, end: Global_DOF) {
	return Global_DOF(sys.field_offsets[int(var)]), Global_DOF(sys.field_offsets[int(var) + 1])
}

// Total number of global DOFs across all variables.
system_total_dofs :: #force_inline proc(sys: System) -> int {
	return sys.field_offsets[int(sys.num_vars)]
}

// Resolve a variable-local element DOF to a flat Global_DOF.
system_global_dof :: #force_inline proc(
	sys: System,
	var: Var_Handle,
	element_id: Entity_ID,
	local_dof: Local_DOF,
) -> Global_DOF {
	begin, _ := system_var_range(sys, var)
	return begin + Global_DOF(sys.dof_layouts[var].local_to_global[element_id][local_dof])
}

// ----------------------------
// Description building helpers
// -----------------------------

description_add_variable :: proc(desc: ^System_Description, var: Variable_Description) -> Var_Handle {
	assert(desc.num_vars < MAX_VARS, "System cannot define more variables; increase MAX_VARS")
	handle := desc.num_vars
	desc.variables[handle] = var
	desc.num_vars += 1
	return handle
}

description_couple :: proc(desc: ^System_Description, var_a: Var_Handle, var_b: Var_Handle) {
	desc.variables[var_a].coupled_to += {int(var_b)}
}

// ----------------
// System lifecycle
// ----------------

system_from_description :: proc(mesh: Mesh, desc: System_Description, allocator: mem.Allocator) -> (sys: System) {
	sys.desc = desc
	sys.allocator = allocator

	sys.field_offsets = make([]int, int(sys.num_vars) + 1, allocator)

	total_size: int
	for var in Var_Handle(0) ..< sys.num_vars {
		sys.dof_layouts[var] = dof_layout(mesh, sys.variables[var].bd, sys.variables[var].components, allocator)
		sys.field_offsets[int(var)] = total_size
		total_size += sys.dof_layouts[var].num_coeffs
	}
	sys.field_offsets[int(sys.num_vars)] = total_size

	sys.sparsity = system_sparsity(mesh, sys, allocator)

	return sys
}

system_destroy :: proc(sys: ^System) {
	context.allocator = sys.allocator

	delete(sys.sparsity.row_ptrs)
	delete(sys.sparsity.columns)
	delete(sys.field_offsets)

	for var in Var_Handle(0) ..< sys.num_vars {
		layout := &sys.dof_layouts[var]
		for &arr in layout.local_to_global {delete(arr)}
		delete(layout.local_to_global)

		for _, &v in layout.boundary_dofs {
			for bnd_dof in v {delete(bnd_dof.local_dofs)}
			delete(v)
		}
		delete(layout.boundary_dofs)
	}
}

// -------------------------
// System allocation helpers
// -------------------------

system_vector :: proc(sys: System, allocator := context.allocator) -> Vector {
	return make(Vector, system_total_dofs(sys), allocator)
}

system_constraint_mask :: proc(sys: System, allocator := context.allocator) -> []bool {
	return make([]bool, system_total_dofs(sys), allocator)
}

system_matrix :: proc(sys: System, allocator := context.allocator) -> Sparse_Matrix {
	return sparse_matrix_from_sparsity(sys.sparsity, allocator)
}

// ------------
// Field access
// ------------

// Returns the sub-slice of data corresponding to var's global DOFs.
system_var_slice :: proc(sys: System, var: Var_Handle, data: Vector) -> Vector {
	b, e := system_var_range(sys, var)
	return data[b:e]
}

// Gathers the coefficients for all local DOFs of var on element_id into a
// freshly allocated slice in local DOF order. Intended for use in weak form
// loops: coeffs[j] aligns with the same j used to index basis values.
system_gather_var_coeffs :: proc(
	sys: System,
	var: Var_Handle,
	element_id: Entity_ID,
	data: Vector,
	allocator := context.allocator,
) -> []f64 {
	coeffs := make([]f64, len(sys.dof_layouts[var].local_to_global[element_id]), allocator)
	for &coeff, i in coeffs {
		coeff = data[system_global_dof(sys, var, element_id, Local_DOF(i))]
	}
	return coeffs
}

evaluate_scalar_var :: proc(space: Grad_Space(.Scalar), coeffs: []f64, allocator := context.allocator) -> []f64 {
	return evaluate_var_generic(f64, space, coeffs, allocator)
}

evaluate_vector_var :: proc(
	space: $T,
	coeffs: []f64,
	allocator := context.allocator,
) -> []Vec3 where T ==
	Grad_Space(.Vector) { 	// where clause here as Div_Space and Curl_Space if ever implemented would also be valid to call here.
	return evaluate_var_generic(Vec3, space, coeffs, allocator)
}

evaluate_var :: proc {
	evaluate_scalar_var,
	evaluate_vector_var,
}

evaluate_var_gradient :: proc {
	evaluate_scalar_var_gradient,
}

@(private)
evaluate_var_generic :: proc(
	$OUTPUT_TYPE: typeid,
	space: $T,
	coeffs: []f64,
	allocator := context.allocator,
) -> []OUTPUT_TYPE {
	out := make([]OUTPUT_TYPE, space_points(space), allocator)
	for p in 0 ..< space_points(space) {
		for dof in 0 ..< space_arity(space) {
			out[p] += coeffs[dof] * space_value(space, p, dof)
		}
	}
	return out
}

@(private)
evaluate_var_gradient_generic :: proc(
	$OUTPUT_TYPE: typeid,
	space: $T,
	mapping: Mapped_Element,
	coeffs: []f64,
	allocator := context.allocator,
) -> []OUTPUT_TYPE {
	out := make([]OUTPUT_TYPE, space_points(space), allocator)
	for p in 0 ..< space_points(space) {
		for dof in 0 ..< space_arity(space) {
			out[p] += coeffs[dof] * space_gradient(space, mapping, p, dof)
		}
	}
	return out
}

evaluate_scalar_var_gradient :: proc(
	space: Grad_Space(.Scalar),
	mapping: Mapped_Element,
	coeffs: []f64,
	allocator := context.allocator,
) -> []Vec3 {
	return evaluate_var_gradient_generic(Vec3, space, mapping, coeffs, allocator)
}

system_var_bd :: proc(sys: System, var: Var_Handle) -> Basis_Descriptor {
	return sys.variables[var].bd
}

// ------------------------------------------------------------------
// Boundary constraint helpers (used for strong dirichlet conditions)
// ------------------------------------------------------------------


Boundary_DOF_Visit :: struct {
	gdof:       Global_DOF,
	element_id: Entity_ID,
	basis_dof:  int,
	component:  int,
}

// Marks all DOFs for (var, id) as constrained. Call once per variable per
// boundary before system_finalize_constraints.
system_mark_boundary :: proc(sys: System, var: Var_Handle, id: Boundary_ID, mask: []bool) {
	bnd_dofs, ok := sys.dof_layouts[var].boundary_dofs[id]
	if !ok {return}

	for bnd_dof in bnd_dofs {
		for local_dof in bnd_dof.local_dofs {
			mask[system_global_dof(sys, var, bnd_dof.element_id, local_dof)] = true
		}
	}
}


Boundary_DOF_Iterator :: struct {
	sys:       System,
	var:       Var_Handle,
	id:        Boundary_ID,
	bnd_idx:   int,
	local_idx: int,
}

boundary_dof_iter_create :: proc(sys: System, var: Var_Handle, id: Boundary_ID) -> Boundary_DOF_Iterator {
	return {sys = sys, var = var, id = id}
}

boundary_dof_iter_next :: proc(it: ^Boundary_DOF_Iterator) -> (visit: Boundary_DOF_Visit, ok: bool) {
	bnd_dofs := it.sys.dof_layouts[it.var].boundary_dofs[it.id] or_return

	for it.bnd_idx < len(bnd_dofs) {
		bnd_dof := bnd_dofs[it.bnd_idx]
		if it.local_idx < len(bnd_dof.local_dofs) {
			local_dof := bnd_dof.local_dofs[it.local_idx]
			components := it.sys.variables[it.var].components

			visit = Boundary_DOF_Visit {
				gdof       = system_global_dof(it.sys, it.var, bnd_dof.element_id, local_dof),
				element_id = bnd_dof.element_id,
				basis_dof  = int(local_dof) / components,
				component  = int(local_dof) % components,
			}

			it.local_idx += 1

			return visit, true
		}
		it.bnd_idx += 1
		it.local_idx = 0
	}
	return {}, false
}

// ---------------------------------------------------------------------------
// Assembly loop calls
// ---------------------------------------------------------------------------

// Element local matrix and source vector
system_local_problem :: proc(sys: System, element_id: Entity_ID, allocator := context.allocator) -> Local_System {
	offsets := make([]int, int(sys.num_vars) + 1, allocator)
	current_size: int

	for var in 0 ..< sys.num_vars {
		offsets[int(var)] = current_size
		current_size += len(sys.dof_layouts[var].local_to_global[element_id])
	}
	offsets[int(sys.num_vars)] = current_size

	return Local_System {
		mat = dense_matrix_create(current_size, current_size, allocator),
		rhs = make(Vector, current_size, allocator),
		offsets = offsets,
	}
}


@(private = "file")
scatter_columns :: #force_inline proc(
	sys: System,
	element: Entity_ID,
	ls: Local_System,
	A_g: Sparse_Matrix,
	constraint_mask: []bool,
	test_var: Var_Handle,
	i: int,
	global_row: Global_DOF,
) {
	for trial_var in Var_Handle(0) ..< sys.num_vars {
		if int(trial_var) not_in sys.variables[test_var].coupled_to {continue}
		trial_dofs := sys.dof_layouts[trial_var].local_to_global[element]
		for j in 0 ..< len(trial_dofs) {
			global_col := system_global_dof(sys, trial_var, element, Local_DOF(j))
			if constraint_mask[global_col] {continue}
			mat_val := local_system_mat_get(ls, test_var, Local_DOF(i), trial_var, Local_DOF(j))
			mat_slot := Sparse_Slot(idx(A_g, int(global_row), int(global_col)))
			A_g.values[mat_slot] += mat_val
		}
	}
}

// Accumulates the local system into the global system.
system_scatter :: proc(
	sys: System,
	element: Entity_ID,
	A_g: Sparse_Matrix,
	b_g: Vector,
	constraint_mask: []bool,
	ls: Local_System,
) {
	for test_var in Var_Handle(0) ..< sys.num_vars {
		test_dofs := sys.dof_layouts[test_var].local_to_global[element]
		for i in 0 ..< len(test_dofs) {
			global_row := system_global_dof(sys, test_var, element, Local_DOF(i))
			if constraint_mask[global_row] {continue}
			b_g[global_row] += local_system_rhs_get(ls, test_var, Local_DOF(i))
			scatter_columns(sys, element, ls, A_g, constraint_mask, test_var, i, global_row)
		}
	}
}

// Threaded scatter, requires calling `system_flush_orphans` after assembly and setup of partitions.
system_thread_scatter :: proc(
	sys: System,
	element: Entity_ID,
	A_g: Sparse_Matrix,
	b_g: Vector,
	constraint_mask: []bool,
	ls: Local_System,
	partition: ^Thread_Partition,
) {
	for test_var in Var_Handle(0) ..< sys.num_vars {
		test_dofs := sys.dof_layouts[test_var].local_to_global[element]
		for i in 0 ..< len(test_dofs) {
			global_row := system_global_dof(sys, test_var, element, Local_DOF(i))
			if constraint_mask[global_row] {continue}
			rhs_val := local_system_rhs_get(ls, test_var, Local_DOF(i))
			if infra.in_range(test_dofs[i], partition.owned[test_var]) {
				b_g[global_row] += rhs_val
				scatter_columns(sys, element, ls, A_g, constraint_mask, test_var, i, global_row)
			} else {
				append(&partition.rhs, Vector_Orphan{global_row, rhs_val})
				for trial_var in Var_Handle(0) ..< sys.num_vars {
					if int(trial_var) not_in sys.variables[test_var].coupled_to {continue}
					trial_dofs := sys.dof_layouts[trial_var].local_to_global[element]
					for j in 0 ..< len(trial_dofs) {
						global_col := system_global_dof(sys, trial_var, element, Local_DOF(j))
						if constraint_mask[global_col] {continue}
						mat_val := local_system_mat_get(ls, test_var, Local_DOF(i), trial_var, Local_DOF(j))
						mat_slot := Sparse_Slot(idx(A_g, int(global_row), int(global_col)))
						append(&partition.lhs, Matrix_Orphan{mat_slot, mat_val})
					}
				}
			}
		}
	}
}

// Call after assembly but before finalizing constraints.
system_flush_orphans :: proc(partitions: Thread_Partitions, A_g: Sparse_Matrix, b_g: Vector) {
	for p in partitions.partitions {
		for o in p.rhs {b_g[o.slot] += o.value}
		for o in p.lhs {A_g.values[o.slot] += o.value}
		clear(&p.rhs)
		clear(&p.lhs)
	}
}

// Must be called after assembly (and after flush_orphans if threaded),
// before solving. Sets constrained rows to identity in A and zero in b.
system_finalize_constraints :: proc(sys: System, A_g: Sparse_Matrix, b_g: Vector, constraint_mask: []bool) {
	for var in Var_Handle(0) ..< sys.num_vars {
		begin, end := system_var_range(sys, var)

		for gdof in begin ..< end {
			if !constraint_mask[gdof] {continue}

			// Residual fomrulation zero residual for constrained rows.
			b_g[gdof] = 0

			row_start := A_g.row_ptrs[gdof]
			row_end := A_g.row_ptrs[gdof + 1]
			for k in row_start ..< row_end {A_g.values[k] = 0}

			// Diagonal = 1 so the system remains non-singular.
			A_g.values[idx(A_g, int(gdof), int(gdof))] = 1
		}
	}
}

// ---------------------------------------------------------------------------
// Thread partition lifecycle
// ---------------------------------------------------------------------------

// Allocator must be thread-safe.
system_create_thread_partitions :: proc(
	sys: System,
	prt: ^infra.Parallel_Runtime,
	allocator: mem.Allocator,
) -> Thread_Partitions {
	parts := make([]^Thread_Partition, prt.total_threads, allocator)
	for t in 0 ..< prt.total_threads {
		parts[t] = new(Thread_Partition, allocator)
		for var in Var_Handle(0) ..< sys.num_vars {
			parts[t].owned[var] = infra.parallel_runtime_partition(prt, {0, sys.dof_layouts[var].num_coeffs}, t)
		}
		parts[t].lhs = make([dynamic]Matrix_Orphan, 0, 1024, allocator)
		parts[t].rhs = make([dynamic]Vector_Orphan, 0, 1024, allocator)
	}
	return Thread_Partitions{partitions = parts, allocator = allocator}
}

system_destroy_thread_partitions :: proc(tp: ^Thread_Partitions) {
	context.allocator = tp.allocator
	for p in tp.partitions {
		delete(p.lhs)
		delete(p.rhs)
		free(p)
	}
	delete(tp.partitions)
}

// ---------------------------------------------------------------------------
// Internally used by system
// ---------------------------------------------------------------------------

system_sparsity :: proc(mesh: Mesh, sys: System, allocator := context.allocator) -> (sp: Sparsity) {
	context.allocator = allocator

	rows := system_total_dofs(sys)
	sp.row_ptrs = make([]i32, rows + 1)

	row_cols := make([]map[int]bool, rows)
	for &m in row_cols {m = make(map[int]bool)}
	defer {for &m in row_cols {delete(m)}; delete(row_cols)}

	for element in 0 ..< len(mesh.elements) {
		for ti in Var_Handle(0) ..< sys.num_vars {
			row_off := sys.field_offsets[int(ti)]
			test_dofs := sys.dof_layouts[ti].local_to_global[element]

			// Always insert diagonal — required for strong constraint rows.
			for row in test_dofs {
				row_cols[row_off + row][row_off + row] = true
			}

			for tj in sys.desc.variables[ti].coupled_to {
				col_off := sys.field_offsets[tj]
				trial_dofs := sys.dof_layouts[tj].local_to_global[element]

				for row in test_dofs {
					for col in trial_dofs {
						row_cols[row_off + row][col_off + col] = true
					}
				}
			}
		}
	}

	for i in 0 ..< rows {
		sp.row_ptrs[i + 1] = sp.row_ptrs[i] + i32(len(row_cols[i]))
	}

	sp.columns = make([]i32, sp.row_ptrs[rows])

	for row in 0 ..< rows {
		start := sp.row_ptrs[row]
		end := sp.row_ptrs[row + 1]
		i := start
		for col in row_cols[row] {
			sp.columns[i] = i32(col)
			i += 1
		}
		slice.sort(sp.columns[start:end])
	}

	return sp
}

dof_layout :: proc(
	mesh: Mesh,
	bd: Basis_Descriptor,
	components: int,
	allocator := context.allocator,
) -> (
	layout: DOF_Layout,
) {
	orient_edge_dof :: proc(element: Mesh_Element, support: Basis_Support, components: int, cmpnt: int) -> int {
		assert(support.entity_dim == .D1)
		if element.edge_orientation[support.entity_index] {
			return support.entity_basis_index * components + cmpnt
		}
		reversed_scalar := support.entity_basis_arity - 1 - support.entity_basis_index
		return reversed_scalar * components + cmpnt
	}

	assign_shared_dof :: proc(m: ^map[Key]int, key: Key, next_global: ^int) -> int {
		if existing, ok := m[key]; ok {return existing}
		idx := next_global^
		m[key] = idx
		next_global^ += 1
		return idx
	}

	context.allocator = allocator

	Key :: struct {
		entity_id:   Entity_ID,
		basis_index: int,
	}
	vertex_map := make(map[Key]int)
	edge_map := make(map[Key]int)
	face_map := make(map[Key]int)
	defer {delete(vertex_map); delete(edge_map); delete(face_map)}

	layout.local_to_global = make([][]int, len(mesh.elements))
	layout.boundary_dofs = make(map[Boundary_ID][]Boundary_DOF)
	next_global: int

	Bnd_Key :: struct {
		element_idx: int,
		id:          Boundary_ID,
	}
	bnd_dof_builder := make(map[Bnd_Key][dynamic]Local_DOF)
	defer {for _, &arr in bnd_dof_builder {delete(arr)}; delete(bnd_dof_builder)}

	for element, element_idx in mesh.elements {
		num_basis := basis_count(element.type, bd)
		num_dofs := num_basis * components
		layout.local_to_global[element_idx] = make([]int, num_dofs)
		conn := element.downward_connectivity

		for support, support_idx in basis_support(element.type, bd) {
			for cmpnt in 0 ..< components {
				local_dof := support_idx * components + cmpnt

				for bnd_facet in element.boundary_facets {
					if support_on_facet(element, support, bnd_facet) {
						id := element.boundary_ids[bnd_facet]
						key := Bnd_Key{element_idx, id}
						if key not_in bnd_dof_builder {bnd_dof_builder[key] = make([dynamic]Local_DOF)}
						append(&bnd_dof_builder[key], Local_DOF(local_dof))
					}
				}

				global_dof := &layout.local_to_global[element_idx][local_dof]
				switch support.entity_dim {
				case element_dim(element.type):
					global_dof^ = next_global
					next_global += 1
				case .D0:
					global_dof^ = assign_shared_dof(
						&vertex_map,
						{conn[.D0][support.entity_index], support.entity_basis_index * components + cmpnt},
						&next_global,
					)
				case .D1:
					global_dof^ = assign_shared_dof(
						&edge_map,
						{conn[.D1][support.entity_index], orient_edge_dof(element, support, components, cmpnt)},
						&next_global,
					)
				case .D2:
					assert(
						support.entity_basis_arity <= 1,
						"Basis with multiple face-supported dofs is currently unimplemented.",
					)
					global_dof^ = assign_shared_dof(
						&face_map,
						{conn[.D2][support.entity_index], support.entity_basis_index * components + cmpnt},
						&next_global,
					)
				case .D3:
					unreachable()
				case:
					unreachable()
				}
			}
		}
	}

	Id_Elem_Map :: map[Boundary_ID][dynamic]Boundary_DOF
	id_map := make(Id_Elem_Map)
	defer {delete(id_map)}

	for key, &arr in bnd_dof_builder {
		if key.id not_in id_map {id_map[key.id] = make([dynamic]Boundary_DOF)}
		append(
			&id_map[key.id],
			Boundary_DOF{element_id = Entity_ID(key.element_idx), local_dofs = slice.clone(arr[:])},
		)
	}
	for id, &arr in id_map {
		layout.boundary_dofs[id] = slice.clone(arr[:])
		delete(arr)
	}

	layout.num_coeffs = next_global
	return
}
