package small_strain

import "../../fem"

Bulk_Response :: struct {
    material: Material,
    source: Source,
    u_coeffs:     []f64,
	u_dot_coeffs: []f64,
	u_ddot_coeffs: []f64,
}

Boundary_Response :: struct {
    variational: Variational,
    displacement: []fem.Vec3,
}

Variational :: struct {
	t:     []fem.Vec3,
	u_ref: []fem.Vec3,
	k_s:   []f64,
}

Source :: struct {
	body: []fem.Vec3,
}

Material :: struct {
	stress:              []fem.Voigt6,
	constitutive_tensor: []fem.Voigt6x6,
	density:             []f64,
}

to_engineering_shear :: proc(eps: fem.Voigt6) -> fem.Voigt6 {
	return {eps[0], eps[1], eps[2], 2 * eps[3], 2 * eps[4], 2 * eps[5]}
}


assemble_bulk :: proc(
	ls:           fem.Local_System,
	U_handle:     fem.Var_Handle,
	space:        fem.Grad_Space(.Vector),
	quad:         fem.Mapped_Element,
	ts:           fem.Timestep,
	bulk:         Bulk_Response,
) {
	for qp in 0 ..< fem.space_points(space) {
	    C := bulk.material.constitutive_tensor[qp]
		for test in 0 ..< fem.space_arity(space) {
			residual: f64
			defer fem.local_system_rhs_add(ls, U_handle, test, residual)

            B_test := to_engineering_shear(fem.space_symmetric_gradient(space, quad, qp, test))
 			N_test := fem.space_value(space, qp, test)

 			// internal virtual work
 			residual -= fem.dot(B_test, bulk.material.stress[qp]) * fem.dV(quad, qp)

 			// body force
 			residual += fem.dot(N_test, bulk.source.body[qp]) * fem.dV(quad, qp)

			for trial in 0 ..< fem.space_arity(space) {
			    jacobian: f64
				defer fem.local_system_mat_add(ls, U_handle, test, U_handle, trial, jacobian)

				B_trial := to_engineering_shear(fem.space_symmetric_gradient(space, quad, qp, trial))
				N_trial := fem.space_value(space, qp, trial)

				// stiffness

				CB_trial := C * B_trial

				jacobian += fem.dot(B_test, CB_trial) * fem.dV(quad, qp)

				// mass
				M_ij := fem.dot(N_test, N_trial) * bulk.material.density[qp] * fem.dV(quad, qp)

				residual -= M_ij * bulk.u_ddot_coeffs[trial]
				jacobian += M_ij * ts.du_du_ddot[U_handle]
			}
		}
	}
}

assemble_boundary :: proc(
    ls:           fem.Local_System,
	U_handle:     fem.Var_Handle,
	space:        fem.Grad_Space(.Vector),
	quad:         fem.Mapped_Element,
	boundary:     Boundary_Response,
) {
    	for qp in 0 ..< fem.space_points(space) {
			for test in 0 ..< fem.space_arity(space) {
				residual: f64
				defer fem.local_system_rhs_add(ls, U_handle, test, residual)

				N_test := fem.space_value(space, qp, test)

				// traction
				residual += fem.dot(N_test, boundary.variational.t[qp]) * fem.dS(quad, qp)

				// spring
				residual += fem.dot(N_test, boundary.variational.k_s[qp] * (boundary.variational.u_ref[qp] - boundary.displacement[qp])) * fem.dS(quad, qp)

				for trial in 0 ..< fem.space_arity(space) {
					jacobian: f64
					defer fem.local_system_mat_add(ls, U_handle, test, U_handle, trial, jacobian)

					N_trial := fem.space_value(space, qp, trial)

					jacobian += boundary.variational.k_s[qp] * fem.dot(N_test, N_trial) * fem.dS(quad, qp)
				}
			}
		}
}


Problem_Data :: struct {
	ics:             map[fem.Section_ID]Initial_Condition_Int,
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][dynamic; 16]Source_Int,
	variational_bcs: map[fem.Boundary_ID][dynamic; 16]Variational_Int,
	fixed_bcs:  map[fem.Boundary_ID]Fixed_Int,
}

empty_problem_data :: proc(mesh: fem.Mesh, allocator := context.allocator) -> (d: Problem_Data) {
	context.allocator = allocator
    d = {
        ics = make(map[fem.Section_ID]Initial_Condition_Int),
	   materials = make(map[fem.Section_ID]Material_Int),
	   sources = make(map[fem.Section_ID][dynamic; 16]Source_Int),
	   variational_bcs = make(map[fem.Boundary_ID][dynamic; 16]Variational_Int),
	   fixed_bcs = make(map[fem.Boundary_ID]Fixed_Int),
	}
	// initialize those fixed capacity guys
	for _, id in mesh.boundary_names {d.variational_bcs[id] = {}}
	for _, id in mesh.section_names {d.sources[id] = {}}

    return d
}


Fixed_Int :: struct {
	procedure: fem.DOF_Functional_Vector_Proc,
	data:      rawptr,
}

Source_Int :: struct {
	procedure: Source_Proc,
	data:      rawptr,
}

Material_Int :: struct {
	procedure: Material_Proc,
	data:      rawptr,
}

Variational_Int :: struct {
	procedure: Variational_Proc,
	data:      rawptr,
}

Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Vector_Proc,
}


Variational_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Variational)

Source_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Source)

Material_Proc :: #type proc(
	mapped: fem.Mapped_Element,
	time: f64,
	current_strain: []fem.Voigt6,
	data: rawptr,
	out: Material,
)


empty_source :: proc(n_p: int) -> Source {return {make([]fem.Vec3, n_p)}}

empty_material :: proc(n_p: int) -> Material {
	return {make([]fem.Voigt6, n_p), make([]fem.Voigt6x6, n_p), make([]f64, n_p)}
}

empty_variational :: proc(n_p: int) -> Variational {
	return {make([]fem.Vec3, n_p), make([]fem.Vec3, n_p), make([]f64, n_p)}
}
