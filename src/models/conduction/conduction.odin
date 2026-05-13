package conduction

import "../../fem"

Bulk_Response :: struct {
    material:     Material,
	source:       Source,
	u_coeffs:     []f64,
	u_dot_coeffs: []f64,
	grad_T: []fem.Vec3,
}

Boundary_Response :: struct {
    bc:       Variational_BC,
	u_coeffs: []f64,
}

Source :: struct {
	Q:   []f64,
	d_Q: []f64,	
}

Material :: struct {
	k, rho, cp:       []f64,
	d_k, d_rho, d_cp: []f64,
}

Variational_BC :: struct {
	q, h, T_amb: []f64,
	d_q, d_h:    []f64,
}

// Assembles the volume integral (bulk) terms.
assemble_bulk :: proc(
	ls:           fem.Local_System,
	T_handle:     fem.Var_Handle,
	space:        fem.Grad_Space(.Scalar),
	quad:         fem.Mapped_Element,
	ts:           fem.Timestep,
	bulk:         Bulk_Response,
) {
	for qp in 0 ..< fem.space_points(space) {
		for test in 0 ..< fem.space_arity(space) {
			residual: f64
			defer fem.local_system_rhs_add(ls, T_handle, test, residual)

			grad_test := fem.space_gradient(space, quad, qp, test)
			val_test  := fem.space_value(space, qp, test)

			// sources (linear part)
			residual += bulk.source.Q[qp] * val_test * fem.dV(quad, qp)

			for trial in 0 ..< fem.space_arity(space) {
				jacobian: f64
				defer fem.local_system_mat_add(ls, T_handle, test, T_handle, trial, jacobian)

				grad_trial := fem.space_gradient(space, quad, qp, trial)
				val_trial  := fem.space_value(space, qp, trial)

				// stiffness
				K_ij       := fem.dot(grad_test, grad_trial) * bulk.material.k[qp] * fem.dV(quad, qp)
				K_ij_deriv := fem.dot(bulk.grad_T[qp], grad_test) * bulk.material.d_k[qp] * val_trial * fem.dV(quad, qp)

				residual -= K_ij * bulk.u_coeffs[trial]
				jacobian += K_ij + K_ij_deriv

				// sources (nonlinear)
				jacobian -= bulk.source.d_Q[qp] * val_test * val_trial * fem.dV(quad, qp)

				// mass (timestepping)
				M_ij := val_trial * val_test * bulk.material.rho[qp] * bulk.material.cp[qp] * fem.dV(quad, qp)
				M_ij_deriv :=
					(val_trial *
						val_test *
						(bulk.material.d_rho[qp] * bulk.material.cp[qp] +
							bulk.material.rho[qp] * bulk.material.d_cp[qp])) *
					bulk.u_dot_coeffs[trial] *
					fem.dV(quad, qp)

				residual -= M_ij * bulk.u_dot_coeffs[trial]
				jacobian += M_ij * ts.du_du_dot[T_handle] - M_ij_deriv
			}
		}
	}
}

// Assembles the boundary terms.
assemble_boundary :: proc(
	ls:       fem.Local_System,
	T_handle: fem.Var_Handle,
	space:    fem.Grad_Space(.Scalar),
	quad:     fem.Mapped_Element,
	boundary: Boundary_Response,
) {
	for qp in 0 ..< fem.space_points(space) {
		for test in 0 ..< fem.space_arity(space) {
			residual: f64
			defer fem.local_system_rhs_add(ls, T_handle, test, residual)

			val_test := fem.space_value(space, qp, test)

			residual +=
				(boundary.bc.q[qp] +
					boundary.bc.h[qp] * boundary.bc.T_amb[qp] -
					boundary.bc.h[qp] * boundary.u_coeffs[qp]) *
				val_test *
				fem.dS(quad, qp)

			for trial in 0 ..< fem.space_arity(space) {
				jacobian: f64
				defer fem.local_system_mat_add(ls, T_handle, test, T_handle, trial, jacobian)

				val_trial := fem.space_value(space, qp, trial)

				jacobian += boundary.bc.h[qp] * val_trial * val_test * fem.dS(quad, qp)
				// nonlinear h
				jacobian +=
					boundary.bc.d_h[qp] *
					(boundary.bc.T_amb[qp] - boundary.u_coeffs[qp]) *
					val_trial * val_test *
					fem.dS(quad, qp)
				// nonlinear q
				jacobian -= boundary.bc.d_q[qp] * val_trial * val_test * fem.dS(quad, qp)
			}
		}
	}
}


// Default interfaces used for conduction as the primary model.
// In a multi-physics model some of these may be superceeded by an interface defined by the coupled model.

Problem_Data :: struct {
	ics:             map[fem.Section_ID]Initial_Condition_Int,
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][dynamic; 16]Source_Int,
	variational_bcs: map[fem.Boundary_ID][dynamic; 16]Variational_Int,
	isothermal_bcs:  map[fem.Boundary_ID]Isothermal_Int,
}

empty_problem_data :: proc(mesh: fem.Mesh, allocator := context.allocator) -> (d: Problem_Data) {
	context.allocator = allocator
    d = {
        ics = make(map[fem.Section_ID]Initial_Condition_Int),
	   materials = make(map[fem.Section_ID]Material_Int),
	   sources = make(map[fem.Section_ID][dynamic; 16]Source_Int),
	   variational_bcs = make(map[fem.Boundary_ID][dynamic; 16]Variational_Int),
	   isothermal_bcs = make(map[fem.Boundary_ID]Isothermal_Int),
	}
	// initialize those fixed capacity guys
	for _, id in mesh.boundary_names {d.variational_bcs[id] = {}}
	for _, id in mesh.section_names {d.sources[id] = {}}

    return d
}


Isothermal_Int :: struct {
	procedure: fem.DOF_Functional_Scalar_Proc,
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


Variational_Proc :: #type proc(
	mapped: fem.Mapped_Element,
	time: f64,
	data: rawptr,
	temperature: []f64,
	out: Variational_BC,
)

Source_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Source)

Material_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Material)

Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Scalar_Proc,
}

empty_source :: proc(n_p: int) -> Source {return {make([]f64, n_p), make([]f64, n_p)}}

empty_material :: proc(n_p: int) -> Material {
	return {make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p)}
}

empty_variational :: proc(n_p: int) -> Variational_BC {
	return {make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p)}
}


