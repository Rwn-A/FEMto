package inc_flow

import "../../fem"

import "core:math/linalg"
import "core:math"
import "core:log"

Bulk_Response :: struct {
    material: Material,
    source: Source,
    v_coeffs: []f64,
    p_coeffs: []f64,
    velocity: []fem.Vec3,
    grad_v: []fem.Mat3,
    grad_p: []fem.Vec3
}

Boundary_Response :: struct {
    v_variational: V_Variational,
    p_variational: P_Variational,
}

Problem_Data :: struct {
	v_ics:             map[fem.Section_ID]V_Initial_Condition_Int,
	p_ics:             map[fem.Section_ID]P_Initial_Condition_Int,
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][dynamic; 16]Source_Int,
	v_var_bcs: map[fem.Boundary_ID][dynamic; 16]V_Variational_Int,
	p_var_bcs: map[fem.Boundary_ID][dynamic; 16]P_Variational_Int,
	v_fixed_bcs:  map[fem.Boundary_ID]V_Fixed_Int,
	p_fixed_bcs:  map[fem.Boundary_ID]P_Fixed_Int,
}

empty_problem_data :: proc(mesh: fem.Mesh, allocator := context.allocator) -> (d: Problem_Data) {
	context.allocator = allocator
    d = {
    v_ics         = make(map[fem.Section_ID]V_Initial_Condition_Int),
	p_ics        =  make(map[fem.Section_ID]P_Initial_Condition_Int),
	materials   =   make(map[fem.Section_ID]Material_Int),
	sources    =    make(map[fem.Section_ID][dynamic; 16]Source_Int),
	v_var_bcs  =    make(map[fem.Boundary_ID][dynamic; 16]V_Variational_Int),
	p_var_bcs  =    make(map[fem.Boundary_ID][dynamic; 16]P_Variational_Int),
	v_fixed_bcs =   make(map[fem.Boundary_ID]V_Fixed_Int),
	p_fixed_bcs  =  make(map[fem.Boundary_ID]P_Fixed_Int),
	}
	// initialize those fixed capacity guys
	for _, id in mesh.boundary_names {d.v_var_bcs[id] = {}; d.p_var_bcs[id] = {}}
	for _, id in mesh.section_names {d.sources[id] = {}}

    return d
}


assemble_bulk :: proc(
	ls:           fem.Local_System,
	V_handle:     fem.Var_Handle,
	P_handle:     fem.Var_Handle,
	V_space:      fem.Grad_Space(.Vector),
	P_space:       fem.Grad_Space(.Scalar),
	quad:         fem.Mapped_Element,
	ts:           fem.Timestep,
	bulk:         Bulk_Response,
) {

    dimension: f64 = 3
    stabilization_param: f64 = 4

    element_volume: f64
    for qp in 0 ..< fem.space_points(P_space) {
        element_volume += fem.dV(quad, qp)
    }

	for qp in 0 ..< fem.space_points(P_space) {
	   density := bulk.material.density[qp]
	   viscosity := bulk.material.dynamic_viscosity[qp]
	    h_e := math.pow(element_volume, 1.0 / dimension)
        vel_norm := max(linalg.length(bulk.velocity[qp]), 1e-8)

        time_term := 4 / math.pow(ts.dt, 2) if ts.dt != 0  else 0
        tau := 1.0 / math.sqrt(
             time_term +
            math.pow(2.0 * vel_norm / h_e, 2) +
            math.pow(stabilization_param * viscosity / (h_e * h_e), 2)
        )

		for test in 0 ..< fem.space_arity(V_space) {
            residual: f64
            defer fem.local_system_rhs_add(ls, V_handle, test, residual)

            w := fem.space_value(V_space, qp, test)
            grad_w := fem.space_gradient(V_space, quad, qp, test)

			for trial in 0 ..< fem.space_arity(V_space) {
			    jacobian: f64
				defer fem.local_system_mat_add(ls, V_handle, test, V_handle, trial, jacobian)

                u := fem.space_value(V_space, qp, trial)
                grad_u := fem.space_gradient(V_space, quad, qp, trial)

                convection_ij := density * fem.dot(grad_u * bulk.velocity[qp], w) * fem.dV(quad, qp)
                diffusion_ij := viscosity * fem.frob(grad_u, grad_w) * fem.dV(quad, qp)

                residual +=  (convection_ij + diffusion_ij) * bulk.v_coeffs[trial]
                jacobian += convection_ij + diffusion_ij

                deriv_ij := density * fem.dot(u, w) * fem.dV(quad, qp)

                // residual += deriv_ij * V_dot_coeffs[trial]
                // jacobian += deriv_ij * du_udot
			}

            for trial in 0 ..< fem.space_arity(P_space) {
                jacobian: f64
				defer fem.local_system_mat_add(ls, V_handle, test, P_handle, trial, jacobian)

                p := fem.space_value(P_space, qp, trial)
				P_V_coupling_ij := -p * linalg.trace(grad_w) * fem.dV(quad, qp)

                residual += P_V_coupling_ij * bulk.p_coeffs[trial]
                jacobian += P_V_coupling_ij
			}

            // chat wrote the new coe for tau, no clue if that's correct

            r_supg := density * bulk.grad_v[qp] * bulk.velocity[qp] + bulk.grad_p[qp]
            SUPG_stabilization := tau * fem.dot( grad_w * bulk.velocity[qp], r_supg) * fem.dV(quad, qp)
            residual += SUPG_stabilization
		}

        for test in 0 ..< fem.space_arity(P_space) { // for u,q
            residual: f64
			defer fem.local_system_rhs_add(ls, P_handle, test, residual)

            q := fem.space_value(P_space, qp, test)
            grad_q := fem.space_gradient(P_space, quad, qp, test)

			for trial in 0 ..< fem.space_arity(V_space) {
                jacobian: f64
				defer fem.local_system_mat_add(ls, P_handle, test, V_handle, trial, jacobian)

                grad_u := fem.space_gradient(V_space, quad, qp, trial)

				continuity_ij := q * linalg.trace(grad_u) * fem.dV(quad, qp)

                pspg_ij := density * tau * linalg.dot(grad_q, grad_u * bulk.velocity[qp]) * fem.dV(quad, qp)

                residual += (continuity_ij + pspg_ij) * bulk.p_coeffs[test]
                jacobian += continuity_ij + pspg_ij
			}

            for trial in 0 ..< fem.space_arity(P_space) {
                jacobian: f64
				defer fem.local_system_mat_add(ls, P_handle, test, P_handle, trial, jacobian)

                grad_p := fem.space_gradient(P_space, quad, qp, trial)

                pspg_ij := density * tau * linalg.dot(grad_q, grad_p) * fem.dV(quad, qp)

                residual += pspg_ij * bulk.p_coeffs[trial]
                jacobian += pspg_ij
			}
		}
	}
}


assemble_boundary :: proc(
    ls:           fem.Local_System,
	V_handle:     fem.Var_Handle,
	P_handle:     fem.Var_Handle,
	v_space:        fem.Grad_Space(.Vector),
	p_space:        fem.Grad_Space(.Scalar),
	quad:         fem.Mapped_Element,
	boundary:     Boundary_Response,
) {

}

V_Variational :: struct {
    shear_stress: []fem.Vec3,
}

P_Variational :: struct {
    ngrad: []f64,
}

Source :: struct {
	body: []fem.Vec3,
}

Material :: struct {
	density: []f64,
	dynamic_viscosity: []f64,
}

V_Fixed_Int :: struct {
	procedure: fem.DOF_Functional_Vector_Proc,
	data:      rawptr,
}

P_Fixed_Int :: struct {
	procedure: fem.DOF_Functional_Scalar_Proc,
	data:      rawptr,
}

Material_Int :: struct {
	procedure: Material_Proc,
	data:      rawptr,
}

Source_Int :: struct {
    procedure: Source_Proc,
    data: rawptr,
}

V_Variational_Int :: struct {
	procedure: V_Variational_Proc,
	data:      rawptr,
}

P_Variational_Int :: struct {
	procedure: P_Variational_Proc,
	data:      rawptr,
}

V_Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Vector_Proc,
}

P_Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Scalar_Proc,
}


V_Variational_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: V_Variational)
P_Variational_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: P_Variational)
Source_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Source)

Material_Proc :: #type proc(
	mapped: fem.Mapped_Element,
	time: f64,
	data: rawptr,
	out: Material,
)


empty_source :: proc(n_p: int) -> Source {return {make([]fem.Vec3, n_p)}}

empty_material :: proc(n_p: int) -> Material {
	return {make([]f64, n_p), make([]f64, n_p)}
}

empty_v_variational :: proc(n_p: int) -> V_Variational {
	return {make([]fem.Vec3, n_p)}
}

empty_p_variational :: proc(n_p: int) -> P_Variational {
	return {make([]f64, n_p)}
}




