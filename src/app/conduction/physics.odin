package conduction

import "core:fmt"
import "core:log"
import "core:mem/virtual"

import "../../fem"
import "../../fem/infra"
import "../../fem/serialization"

Model_Parameters :: struct {
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][]Source_Int,
	variational_bcs: map[fem.Boundary_ID][]Variational_Int,
	isothermal_bcs:  map[fem.Boundary_ID]Isothermal_Int,
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

Variational_BC :: struct {
	q, h, T_amb: []f64,
	d_q, d_h:    []f64,
}

Source :: struct {
	Q:   []f64,
	d_Q: []f64,
}

Material :: struct {
	k, rho, cp:       []f64,
	d_k, d_rho, d_cp: []f64,
}

Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Scalar_Proc,
}

@(private = "file")
empty_source :: proc(n_p: int) -> Source {return {make([]f64, n_p), make([]f64, n_p)}}

@(private = "file")
empty_material :: proc(n_p: int) -> Material {
	return {make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p)}
}

@(private = "file")
empty_variational :: proc(n_p: int) -> Variational_BC {
	return {make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p), make([]f64, n_p)}
}

apply_constraints :: proc(
	params: Model_Parameters,
	system: fem.System,
	T_handle: fem.Var_Handle,
	mesh: fem.Mesh,
	time: f64,
	u: fem.Vector,
	cm: []bool,
	allocator := context.allocator,
) {
	for _, id in mesh.boundary_names {
		iso := params.isothermal_bcs[id] or_continue
		fem.system_apply_boundary_functional(system, T_handle, mesh, id, time, u, cm, iso.procedure, iso.data, allocator)
	}
}

apply_ics :: proc(
	ics: map[fem.Section_ID]Initial_Condition_Int,
	system: fem.System,
	T_handle: fem.Var_Handle,
	mesh: fem.Mesh,
	u: fem.Vector,
	allocator := context.allocator,
) {
	for section_id, ic in ics {
		// time 0 because time dependent ics dont really make sense.
		fem.system_project_dofs(system, T_handle, mesh, section_id, 0, u, ic.procedure, ic.data, allocator)
	}
}



// Written assuming arena backed allocator, this function does not attempt to free any allocated memory.
weak_form :: proc(
	params: Model_Parameters,
	system: fem.System,
	ls: fem.Local_System,
	T_handle: fem.Var_Handle,
	u: fem.Vector,
	u_dot: fem.Vector,
	ts: fem.Timestep,
	element: fem.Mesh_Element,
	allocator := context.allocator,
) {
	context.allocator = allocator

	u_coeffs := fem.system_gather_var_coeffs(system, T_handle, element.id, u)
	u_dot_coeffs := fem.system_gather_var_coeffs(system, T_handle, element.id, u_dot)

	bd := fem.system_var_bd(system, T_handle)

	quad := fem.map_quadrature(element, {.Interior, .Quad_3, 0})
	space := fem.basis_grad_space(quad, bd, .Scalar)

	source := empty_source(fem.space_points(space))
	material := empty_material(fem.space_points(space))

	T := fem.evaluate_var(space, u_coeffs)
	grad_T := fem.evaluate_var_gradient(space, quad, u_coeffs)

	mat := params.materials[element.section]

	mat.procedure(quad, ts.time, mat.data, T, material)

	for source_int in params.sources[element.section]{
	   source_int.procedure(quad, ts.time, source_int.data, T, source)
	}


	for qp in 0 ..< fem.space_points(space) {
		for test in 0 ..< fem.space_arity(space) {
			residual: f64
			defer fem.local_system_rhs_add(ls, T_handle, test, residual)

			grad_test := fem.space_gradient(space, quad, qp, test)
			val_test := fem.space_value(space, qp, test)

			// sources (linear part)
			residual += source.Q[qp] * val_test * fem.dV(quad, qp)


			for trial in 0 ..< fem.space_arity(space) {
				jacobian: f64
				defer fem.local_system_mat_add(ls, T_handle, test, T_handle, trial, jacobian)

				grad_trial := fem.space_gradient(space, quad, qp, trial)
				val_trial := fem.space_value(space, qp, trial)

				// stiffness
				K_ij := fem.dot(grad_test, grad_trial) * material.k[qp] * fem.dV(quad, qp)
				K_ij_deriv := fem.dot(grad_T[qp], grad_test) * material.d_k[qp] * val_trial * fem.dV(quad, qp)

				residual -= K_ij * u_coeffs[trial]
				jacobian += K_ij + K_ij_deriv

				// sources (non linear)
				jacobian -= source.d_Q[qp] * val_test * val_trial * fem.dV(quad, qp)

				//mass (timestepping)

				M_ij := val_trial * val_test * material.rho[qp] * material.cp[qp] * fem.dV(quad, qp)
				M_ij_deriv :=
					(val_trial *
						val_test *
						(material.d_rho[qp] * material.cp[qp] + material.rho[qp] * material.d_cp[qp])) *
					u_dot_coeffs[trial] *
					fem.dV(quad, qp)

				residual -= M_ij * u_dot_coeffs[trial]

				jacobian += M_ij * ts.du_du_dot[T_handle] - M_ij_deriv
			}
		}
	}

    // variational bcs.
	for facet in element.boundary_facets {
		quad := fem.map_quadrature(element, {.Surface, fem.infer_quadrature(bd.order), facet})
		space := fem.basis_grad_space(quad, bd, .Scalar)

		u_coeffs := fem.system_gather_var_coeffs(system, T_handle, element.id, u)
		T := fem.evaluate_var(space, u_coeffs)

		id := element.boundary_ids[facet]
		bcs_ints := params.variational_bcs[id] or_continue

		bcs := empty_variational(fem.space_points(space))

		for bc in bcs_ints{
		  bc.procedure(quad, ts.time, bc.data, T, bcs)
		}

		for qp in 0..<fem.space_points(space) {
		  for test in 0..<fem.space_arity(space) {
		      residual: f64
			  defer fem.local_system_rhs_add(ls, T_handle, test, residual)

		      val_test := fem.space_value(space, qp, test)

		      residual += (bcs.q[qp] + bcs.h[qp] * bcs.T_amb[qp] - bcs.h[qp] * u_coeffs[qp]) * val_test * fem.dS(quad, qp)

		      for trial in 0..<fem.space_arity(space) {
		            jacobian: f64
				    defer fem.local_system_mat_add(ls, T_handle, test, T_handle, trial, jacobian)

					val_trial := fem.space_value(space, qp, trial)

					jacobian += bcs.h[qp] * val_trial * val_test * fem.dS(quad, qp)
					// non linear h
					jacobian += bcs.d_h[qp] * (bcs.T_amb[qp] - u_coeffs[qp]) * val_trial * val_test * fem.dS(quad, qp)
					// non linear q
					jacobian -= bcs.d_q[qp] * val_trial * val_test * fem.dS(quad, qp)
		      }

		  }
		}

	}
}
