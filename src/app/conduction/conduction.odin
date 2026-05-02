package conduction

import "core:fmt"
import "core:mem/virtual"

import "../../fem"
import "../../fem/infra"
import "../../fem/serialization"

Model_Parameters :: struct {
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][]Source_Int,
	variational_bcs: map[fem.Boundary_ID][]Variational_Int,
	isothermal_bcs:  map[fem.Boundary_ID]Isothermal_Int,
	ics:             map[fem.Section_ID]f64,
}

Isothermal_Int :: struct {
	procedure: Isothermal_Proc,
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

Isothermal_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr) -> Isothermal_BC

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

Isothermal_BC :: f64

Source :: struct {
	Q:   []f64,
	d_Q: []f64,
}

Material :: struct {
	k, rho, cp:       []f64,
	d_k, d_rho, d_cp: []f64,
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
	ts: fem.Timestep,
	u: fem.Vector,
	cm: []bool,
	allocator := context.allocator,
) {
	context.allocator = allocator

	for _, id in mesh.boundary_names {
		iso_int := params.isothermal_bcs[id] or_continue

		fem.system_mark_boundary(system, T_handle, id, cm)

		bdi := fem.boundary_dof_iter_create(system, T_handle, id)

		for bdof in fem.boundary_dof_iter_next(&bdi) {
			element := mesh.elements[bdof.element_id]
			rule := fem.basis_dof_functional_rule(element.type, fem.system_var_bd(system, T_handle), bdof.basis_dof)

			mapped := fem.map_element(element, &rule)
			defer fem.mapped_destroy(&mapped)

			value := iso_int.procedure(mapped, ts.time, iso_int.data)

			u[bdof.gdof] = fem.basis_dof_functional(
				element.type,
				fem.system_var_bd(system, T_handle),
				bdof.basis_dof,
				{value},
			)
		}
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

	quad := fem.map_quadrature(element, {.Interior, .Quad_3, 0})
	space := fem.basis_grad_space(quad, fem.system_var_bd(system, T_handle), .Scalar)

	source := empty_source(fem.space_points(space))
	material := empty_material(fem.space_points(space))

	T := fem.evaluate_var(space, u_coeffs)
	grad_T := fem.evaluate_var_gradient(space, quad, u_coeffs)

	mat := params.materials[element.section]

	mat.procedure(quad, 0, mat.data, T, material)

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

	// bcs := empty_variational(fem.space_points(space))
}


solve :: proc(model_params: Model_Parameters, mesh: fem.Mesh, order: fem.Basis_Order, arena: ^virtual.Arena) {
	context.allocator = virtual.arena_allocator(arena)

	sys_desc := fem.System_Description{}

	T_handle := fem.description_add_variable(&sys_desc, {bd = {.Lagrange, order}, components = 1})

	fem.description_couple(&sys_desc, T_handle, T_handle)

	system := fem.system_from_description(mesh, sys_desc, context.allocator)

	cm := fem.system_constraint_mask(system)
	ics := fem.system_vector(system)


	ts, ts_state := fem.timestepper_create(0, 1, 0.1, ics)
	fem.timestepper_set_scheme(&ts, &ts_state, T_handle, .BDF2)

	// ts, ts_state := fem.timestepper_create_steady(ics, system)

	ni_state := fem.nli_create_state(system)

	ca: infra.Checkpoint_Allocator
	infra.ca_init(&ca)

	context.allocator = infra.ca_allocator(&ca)


	for step in fem.timestepper_step(&ts, ts_state, ni_state.solution, system) {
		apply_constraints(model_params, system, T_handle, mesh, step, ni_state.solution, cm)

		ni := fem.nli_create(1e-6, 10)

		lin_tol := fem.nli_linear_tolerance(ni)

		for iter, should_solve in fem.nli_step(&ni, ni_state) {
			infra.ca_check(&ca); defer infra.ca_rewind(&ca)

			if should_solve {
				fem.sparse_solve(ni_state.tangent, ni_state.update, ni_state.residual, tol = lin_tol, kind = .CG_SA)
				fem.nli_update(&ni, ni_state)
			}

			u_dot, _ := fem.timestepper_derivatives(&ts, ts_state, ni_state.solution, system)

			for element, id in mesh.elements {
				infra.ca_check(&ca); defer infra.ca_rewind(&ca)

				ls := fem.system_local_problem(system, fem.Entity_ID(id))

				weak_form(model_params, system, ls, T_handle, ni_state.solution, u_dot, step, element)

				fem.system_scatter(system, fem.Entity_ID(id), ni_state.tangent, ni_state.residual, cm, ls)

			}

			fem.system_finalize_constraints(system, ni_state.tangent, ni_state.residual, cm)

		}
		od: serialization.Output_Variable_Data = {system, T_handle, ni_state.solution}
		of := serialization.output_field_from_system_variable(fem.Grad_Space(.Scalar), &od, "temperature")

		viz_mesh := serialization.vtk_create_visualization_mesh(mesh, order)

		serialization.write_vtu(fmt.aprintf("t_%d.vtu", step.step), mesh, viz_mesh, {of})

	}
}




/*
iter 0:
	- build, set initial residual.

iter 1:
	- solve first, update, rebuild, check convergence.
*/
