package linear_elasticity

import "core:fmt"
import "core:log"
import "core:math"
import "core:mem/virtual"

import "../../fem"
import "../../fem/infra"
import fio "../../fem/serialization"

Model_Parameters :: struct {
	materials:       map[fem.Section_ID]Material_Int,
	sources:         map[fem.Section_ID][]Source_Int,
	variational_bcs: map[fem.Boundary_ID][]Variational_Int,
	fixed_bcs:       map[fem.Boundary_ID]Fixed_Int,
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


Variational_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Variational_BC)

Source_Proc :: #type proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Source)

Material_Proc :: #type proc(
	mapped: fem.Mapped_Element,
	time: f64,
	current_strain: []fem.Voigt6,
	data: rawptr,
	out: Material,
)

Variational_BC :: struct {
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

Initial_Condition_Int :: struct {
	data:      rawptr,
	procedure: fem.DOF_Functional_Vector_Proc,
}

BC_Int :: union {
	Fixed_Int,
	Variational_Int,
}

@(private = "file")
empty_source :: proc(n_p: int) -> Source {return {make([]fem.Vec3, n_p)}}

@(private = "file")
empty_material :: proc(n_p: int) -> Material {
	return {make([]fem.Voigt6, n_p), make([]fem.Voigt6x6, n_p), make([]f64, n_p)}
}

@(private = "file")
empty_variational :: proc(n_p: int) -> Variational_BC {
	return {make([]fem.Vec3, n_p), make([]fem.Vec3, n_p), make([]f64, n_p)}
}

apply_constraints :: proc(
	params: Model_Parameters,
	system: fem.System,
	U_handle: fem.Var_Handle,
	mesh: fem.Mesh,
	time: f64,
	u: fem.Vector,
	cm: []bool,
	allocator := context.allocator,
) {
	for _, id in mesh.boundary_names {
		fixed := params.fixed_bcs[id] or_continue
		fem.system_apply_boundary_functional(
			system,
			U_handle,
			mesh,
			id,
			time,
			u,
			cm,
			fixed.procedure,
			fixed.data,
			allocator,
		)
	}
}

apply_ics :: proc(
	ics: map[fem.Section_ID]Initial_Condition_Int,
	system: fem.System,
	U_handle: fem.Var_Handle,
	start_time: f64,
	mesh: fem.Mesh,
	u: fem.Vector,
	allocator := context.allocator,
) {
	for section_id, ic in ics {
		fem.system_project_dofs(system, U_handle, mesh, section_id, start_time, u, ic.procedure, ic.data, allocator)
	}
}


to_engineering_shear :: proc(eps: fem.Voigt6) -> fem.Voigt6 {
	return {eps[0], eps[1], eps[2], 2 * eps[3], 2 * eps[4], 2 * eps[5]}
}

// Written assuming arena backed allocator, this function does not attempt to free any allocated memory.
weak_form :: proc(
	params: Model_Parameters,
	system: fem.System,
	ls: fem.Local_System,
	U_handle: fem.Var_Handle,
	u: fem.Vector,
	u_dot: fem.Vector,
	u_ddot: fem.Vector,
	ts: fem.Timestep,
	element: fem.Mesh_Element,
	allocator := context.allocator,
) {
	context.allocator = allocator

	u_coeffs := fem.system_gather_var_coeffs(system, U_handle, element.id, u)
	u_ddot_coeffs := fem.system_gather_var_coeffs(system, U_handle, element.id, u_ddot)

	bd := fem.system_var_bd(system, U_handle)
	quad_rule := fem.infer_quadrature(bd.order)

	quad := fem.map_quadrature(element, {.Interior, quad_rule, 0})
	space := fem.basis_grad_space(quad, bd, .Vector)

	source := empty_source(fem.space_points(space))
	material := empty_material(fem.space_points(space))

	current_strain := fem.evaluate_var_symmetric_gradient(space, quad, u_coeffs)
	for qp in 0 ..< fem.space_points(space) {
		current_strain[qp] = to_engineering_shear(current_strain[qp])
	}

	mat := params.materials[element.section]
	mat.procedure(quad, ts.time, current_strain, mat.data, material)

	for source_int in params.sources[element.section] {
		source_int.procedure(quad, ts.time, source_int.data, source)
	}

	for qp in 0 ..< fem.space_points(space) {
		C := material.constitutive_tensor[qp]
		for test in 0 ..< fem.space_arity(space) {
			residual: f64
			defer fem.local_system_rhs_add(ls, U_handle, test, residual)

			B_test := to_engineering_shear(fem.space_symmetric_gradient(space, quad, qp, test))
			N_test := fem.space_value(space, qp, test)

			// internal virtual work
			residual -= fem.dot(B_test, material.stress[qp]) * fem.dV(quad, qp)

			// body force
			residual += fem.dot(N_test, source.body[qp]) * fem.dV(quad, qp)

			for trial in 0 ..< fem.space_arity(space) {
				jacobian: f64
				defer fem.local_system_mat_add(ls, U_handle, test, U_handle, trial, jacobian)

				B_trial := to_engineering_shear(fem.space_symmetric_gradient(space, quad, qp, trial))
				N_trial := fem.space_value(space, qp, trial)

				// stiffness

				CB_trial := C * B_trial

				jacobian += fem.dot(B_test, CB_trial) * fem.dV(quad, qp)

				// mass
				M_ij := fem.dot(N_test, N_trial) * material.density[qp] * fem.dV(quad, qp)

				residual -= M_ij * u_ddot_coeffs[trial]
				jacobian += M_ij * ts.du_du_ddot[U_handle]
			}
		}
	}

	// variational BCs
	for facet in element.boundary_facets {
		quad := fem.map_quadrature(element, {.Surface, quad_rule, facet})
		space := fem.basis_grad_space(quad, bd, .Vector)

		u_coeffs := fem.system_gather_var_coeffs(system, U_handle, element.id, u)

		id := element.boundary_ids[facet]
		bcs_ints := params.variational_bcs[id] or_continue


		bcs := empty_variational(fem.space_points(space))
		for bc in bcs_ints {
			bc.procedure(quad, ts.time, bc.data, bcs)
		}

		u_at_qp := fem.evaluate_var(space, u_coeffs)

		for qp in 0 ..< fem.space_points(space) {
			for test in 0 ..< fem.space_arity(space) {
				residual: f64
				defer fem.local_system_rhs_add(ls, U_handle, test, residual)

				N_test := fem.space_value(space, qp, test)

				// traction
				residual += fem.dot(N_test, bcs.t[qp]) * fem.dS(quad, qp)

				// spring
				residual += fem.dot(N_test, bcs.k_s[qp] * (bcs.u_ref[qp] - u_at_qp[qp])) * fem.dS(quad, qp)

				for trial in 0 ..< fem.space_arity(space) {
					jacobian: f64
					defer fem.local_system_mat_add(ls, U_handle, test, U_handle, trial, jacobian)

					N_trial := fem.space_value(space, qp, trial)

					// spring stiffness
					jacobian += bcs.k_s[qp] * fem.dot(N_test, N_trial) * fem.dS(quad, qp)
				}
			}
		}
	}
}


Von_Mises_Data :: struct {
	system:       fem.System,
	var:          fem.Var_Handle,
	displacement: fem.Vector,
	materials:    map[fem.Section_ID]Material_Int,
}

von_mises_output :: proc(data: ^Von_Mises_Data) -> (of: fio.Output_Field) {
	of.friendly_name = "Von Mises Stress"
	of.components = 1
	of.data = data

	of.value_provider = proc(
		mapped: fem.Mapped_Element,
		time: f64,
		data: rawptr,
		out: [][fio.MAX_OUTPUT_FIELD_COMPONENTS]f64,
	) {
		od := cast(^Von_Mises_Data)data

		bd := fem.system_var_bd(od.system, od.var)

		space := fem.basis_grad_space(mapped, bd, .Vector)
		u_coeffs := fem.system_gather_var_coeffs(od.system, od.var, mapped.element.id, od.displacement)
		strain := fem.evaluate_var_symmetric_gradient(space, mapped, u_coeffs)

		material := empty_material(fem.space_points(space))
		mat := od.materials[mapped.element.section]
		mat.procedure(mapped, time, strain, mat.data, material)

		for p in 0 ..< fem.space_points(space) {
			s := to_engineering_shear(material.stress[p])
			sxx, syy, szz := s[0], s[1], s[2]
			sxy, syz, sxz := s[3], s[4], s[5]
			vm := math.sqrt(
				0.5 *
				((sxx - syy) * (sxx - syy) +
						(syy - szz) * (syy - szz) +
						(szz - sxx) * (szz - sxx) +
						6.0 * (sxy * sxy + syz * syz + sxz * sxz)),
			)
			out[p][0] = vm
		}
	}

	return of

}
