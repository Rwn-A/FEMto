package models


import "core:math/linalg"
import "core:mem"
import "core:log"

import fem "../fe_core"
import "../la"
import fio "../serialization"

import "core:slice"

Conduction_Material_Int :: struct {
	procedure: Conduction_Material_Proc,
	data: rawptr,
}

// out is safe to ovewrite, no need to accumulate.
Conduction_Material_Proc :: #type proc(ctx: fem.Element_Context, current_time: f64, current_soln: []f64, data: rawptr, out: Conduction_Material)

Conduction_Material :: struct {
	k, rho, cp: []f64,
	d_k, d_rho, d_cp: []f64, // derivatives w/ respect to T.
}

Conduction_Source_Int :: struct {
	procedure: Conduction_Source_Proc,
	data: rawptr,
}

// NOTE: out is not to be overwritten but accumulated into, this allows multiple sources to contribute.
Conduction_Source_Proc :: #type proc(ctx: fem.Element_Context, current_time: f64, current_soln: []f64, data: rawptr, out: Conduction_Source)

Conduction_Source :: struct {
	Q: []f64,
	d_Q: []f64, // derivatives w/ respect to T.
}

// specifically variatonial bc's, dirichlet is handled elsewhere.
Conduction_BC_Int :: struct {
	procedure: Conduction_BC_Proc,
	data: rawptr,
}

// out is safe to overwrite. out is only as big as `point_end` - `point_start` as not all facet quad points make up the boundary.
// current solution is sized as normal.
Conduction_BC_Proc :: #type proc(ctx: fem.Element_Context, point_start, point_end: int, current_time: f64, current_soln: []f64, data: rawptr, out: Conduction_BC)

Conduction_BC :: struct {
	q, h, T_amb: []f64,
	d_q, d_h, d_T_amb:    []f64,
}


blank_conduction_material :: proc(num_points: int, allocator := context.allocator) -> Conduction_Material {
	context.allocator = allocator
	return {
		k = make([]f64, num_points),
		rho = make([]f64, num_points),
		cp = make([]f64, num_points),
		d_k = make([]f64, num_points),
		d_rho = make([]f64, num_points),
		d_cp = make([]f64, num_points)
	}
}


blank_conduction_source :: proc(num_points: int, allocator := context.allocator) -> Conduction_Source {
	context.allocator = allocator
	return {
		Q = make([]f64, num_points),
		d_Q = make([]f64, num_points),
	}
}

blank_conduction_bc :: proc(num_points: int, allocator := context.allocator) -> Conduction_BC {
	return {
		q = make([]f64, num_points),
		h = make([]f64, num_points),
		T_amb = make([]f64, num_points),
		d_q = make([]f64, num_points),
		d_h = make([]f64, num_points),
		d_T_amb = make([]f64, num_points)
	}
}

Conduction_Params :: struct {
	soln_order:      fem.Order,
	isothermal_bnds: map[fem.Boundary_ID]f64,

	materials: map[fem.Section_ID]Conduction_Material_Int,
	sources: map[fem.Section_ID][]Conduction_Source_Int,
	variational_bcs: map[fem.Boundary_ID]Conduction_BC_Int,

	// state.
	temp_fh:         fem.Field_Handle,
}

model_conduction :: proc(p: ^Conduction_Params) -> Model {
	define_layout :: proc(m: ^Model, mesh: fem.Mesh, allocator: mem.Allocator) -> (fem.Layout, fem.Constraint_Mask) {
		params := cast(^Conduction_Params)m.data
		context.allocator = allocator

		temp_fd := fem.Field_Descriptor{.LS, params.soln_order}

		layout := fem.layout_create(allocator)

		params.temp_fh = fem.layout_add_field(&layout, mesh, temp_fd)
		fem.layout_couple(&layout, mesh, params.temp_fh, params.temp_fh)

		cm := fem.layout_make_empty_constraint_mask(&layout)

		iso_bnds, _ := slice.map_keys(params.isothermal_bnds)
		defer delete(iso_bnds)

		field_mark_constraints(mesh, layout, cm, params.temp_fh, iso_bnds)

		return layout, cm
	}

	respect_constraints :: proc(
		m: ^Model,
		layout: fem.Layout,
		mesh: fem.Mesh,
		mask: fem.Constraint_Mask,
		current_iterate: la.Block_Vector,
		tc: Time_Context,
		allocator: mem.Allocator,
	) {
		params := cast(^Conduction_Params)m.data

		context.allocator = allocator

		for id, dofs in layout.bound_maps[params.temp_fh] {
			iso_bound := params.isothermal_bnds[id] or_continue
			for d in dofs {
				global := fem.layout_global_pos(layout, params.temp_fh, d.element, d.local)
				current_iterate.values[global] = iso_bound
			}
		}
	}

	build_local_problem :: proc(
		m: ^Model,
		layout: fem.Layout,
		element: fem.Mesh_Element,
		current_iterate: la.Block_Vector,
		previous_iterate: la.Block_Vector, //from timestepping loop.
		tc: Time_Context,
		R: la.Block_Vector,
		J: la.Block_Dense_Matrix,
		allocator: mem.Allocator,
	) {
		params := cast(^Conduction_Params)m.data
		context.allocator = allocator

		basis := fem.basis_create(fem.Basis_LS, element, params.soln_order)

		quad, count := fem.quadrature_for(
			element,
			fem.infer_quadrature_rule(element, params.soln_order, .Interior),
			basis.geometry_required,
		)

		current_soln := make([]f64, len(quad.points))
		material := blank_conduction_material(len(quad.points))
		sources := blank_conduction_source(len(quad.points))

		evaluate_lagrange_scalar(layout, quad, params.temp_fh, basis, current_iterate, current_soln)
		mat_int := params.materials[element.section]
		mat_int.procedure(quad, tc.current_time, current_soln, mat_int.data, material)

		source_ints := params.sources[element.section]
		for source_int in source_ints {
			source_int.procedure(quad, tc.current_time, current_soln, source_int.data, sources)
		}


		dt_inv := (1 / tc.timestep) if tc.is_transient else 0 // zero if steady, 1/dt if transient

		for qp in 0 ..< count {
			for test in 0 ..< basis.arity {
				grad_v := fem.basis_gradient(basis, quad, qp, test)
				v := fem.basis_value(basis, quad, qp, test)
				for trial in 0 ..< basis.arity {

					grad_u := fem.basis_gradient(basis, quad, qp, trial)
					u := fem.basis_value(basis, quad, qp, trial)
					u_j := field_coeff(layout, params.temp_fh, current_iterate, element.id, trial)
					u_j_old := field_coeff(layout, params.temp_fh, previous_iterate, element.id, trial)

					// stiffness
					J_ij := linalg.dot(grad_u, grad_v) * material.k[qp] * fem.dV(quad, qp)
					R.values[la.idx(R, int(params.temp_fh), test)] += -J_ij * u_j
					J.values[la.idx(J, int(params.temp_fh), int(params.temp_fh), test, trial)] += J_ij

					// mass (timestepping)
					M_ij := u * v * material.rho[qp] * material.cp[qp] * fem.dV(quad, qp)
					R.values[la.idx(R, int(params.temp_fh), test)] += -M_ij * (u_j - u_j_old) * dt_inv
					J.values[la.idx(J, int(params.temp_fh), int(params.temp_fh), test, trial)] += M_ij * dt_inv
				}
				R.values[la.idx(R, int(params.temp_fh), test)] += sources.Q[qp] * fem.basis_value(basis, quad, qp, test) * fem.dV(quad, qp)
			}
		}

		bound_facets := fem.boundary_set(element)
		if bound_facets == {} {return} // conduction has no surface terms other than bc's.

		quad, count = fem.quadrature_for(
			element,
			fem.infer_quadrature_rule(element, params.soln_order, .Facet),
			basis.geometry_required,
		)

		// quad points different so new solution
		current_soln = make([]f64, len(quad.points))
		evaluate_lagrange_scalar(layout, quad, params.temp_fh, basis, current_iterate, current_soln)


		for facet in bound_facets {
			id := element.boundaries[facet].?
			bc := params.variational_bcs[id] or_continue

			start, end := fem.facet_quad(quad, facet)
			bc_coeff := blank_conduction_bc(end - start)
			bc.procedure(quad, start, end, tc.current_time, current_soln, bc.data, bc_coeff)

			for qp in start..<end {
				q := bc_coeff.q[qp - start]
				h := bc_coeff.h[qp - start]
				T_amb := bc_coeff.T_amb[qp - start]
				for test in 0..<basis.arity {
					v := fem.basis_value(basis, quad, qp, test)
			        R.values[la.idx(R, int(params.temp_fh), test)] += (q + h * T_amb - h * current_soln[qp]) * v * fem.dS(quad, qp)
					for trial in 0..<basis.arity {
						 u := fem.basis_value(basis, quad, qp, trial)
            			 J.values[la.idx(J, int(params.temp_fh), int(params.temp_fh), test, trial)] += h * u * v * fem.dS(quad, qp)
					}
				}
			}
		}
	}

	output_fields :: proc(
		m: ^Model,
		layout: fem.Layout,
		current_iterate: la.Block_Vector,
		output_fields: ^[dynamic]fio.Output_Field,
	) {
		params := cast(^Conduction_Params)m.data
		append(output_fields, field_to_output_field(layout, params.temp_fh, current_iterate, "Temperature"))
	}

	return Model {
		data = p,
		define_layout = define_layout,
		respect_constraints = respect_constraints,
		build_local_problem = build_local_problem,
		output_fields = output_fields,
	}

}
