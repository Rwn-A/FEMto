package models


import "core:math/linalg"
import "core:mem"

import fem "../fe_core"
import "../la"
import fio "../serialization"

import "core:slice"

Conduction_Params :: struct {
	soln_order:      fem.Order,
	isothermal_bnds: map[fem.Boundary_ID]f64,

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

		for qp in 0 ..< count {
			for test in 0 ..< basis.arity {
				grad_v := fem.basis_gradient(basis, quad, qp, test)
				for trial in 0 ..< basis.arity {
					grad_u := fem.basis_gradient(basis, quad, qp, trial)

					J_ij := linalg.dot(grad_u, grad_v) * 10 * fem.dV(quad, qp)
					u_j := field_coeff(layout, params.temp_fh, current_iterate, element.id, trial)

					R.values[la.idx(R, int(params.temp_fh), test)] += -J_ij * u_j
					J.values[la.idx(J, int(params.temp_fh), int(params.temp_fh), test, trial)] += J_ij
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
