// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
/*
PDE model for 3D elastic materials. Non-linear support is limited as it seems less useful
under the small strain assumption which is what this model assumes, for now.

This model may evolve to include plasticity, hyperelasticity etc or that may move into a seperate model.

Sources, bcs etc may change signature as time goes on, not sure if current_strain is a safer bet to use everywhere.
*/
package models

import "core:log"
import "core:math/linalg"
import "core:math"
import "core:mem"

import fem "../fe_core"
import "../la"
import fio "../serialization"

import "core:slice"

Elastic_Material_Int :: struct {
	procedure: Elastic_Material_Proc,
	data:      rawptr,
}

// out is safe to overwrite, no need to accumulate.
Elastic_Material_Proc :: #type proc(
	ctx: fem.Element_Context,
	current_strain: []fem.Voigt6,
	data: rawptr,
	out: Elastic_Material,
)

Elastic_Material :: struct {
	stress:              []fem.Voigt6,
	constitutive_tensor: []fem.Voigt6x6,
}

Elastic_Source_Int :: struct {
	procedure: Elastic_Source_Proc,
	data:      rawptr,
}

// accumulate into out
Elastic_Source_Proc :: #type proc(
	ctx: fem.Element_Context,
	current_time: f64,
	current_soln: []fem.Vec3,
	data: rawptr,
	out: Elastic_Source,
)

Elastic_Source :: struct {
	F: []fem.Vec3,
}

Elastic_BC_Int :: struct {
	procedure: Elastic_BC_Proc,
	data:      rawptr,
}

// out should be accumulated into not overwritten to allow multiple variational bc's per facet.
Elastic_BC_Proc :: #type proc(
	ctx: fem.Element_Context,
	point_start, point_end: int,
	current_time: f64,
	current_soln: []fem.Vec3,
	data: rawptr,
	out: Elastic_BC,
)

Elastic_BC :: struct {
	t:     []fem.Vec3, //prescribed traction
	k_s:   []f64, //spring/foundation stiffness
	u_ref: []fem.Vec3, //reference displacement for spring
}

blank_elastic_material :: proc(num_points: int, allocator := context.allocator) -> Elastic_Material {
	context.allocator = allocator
	return {stress = make([]fem.Voigt6, num_points), constitutive_tensor = make([]fem.Voigt6x6, num_points)}
}

blank_elastic_source :: proc(num_points: int, allocator := context.allocator) -> Elastic_Source {
	context.allocator = allocator
	return {F = make([]fem.Vec3, num_points)}
}

blank_elastic_bc :: proc(num_points: int, allocator := context.allocator) -> Elastic_BC {
	context.allocator = allocator
	return {t = make([]fem.Vec3, num_points), k_s = make([]f64, num_points), u_ref = make([]fem.Vec3, num_points)}
}

Elasticity_Params :: struct {
	soln_order:              fem.Order,
	dirichlet_displacements: map[fem.Boundary_ID]fem.Vec3,
	materials:               map[fem.Section_ID]Elastic_Material_Int,
	sources:                 map[fem.Section_ID][]Elastic_Source_Int,
	variational_bcs:         map[fem.Boundary_ID][]Elastic_BC_Int,
	ics:                     map[fem.Section_ID]fem.Vec3,
	displ_fh:                fem.Field_Handle,
}


model_elasticity :: proc(p: ^Elasticity_Params) -> Model {
	define_layout :: proc(m: ^Model, mesh: fem.Mesh, allocator: mem.Allocator) -> (fem.Layout, fem.Constraint_Mask) {
		params := cast(^Elasticity_Params)m.data
		context.allocator = allocator

		displ_fd := fem.Field_Descriptor{.LV, params.soln_order}

		layout := fem.layout_create(allocator)

		params.displ_fh = fem.layout_add_field(&layout, mesh, displ_fd)
		fem.layout_couple(&layout, mesh, params.displ_fh, params.displ_fh)

		cm := fem.layout_make_empty_constraint_mask(&layout)

		iso_bnds, _ := slice.map_keys(params.dirichlet_displacements)
		defer delete(iso_bnds)

		field_mark_constraints(mesh, layout, cm, params.displ_fh, iso_bnds)

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
		params := cast(^Elasticity_Params)m.data

		context.allocator = allocator

		for id, dofs in layout.bound_maps[params.displ_fh] {
			displacement := params.dirichlet_displacements[id] or_continue
			for d in dofs {
				basis := fem.basis_create(fem.Basis_LV, mesh.elements[d.element], params.soln_order)
				_, cmpnt := fem.basis_decompose_component(basis, d.local)
				global := fem.layout_global_pos(layout, params.displ_fh, d.element, d.local)
				current_iterate.values[global] = displacement[cmpnt]
			}
		}
	}

	set_initial_conditions :: proc(m: ^Model, layout: fem.Layout, mesh: fem.Mesh, initial_solution: la.Block_Vector) {
		params := cast(^Elasticity_Params)m.data

		for element in mesh.elements {
			ic_val := params.ics[element.section] or_else 0
			basis := fem.basis_create(fem.Basis_LV, element, params.soln_order)
			for local_dof in 0 ..< basis.arity {
				_, cmpnt := fem.basis_decompose_component(basis, local_dof)
				global := fem.layout_global_pos(layout, params.displ_fh, element.id, local_dof)
				initial_solution.values[global] = ic_val[cmpnt]
			}
		}
	}

	to_engineering_shear :: proc(eps: fem.Voigt6) -> fem.Voigt6 {
		return {eps[0], eps[1], eps[2], 2 * eps[3], 2 * eps[4], 2 * eps[5]}
	}

	build_local_problem :: proc(
		m: ^Model,
		layout: fem.Layout,
		element: fem.Mesh_Element,
		current_iterate: la.Block_Vector,
		previous_iterate: la.Block_Vector,
		tc: Time_Context,
		R: la.Block_Vector,
		J: la.Block_Dense_Matrix,
		allocator: mem.Allocator,
	) {
		params := cast(^Elasticity_Params)m.data
		context.allocator = allocator

		basis := fem.basis_create(fem.Basis_LV, element, params.soln_order)

		quad, count := fem.quadrature_for(
			element,
			fem.infer_quadrature_rule(element, params.soln_order, .Interior),
			basis.geometry_required,
		)

		current_soln := make([]fem.Vec3, len(quad.points))
		current_strain := make([]fem.Voigt6, len(quad.points))
		material := blank_elastic_material(len(quad.points))
		sources := blank_elastic_source(len(quad.points))

		field_evaluate(layout, quad, params.displ_fh, basis, current_iterate, current_soln)
		evaluate_symmetric_gradient(layout, quad, params.displ_fh, basis, current_iterate, current_strain)

		for i in 0 ..< len(current_strain) {
			current_strain[i] = to_engineering_shear(current_strain[i])
		}

		mat_int := params.materials[element.section]
		mat_int.procedure(quad, current_strain, mat_int.data, material)

		source_ints := params.sources[element.section]
		for source_int in source_ints {
			source_int.procedure(quad, tc.current_time, current_soln, source_int.data, sources)
		}

		for qp in 0 ..< count {
			stress := material.stress[qp]
			C := material.constitutive_tensor[qp]
			strain := current_strain[qp]

			for test in 0 ..< basis.arity {
				eps_v := to_engineering_shear(fem.basis_sym_grad(basis, quad, qp, test))
				v := fem.basis_value(basis, quad, qp, test)
				r_idx := la.idx(R, int(params.displ_fh), test)

				// internal force
				R.values[r_idx] += -fem.voigt_dot(stress, eps_v) * fem.dV(quad, qp)

				// body force
				R.values[r_idx] += (linalg.dot(sources.F[qp], v) * fem.dV(quad, qp))

				for trial in 0 ..< basis.arity {
					eps_u := to_engineering_shear(fem.basis_sym_grad(basis, quad, qp, trial))

					j_idx := la.idx(J, int(params.displ_fh), int(params.displ_fh), test, trial)

					C_eps_u: fem.Voigt6
					fem.voigt_gemv(C, eps_u, &C_eps_u, 1, 0)
					J_ij := fem.voigt_dot(eps_v, C_eps_u) * fem.dV(quad, qp)

					J.values[j_idx] += J_ij
				}
			}
		}
		// boundary terms
		bound_facets := fem.boundary_set(element)
		if bound_facets == {} {
			return
		}

		quad, count = fem.quadrature_for(
			element,
			fem.infer_quadrature_rule(element, params.soln_order, .Facet),
			basis.geometry_required,
		)

		current_soln = make([]fem.Vec3, len(quad.points))
		field_evaluate(layout, quad, params.displ_fh, basis, current_iterate, current_soln)
		for facet in bound_facets {
			id := element.boundaries[facet].?
			bcs := params.variational_bcs[id] or_continue

			start, end := fem.facet_quad(quad, facet)

			bc_coeff := blank_elastic_bc(end - start)
			for bc in bcs {
				bc.procedure(quad, start, end, tc.current_time, current_soln, bc.data, bc_coeff)
			}

			for qp in start ..< end {
				t := bc_coeff.t[qp - start]
				k_s := bc_coeff.k_s[qp - start]
				n := fem.ctx_facet_normal(quad, qp)
				u_n := linalg.dot(current_soln[qp], n)
				u_ref_n := linalg.dot(bc_coeff.u_ref[qp - start], n)

				for test in 0 ..< basis.arity {
					v := fem.basis_value(basis, quad, qp, test)
					v_n := linalg.dot(v, n)
					r_idx := la.idx(R, int(params.displ_fh), test)

					// traction
					R.values[r_idx] += linalg.dot(t, v) * fem.dS(quad, qp)

					// spring
					R.values[r_idx] += k_s * (u_ref_n - u_n) * v_n * fem.dS(quad, qp)

					for trial in 0 ..< basis.arity {
						u_v := fem.basis_value(basis, quad, qp, trial)
						u_v_n := linalg.dot(u_v, n)
						j_idx := la.idx(J, int(params.displ_fh), int(params.displ_fh), test, trial)
						J.values[j_idx] += k_s * u_v_n * v_n * fem.dS(quad, qp)
					}
				}
			}
		}
	}

	output_fields :: proc(
		m: ^Model,
		layout: fem.Layout,
		current_iterate: la.Block_Vector,
		tc: Time_Context,
		output_fields: ^[dynamic]fio.Output_Field,
	) {
		params := cast(^Elasticity_Params)m.data
		append(output_fields, field_to_output_field(layout, params.displ_fh, current_iterate, "Displacement"))

		Von_Mises_Info :: struct {
			displacement:   []f64,
			layout: fem.DOF_Layout,
			order:  fem.Order,
			material: map[fem.Section_ID]Elastic_Material_Int,
		}

		von_mises_info := new(Von_Mises_Info)

		von_mises_info^ = {
			displacement = la.block_vector_view(current_iterate, int(params.displ_fh)),
			layout = layout.dof_layouts[params.displ_fh],
			order = params.soln_order,
			material = params.materials,
		}

		von_mises_proc :: proc(ctx: fem.Element_Context, data: rawptr, out: [][fio.MAX_OUTPUT_FIELD_COMPONENTS]f64) {
			info := cast(^Von_Mises_Info)data

			current_strain := make([]fem.Voigt6, len(ctx.points))
			material := blank_elastic_material(len(ctx.points))
			basis := fem.basis_create(fem.Basis_LV, ctx.element, info.order)

			for pi in 0 ..< len(ctx.points) {
				for dof in 0 ..< basis.arity {
					coeff := info.displacement[info.layout.mapping[ctx.element.id][dof]]
					current_strain[pi] += fem.basis_sym_grad(basis, ctx, pi, dof) * coeff
				}
			}

			for i in 0 ..< len(current_strain) {
				current_strain[i] = to_engineering_shear(current_strain[i])
			}

			mat_int := info.material[ctx.element.section]
			mat_int.procedure(ctx, current_strain, mat_int.data, material)


			for pi in 0..<len(ctx.points) {
				stress := material.stress[pi]
				xx, yy, zz := stress[0], stress[1], stress[2]
	    		xy, yz, xz := stress[3], stress[4], stress[5]
	    		von_mises := math.sqrt(
	        	   0.5 * ((xx-yy)*(xx-yy) + (yy-zz)*(yy-zz) + (zz-xx)*(zz-xx) +
	               6.0 * (xy*xy + yz*yz + xz*xz)),
	    		)
	    		out[pi][0] = von_mises
			}
			}

			output_field := fio.Output_Field{
				friendly_name = "Von Mises Stress",
				components = 1,
				data = von_mises_info,
				value_provider = von_mises_proc,
				geometry_required = fem.Basis_LV_Geometry,
			}

			append(output_fields, output_field)

	}

	return Model {
		data = p,
		define_layout = define_layout,
		set_initial_conditions = set_initial_conditions,
		respect_constraints = respect_constraints,
		build_local_problem = build_local_problem,
		output_fields = output_fields,
	}

}
