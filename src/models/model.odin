package models

import "core:slice"

import fem "../fe_core"
import "../infra"
import "../la"
import fio "../serialization"
import "core:log"
import "core:path/filepath"
import "core:strings"

import "core:fmt"
import "core:os/os2"

import "core:mem"

Time_Context :: struct {
	is_transient:                     bool,
	current_time, timestep, end_time: f64,
}

// which fields are outputted is configuration on the specific model, whereas this is solver config
Output_Config :: struct {
	output_dir:  string,
	frequency:   int,
	prefix:      string,
	output_rule: fem.Subcell_Rule,
}

Solver_Config :: struct {
	tc:                  Time_Context,
	output:              Output_Config,
	linsolve_rtol:       f64,
	linsolve_max_iter:   int,
	non_linear_rtol:     f64,
	non_linear_max_iter: int,
}

Solve_Result :: struct {
	converged:   bool,
	reason:      enum {
		Linear_Solver,
		Non_Linear_Solver,
	},
	final_resid: f64,
	final_iters: int,
}

Model :: struct {
	// arbitray model data
	data:                rawptr,

	// allocate a layout that is sized correctly for the pde to be solved.
	// also return which global dofs are constrained.
	define_layout:       proc(m: ^Model, mesh: fem.Mesh, allocator: mem.Allocator) -> (fem.Layout, fem.Constraint_Mask),

	set_initial_conditions: proc(m: ^Model, layout: fem.Layout, mesh:fem.Mesh, initial_solution: la.Block_Vector),

	// update the current iterate such that it respects the constrained dofs.
	// mask is passed in to make it quicker to identify constrained dofs without duplicating work from `define_layout`.
	// allocator is provided in the case a temporary allocation is needed to compute the constrained value.
	respect_constraints: proc(
		m: ^Model,
		layout: fem.Layout,
		mesh: fem.Mesh,
		mask: fem.Constraint_Mask,
		current_iterate: la.Block_Vector,
		tc: Time_Context,
		allocator: mem.Allocator,
	),

	// layout is provided to READ from current iterate for non linear terms and to WRITE to R and J.
	build_local_problem: proc(
		m: ^Model,
		layout: fem.Layout,
		element: fem.Mesh_Element,
		current_iterate: la.Block_Vector,
		previous_iterate: la.Block_Vector,
		tc: Time_Context,
		R: la.Block_Vector,
		J: la.Block_Dense_Matrix,
		allocator: mem.Allocator,
	),

	// model decides based on its configuration which fields to output.
	// current iterate should be read to actually compute the field data for the global fields.
	// append into output fields.
	output_fields:       proc(
		m: ^Model,
		layout: fem.Layout,
		current_iterate: la.Block_Vector,
		output_fields: ^[dynamic]fio.Output_Field,
	),
}

solve_model :: proc(config: Solver_Config, model: ^Model, mesh: fem.Mesh) -> Solve_Result {
	ca := infra.Checkpoint_Allocator{}

	// this needs to use previous allocator so its not wiped on a reset.
	pvd_paths := make([dynamic]string)
	pvd_times := make([dynamic]f64)
	defer delete(pvd_paths)
	defer delete(pvd_times)

	infra.ca_init(&ca)
	defer infra.ca_deinit(&ca)

	perm_alloc := context.allocator
	context.allocator = infra.ca_allocator(&ca)

	// setup output
	viz_mesh := fio.vtk_create_visualization_mesh(mesh, config.output.output_rule)
	if !os2.exists(config.output.output_dir) {
		os2.make_directory(config.output.output_dir)
	}

	// solver state
	layout, constraint_mask := model->define_layout(mesh, context.allocator)

	u_k := fem.layout_make_vec(&layout) // current solution
	u_prev := fem.layout_make_vec(&layout) //solution at previous time
	R := fem.layout_make_vec(&layout) // working residual vector
	J := fem.layout_make_matrix(&layout) // working jacobian matrix
	du := fem.layout_make_vec(&layout) // incremental update to u



	// time loop
	tc := config.tc

	if tc.is_transient {
		model->set_initial_conditions(layout, mesh, u_prev)
		model->respect_constraints(layout, mesh, constraint_mask, u_prev, tc, context.allocator)
		copy(u_k.values, u_prev.values)
	}else{
		model->set_initial_conditions(layout, mesh, u_k)
	}

	for steps_taken := 0;; steps_taken += 1 {
		loop_start := infra.ca_check(&ca)
		defer infra.ca_rewind_to(&ca, loop_start)

		// constraint values, possibly changer per timestep, but are required to be linear.
		model->respect_constraints(layout, mesh, constraint_mask, u_k, tc, context.allocator)

		// compute initial J and R
		la.scal(R, 0)
		la.scal(J, 0)
		for element in mesh.elements {
			element_loop := infra.ca_check(&ca)
			defer infra.ca_rewind_to(&ca, element_loop)

			local_r, local_j := fem.layout_local_problem(&layout, element)

			model->build_local_problem(layout, element, u_k, u_prev, tc, local_r, local_j, context.allocator)

			fem.layout_scatter_local(layout, element.id, J, R, local_j, local_r, constraint_mask)
		}
		fem.layout_finalize_constraints(layout, J, R, constraint_mask)

		R0_norm := la.nrm2(R)

		// non linear loop.
		// for purely linear problems this is a slightly unoptimal as we build R and J twice.
		non_lin_iters: int
		for {
			non_lin_iters += 1
			converged, iters, resid := la.sparse_cg_solve(J, du, R, config.linsolve_max_iter, config.linsolve_rtol)

			if converged == false {
				return {converged = false, reason = .Linear_Solver, final_resid = resid, final_iters = iters}
			} else {
				log.debugf("Linsolve converged after %d iterations. final residual %f", iters, resid)
			}

			la.axpy(du, u_k, 1)

			la.scal(R, 0)
			la.scal(J, 0)
			for element in mesh.elements {
				element_loop := infra.ca_check(&ca)
				defer infra.ca_rewind_to(&ca, element_loop)

				local_r, local_j := fem.layout_local_problem(&layout, element)

				model->build_local_problem(layout, element, u_k, u_prev, tc, local_r, local_j, context.allocator)

				fem.layout_scatter_local(layout, element.id, J, R, local_j, local_r, constraint_mask)
			}
			fem.layout_finalize_constraints(layout, J, R, constraint_mask)

			R_norm := la.nrm2(R)
			if R_norm / R0_norm < config.non_linear_rtol {
				log.debugf("Solution converged after %d iterations with final residual %f", non_lin_iters, R_norm / R0_norm)
				break
			}
			if non_lin_iters >= config.non_linear_max_iter {
				return {
					converged = false,
					reason = .Non_Linear_Solver,
					final_resid = R_norm / R0_norm,
					final_iters = non_lin_iters,
				}
			}
		}

		copy(u_prev.values, u_k.values)

		tc.current_time += tc.timestep

		if (steps_taken + 1) % config.output.frequency == 0 {
			output_field_buffer := make([dynamic]fio.Output_Field, 0, 16)

			model->output_fields(layout, u_k, &output_field_buffer)

			filename := fmt.aprintf("%s_%d.vtu", config.output.prefix, steps_taken)
			path, _ := os2.join_path({config.output.output_dir, filename}, context.allocator)

			fio.write_vtu(path, mesh, viz_mesh, output_field_buffer[:])

			if tc.is_transient {
				append(&pvd_paths, strings.clone(path, perm_alloc))
				append(&pvd_times, tc.current_time)
			}
		}

		if tc.current_time >= tc.end_time || !tc.is_transient {
			break
		}
	}

	if tc.is_transient {
		filename := fmt.aprintf("%s.pvd", config.output.prefix)
		path, _ := os2.join_path({config.output.output_dir, filename}, context.allocator)
		fio.write_pvd(path, pvd_paths[:], pvd_times[:])
	}

	return {converged = true}
}


//TODO: consider moving this stuff into layout

field_to_output_field :: proc(
	layout: fem.Layout,
	field: fem.Field_Handle,
	u: la.Block_Vector,
	name: string,
) -> fio.Output_Field {
	fd := layout.fds[field]


	// TODO: this would be a mem leak if we werent using checkpoint allocator.
	info := new(Output_Field_Info)

	info.layout = layout.dof_layouts[field]
	info.order = fd.order
	info.data = la.block_vector_view(u, int(field))

	cmpnts: int
	fn: fio.Output_Proc
	switch fd.family {
	case .LS:
		cmpnts = 1
		fn = output_lagrange_scalar
	case .LV:
		cmpnts = 3
		fn = output_lagrange_vector
	}
	return {friendly_name = name, components = cmpnts, value_provider = fn, data = info}
}

field_coeff :: #force_inline proc(
	layout: fem.Layout,
	field: fem.Field_Handle,
	u: la.Block_Vector,
	element_id: fem.Entity_ID,
	dof: int,
) -> f64 {
	return u.values[fem.layout_global_pos(layout, field, element_id, dof)]
}


field_mark_constraints :: proc(
	mesh: fem.Mesh,
	layout: fem.Layout,
	cm: fem.Constraint_Mask,
	fh: fem.Field_Handle,
	constrained_facets: []fem.Boundary_ID,
) {
	for id in constrained_facets {
		for dof in layout.bound_maps[fh][id] {
			global := fem.layout_global_pos(layout, fh, dof.element, dof.local)
			cm.is_constrained[global] = true

		}
	}
}

Output_Field_Info :: struct {
	data:   []f64,
	layout: fem.DOF_Layout,
	order:  fem.Order,
}

evaluate_lagrange_scalar :: proc(
	layout: fem.Layout,
	ctx: fem.Element_Context,
	field: fem.Field_Handle,
	basis: fem.Basis_LS,
	u: la.Block_Vector,
	out: []f64,
) {
	for pi in 0..<len(ctx.points) {
		for dof in 0..<basis.arity {
			coeff := field_coeff(layout, field, u, ctx.element.id, dof)
			out[pi] += fem.basis_value(basis, ctx, pi, dof) * coeff
		}
	}
}

evaluate_lagrange_scalar_gradient :: proc(
	layout: fem.Layout,
	ctx: fem.Element_Context,
	field: fem.Field_Handle,
	basis: fem.Basis_LS,
	u: la.Block_Vector,
	out: []fem.Vec3,
) {
	for pi in 0..<len(ctx.points) {
		for dof in 0..<basis.arity {
			coeff := field_coeff(layout, field, u, ctx.element.id, dof)
			out[pi] += fem.basis_gradient(basis, ctx, pi, dof) * coeff
		}
	}
}

output_lagrange_scalar :: proc(ctx: fem.Element_Context, point: int, data: rawptr) -> (r: [fio.MAX_OUTPUT_FIELD_COMPONENTS]f64) {
	info := cast(^Output_Field_Info)data

	basis := fem.basis_create(fem.Basis_LS, ctx.element, info.order)

	for i in 0 ..< basis.arity {
		r[0] += fem.basis_value(basis, ctx, point, i) * info.data[info.layout.mapping[ctx.element.id][i]]
	}

	return r
}


output_lagrange_vector :: proc(ctx: fem.Element_Context, point: int, data: rawptr) -> (r: [fio.MAX_OUTPUT_FIELD_COMPONENTS]f64) {
	info := cast(^Output_Field_Info)data

	basis := fem.basis_create(fem.Basis_LV, ctx.element, info.order)

	for i in 0 ..< basis.arity {
		v := fem.basis_value(basis, ctx, point, i) * info.data[info.layout.mapping[ctx.element.id][i]]
		r[0] += v[0]; r[1] += v[1]; r[2] += v[2]
	}

	return r
}
