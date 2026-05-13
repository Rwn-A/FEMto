package inc_flow

import "core:log"
import "core:slice"
import "core:math"

import "../../cfg"
import "../../fem"
import "../../fem/infra"
import fio"../../fem/serialization"


System_Context :: struct {
    problem: Problem_Data,
    sys: fem.System,
    v_handle: fem.Var_Handle,
    p_handle: fem.Var_Handle,
    solution: fem.Vector,
}

Solver_Options :: struct {
    using linsolve_options: fem.Solver_Options,
    max_nonlinear_iters: int,
    tolerance: f64,
}

Inc_Flow_Driver :: struct {
	using sys_ctx: System_Context,
    ts: fem.Timestepper,
    out_h: ^cfg.Output_Handler,
    solver_options: Solver_Options,
}

Assembly_Data :: struct {
	allocators: []infra.Checkpoint_Allocator,
	partitions: fem.Thread_Partitions,
	sys_ctx: System_Context,
	step_state: fem.Timestep_State,
	mesh: fem.Mesh,
	residual, update: fem.Vector,
	tangent: fem.Sparse_Matrix,
	cm: []bool,
}

configure_driver :: proc(mesh: fem.Mesh, out_h: ^cfg.Output_Handler, plugins: ^cfg.Plugin_Registry, schema: cfg.Root_Schema, allocator := context.allocator) -> (cd: Inc_Flow_Driver, ok: bool) {
	SUPPORTED_FAMILIES :: bit_set[fem.Basis_Family]{.Lagrange}
	SUPPORTED_TIME_SCHEMES :: bit_set[fem.Time_Scheme]{.BE, .BDF2}
	SUPPORTED_LINEAR_SOLVERS: bit_set[fem.Solver_Kind]  : ~{} // anything goes

    context.allocator = allocator

    system_desc: fem.System_Description

    bd := cfg.load_discretization(SUPPORTED_FAMILIES, .Lagrange, schema.discretization) or_return
	cd.v_handle = fem.description_add_variable(&system_desc, {bd = bd, components = 3})
	cd.p_handle = fem.description_add_variable(&system_desc, {bd = bd, components = 1})

	fem.description_couple(&system_desc, cd.v_handle, cd.v_handle)
	fem.description_couple(&system_desc, cd.v_handle, cd.p_handle)
	fem.description_couple(&system_desc, cd.p_handle, cd.v_handle)
	fem.description_couple(&system_desc, cd.p_handle, cd.p_handle)


	cd.sys = fem.system_from_description(mesh, system_desc, context.allocator)

	if tc := schema.time_control; tc != {}{
	    cd.ts = fem.timestepper_create(tc.start, tc.end, tc.timestep)
	    time_scheme := cfg.load_time_scheme(SUPPORTED_TIME_SCHEMES, .Newmark, tc.scheme) or_return
		fem.timestepper_set_scheme(&cd.ts, cd.v_handle, time_scheme)
	} else{
	   cd.ts = fem.timestepper_create_steady()
	}

	cd.solution = fem.system_vector(cd.sys)

	velocity_data :=  new(fio.Output_Variable_Data, out_h.allocator)
    velocity_data^ = {cd.sys, cd.v_handle, cd.solution}
	velocity_field := fio.output_field_from_system_variable(fem.Grad_Space(.Vector), velocity_data, "velocity")
	cfg.output_handler_add_field(out_h, velocity_field)

	pressure_data :=  new(fio.Output_Variable_Data, out_h.allocator)
    pressure_data^ = {cd.sys, cd.p_handle, cd.solution}
	pressure_field := fio.output_field_from_system_variable(fem.Grad_Space(.Scalar), pressure_data, "pressure")
	cfg.output_handler_add_field(out_h, pressure_field)


	cd.solver_options.max_nonlinear_iters = schema.solver.max_nonlinear_iters
	cd.solver_options.tolerance = schema.solver.tolerance
	cd.solver_options.linsolve_options = {
        max_iters = schema.solver.max_linear_iters,
        kind = cfg.load_solver_kind(SUPPORTED_LINEAR_SOLVERS, .FGMRES_ILU0, schema.solver.override) or_return,
        block_size = 1,
        tol = schema.solver.tolerance
	}

	cd.out_h = out_h

    cd.problem = empty_problem_data(mesh)

    ctx := cfg.Plugin_Context {problem_data = &cd.problem}

    builtin_register(plugins)
    defer builtin_deregister(plugins)

    if bcs, has := cfg.get_field_table(schema.boundaries, "momentum") or_return; has {
        cfg.load_bc_plugins(bcs, mesh, plugins^, &ctx) or_return
    }else{
        log.errorf("Atleast one boundary must be defined.")
        return {}, false
    }

    if sources, has := cfg.get_field_table(schema.sources, "momentum") or_return; has {
        cfg.load_section_plugins(sources, mesh, plugins^, &ctx, .Source) or_return
    }

    if mats, has := cfg.get_field_table(schema.materials, "momentum") or_return; has {
        cfg.load_section_plugins(mats, mesh, plugins^, &ctx, .Material) or_return
    }else{
        log.errorf("Atleast one material must be defined.")
        return {}, false
    }

    if ics, has := cfg.get_field_table(schema.initial_conditions, "velocity") or_return; has {
        cfg.load_section_plugins(ics, mesh, plugins^, &ctx, .IC) or_return
    }

    if ics, has := cfg.get_field_table(schema.initial_conditions, "pressure") or_return; has {
        cfg.load_section_plugins(ics, mesh, plugins^, &ctx, .IC) or_return
    }

    return cd, true
}


drive :: proc(cd: Inc_Flow_Driver, mesh: fem.Mesh, prt: ^infra.Parallel_Runtime, allocator := context.allocator) -> bool{
    context.allocator = allocator
    cd := cd

    log.info("Simulation starting...")

	ics := fem.system_vector(cd.sys)

    ad: Assembly_Data = {
		cm = fem.system_constraint_mask(cd.sys),
		residual = fem.system_vector(cd.sys),
		allocators = make([]infra.Checkpoint_Allocator, prt.total_threads),
		tangent = fem.system_matrix(cd.sys),
		sys_ctx = cd.sys_ctx,
		update = fem.system_vector(cd.sys),
		mesh = mesh,
		partitions = fem.system_create_thread_partitions(cd.sys, prt, allocator),
		step_state = {
			ddt = fem.system_vector(cd.sys),
			d2dt2 = fem.system_vector(cd.sys),
		}
    }
    for &alloc in ad.allocators {infra.ca_init(&alloc)}
	defer { for &alloc in ad.allocators {infra.ca_deinit(&alloc)} }


    // load ics
	apply_ics(cd.sys_ctx, mesh, cd.ts.start, ics)
	apply_constraints(cd.sys_ctx, mesh, cd.ts.start, ics, ad.cm)

	copy(cd.solution, ics)

    // setup timestepper state
    fem.timestepper_init_history(&cd.ts, cd.sys, ics)

    ca: infra.Checkpoint_Allocator
	infra.ca_init(&ca)
	context.allocator = infra.ca_allocator(&ca)
	defer infra.ca_deinit(&ca)


    cfg.output_handler_step(cd.out_h, mesh, fem.timestepper_initial_step(cd.ts), force = true)

    for step in fem.timestepper_step(&cd.ts, cd.sys, cd.solution, &ad.step_state) {
        infra.ca_check(&ca); defer infra.ca_rewind(&ca)
        apply_constraints(cd.sys_ctx, mesh, cd.ts.start, cd.solution, ad.cm)

        ni := fem.nli_create(cd.solver_options.tolerance, cd.solver_options.max_nonlinear_iters)
        cd.solver_options.linsolve_options.tol = fem.nli_linear_tolerance(ni)

        for should_solve in fem.nli_step(&ni, cd.solution, ad.tangent, ad.residual, ad.update) {
        	infra.ca_check(&ca); defer infra.ca_rewind(&ca)
        	if should_solve {
				fem.sparse_solve(ad.tangent, ad.update, ad.residual, cd.solver_options)
				fem.nli_do_update(cd.solution, ad.tangent, residual = ad.residual, update = ad.update)
			}

			fem.timestepper_update_derivatives(&cd.ts, &ad.step_state, cd.sys, cd.solution)


			infra.parallel_for(prt, {0, len(mesh.elements)}, ad, assembly_proc)

			fem.system_flush_orphans(ad.partitions, ad.tangent, ad.residual)
			fem.system_finalize_constraints(cd.sys, ad.tangent, ad.residual, ad.cm)
		}

		if fem.nli_reached_max(ni) {
			log.errorf("Non linear system did not converge to a result in %d iterations.", ni.current_iter)
			//return false
		}

		// -2 because it already incremented past when it was done, then because the first iter doenst solve.
		log.debugf("Nonlinear iterations complete after %d iterations", ni.current_iter - 2)

        cfg.output_handler_step(cd.out_h, mesh, step)
    }

    cfg.output_handler_flush_pvd(cd.out_h)

    log.info("Simulation Complete!")

    return true
}

assembly_proc :: proc(data: Assembly_Data, range: infra.Range) {
	ca := &data.allocators[infra.tid]
	infra.ca_check(ca); defer infra.ca_rewind(ca)
	context.allocator = infra.ca_allocator(ca)

	for i in range.start ..< range.end {
		element := data.mesh.elements[i]
		partition := data.partitions.partitions[infra.tid]

		ls := fem.system_local_problem(data.sys_ctx.sys, element.id)

		weak_form(data.sys_ctx, data.step_state, ls, element)

		fem.system_thread_scatter(data.sys_ctx.sys, element.id, data.tangent, data.residual, data.cm, ls, partition)
	}
}


weak_form :: proc(
	sys_ctx: System_Context,
	step_state: fem.Timestep_State,
	ls:       fem.Local_System,
	element:  fem.Mesh_Element,
) {
	vbd        := fem.system_var_bd(sys_ctx.sys, sys_ctx.v_handle)
	pbd        := fem.system_var_bd(sys_ctx.sys, sys_ctx.p_handle)
	quad_rule := fem.infer_quadrature(vbd.order)

	quad  := fem.map_quadrature(element, {.Interior, quad_rule, 0})
	v_space := fem.basis_grad_space(quad, vbd, .Vector)
	p_space := fem.basis_grad_space(quad, pbd, .Scalar)

	interior := bulk_response(sys_ctx, step_state, quad)
	assemble_bulk(ls, sys_ctx.v_handle, sys_ctx.p_handle, v_space, p_space, quad, step_state.step, interior)

	for facet in element.boundary_facets {
		quad  := fem.map_quadrature(element, {.Surface, quad_rule, facet})
		v_space := fem.basis_grad_space(quad, vbd, .Vector)
		p_space := fem.basis_grad_space(quad, pbd, .Scalar)
		id    := element.boundary_ids[facet]

		if id not_in sys_ctx.problem.v_var_bcs || id not_in sys_ctx.problem.p_var_bcs {continue}

		boundary := boundary_response(sys_ctx, step_state, quad, id)
		assemble_boundary(ls, sys_ctx.v_handle, sys_ctx.p_handle, v_space, p_space, quad, boundary)
	}

	bulk_response :: proc(sys_ctx: System_Context, time_ctx: fem.Timestep_State, mapped: fem.Mapped_Element) -> Bulk_Response {
	    br: Bulk_Response

        v_space := fem.basis_grad_space(mapped, fem.system_var_bd(sys_ctx.sys, sys_ctx.v_handle), .Vector)
        p_space := fem.basis_grad_space(mapped, fem.system_var_bd(sys_ctx.sys, sys_ctx.p_handle), .Scalar)

        br.v_coeffs = fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.v_handle, mapped.element.id, sys_ctx.solution)
        br.p_coeffs = fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.p_handle, mapped.element.id, sys_ctx.solution)

        br.velocity = fem.evaluate_var(v_space, br.v_coeffs)
        br.grad_v = fem.evaluate_var_gradient(v_space, mapped, br.v_coeffs)
        br.grad_p = fem.evaluate_var_gradient(p_space, mapped, br.p_coeffs)


        br.material = empty_material(fem.space_points(v_space))
        br.source  = empty_source(fem.space_points(v_space))

        mat := sys_ctx.problem.materials[mapped.element.section]
        mat.procedure(mapped, time_ctx.step.time, mat.data, br.material)


        for source_int in sys_ctx.problem.sources[mapped.element.section] {
        	source_int.procedure(mapped, time_ctx.step.time, source_int.data, br.source)
        }

		return br
	}

	boundary_response :: proc(sys_ctx: System_Context, time_ctx: fem.Timestep_State, mapped: fem.Mapped_Element, id: fem.Boundary_ID ) -> Boundary_Response {
	    br: Boundary_Response

	     v_space := fem.basis_grad_space(mapped, fem.system_var_bd(sys_ctx.sys, sys_ctx.v_handle), .Vector)

        br.v_variational = empty_v_variational(fem.space_points(v_space))
        for var_int in sys_ctx.problem.v_var_bcs[id] {
           var_int.procedure(mapped, time_ctx.step.time, var_int.data, br.v_variational)
        }

        br.p_variational = empty_p_variational(fem.space_points(v_space))
        for var_int in sys_ctx.problem.p_var_bcs[id] {
           var_int.procedure(mapped, time_ctx.step.time, var_int.data, br.p_variational)
        }

		return br
	}
}


apply_constraints :: proc(
    sys_ctx: System_Context,
	mesh: fem.Mesh,
	time: f64,
	vec: fem.Vector,
	cm: []bool,
	allocator := context.allocator,
) {
	context.allocator = allocator
	for _, id in mesh.boundary_names {
		fixed_v := sys_ctx.problem.v_fixed_bcs[id] or_continue
		fem.system_apply_boundary_functional(sys_ctx.sys, sys_ctx.v_handle, mesh, id, time, vec, cm, fixed_v.procedure, fixed_v.data)
	}
	for name, id in mesh.boundary_names {
		fixed_p := sys_ctx.problem.p_fixed_bcs[id] or_continue
		fem.system_apply_boundary_functional(sys_ctx.sys, sys_ctx.p_handle, mesh, id, time, vec, cm, fixed_p.procedure, fixed_p.data)
	}
}

apply_ics :: proc(
    sys_ctx: System_Context,
	mesh: fem.Mesh,
	start_time: f64,
	vec: fem.Vector,
	allocator := context.allocator,
) {
	context.allocator = allocator
	for section_id, ic in sys_ctx.problem.v_ics {
		fem.system_project_dofs(sys_ctx.sys, sys_ctx.v_handle, mesh, section_id, start_time, vec, ic.procedure, ic.data)
	}
	for section_id, ic in sys_ctx.problem.p_ics {
		fem.system_project_dofs(sys_ctx.sys, sys_ctx.p_handle, mesh, section_id, start_time, vec, ic.procedure, ic.data)
	}
}