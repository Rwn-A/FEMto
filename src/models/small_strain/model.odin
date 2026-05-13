package small_strain

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
    handle: fem.Var_Handle,
    solution: fem.Vector,
}

Solver_Options :: struct {
    using linsolve_options: fem.Solver_Options,
    max_nonlinear_iters: int,
    tolerance: f64,
}

Small_Strain_Driver :: struct {
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

configure_driver :: proc(mesh: fem.Mesh, out_h: ^cfg.Output_Handler, plugins: ^cfg.Plugin_Registry, schema: cfg.Root_Schema, allocator := context.allocator) -> (cd: Small_Strain_Driver, ok: bool) {
	SUPPORTED_FAMILIES :: bit_set[fem.Basis_Family]{.Lagrange}
	SUPPORTED_TIME_SCHEMES :: bit_set[fem.Time_Scheme]{.Newmark, .Steady}
	SUPPORTED_LINEAR_SOLVERS: bit_set[fem.Solver_Kind]  : ~{} // anything goes

    context.allocator = allocator

    system_desc: fem.System_Description

    bd := cfg.load_discretization(SUPPORTED_FAMILIES, .Lagrange, schema.discretization) or_return
	cd.handle = fem.description_add_variable(&system_desc, {bd = bd, components = 3})
	fem.description_couple(&system_desc, cd.handle, cd.handle)
	cd.sys = fem.system_from_description(mesh, system_desc, context.allocator)

	if tc := schema.time_control; tc != {}{
	    cd.ts = fem.timestepper_create(tc.start, tc.end, tc.timestep)
	    time_scheme := cfg.load_time_scheme(SUPPORTED_TIME_SCHEMES, .Newmark, tc.scheme) or_return
		fem.timestepper_set_scheme(&cd.ts, cd.handle, time_scheme)
	} else{
	   cd.ts = fem.timestepper_create_steady()
	}

	cd.solution = fem.system_vector(cd.sys)

	displ_data :=  new(fio.Output_Variable_Data, out_h.allocator)
    displ_data^ = {cd.sys, cd.handle, cd.solution}
	displ_field := fio.output_field_from_system_variable(fem.Grad_Space(.Vector), displ_data, "displacement")
	cfg.output_handler_add_field(out_h, displ_field)

	cd.solver_options.max_nonlinear_iters = schema.solver.max_nonlinear_iters
	cd.solver_options.tolerance = schema.solver.tolerance
	cd.solver_options.linsolve_options = {
        max_iters = schema.solver.max_linear_iters,
        kind = cfg.load_solver_kind(SUPPORTED_LINEAR_SOLVERS, .CG_SA, schema.solver.override) or_return,
        block_size = 3,
        tol = schema.solver.tolerance
	}

	cd.out_h = out_h

    cd.problem = empty_problem_data(mesh)

    ctx := cfg.Plugin_Context {problem_data = &cd.problem}

    builtin_register(plugins)
    defer builtin_deregister(plugins)

    if bcs, has := cfg.get_field_table(schema.boundaries, "small_strain") or_return; has {
        cfg.load_bc_plugins(bcs, mesh, plugins^, &ctx) or_return
    }else{
        log.errorf("Atleast one boundary must be defined.")
        return {}, false
    }

    if sources, has := cfg.get_field_table(schema.sources, "small_strain") or_return; has {
        cfg.load_section_plugins(sources, mesh, plugins^, &ctx, .Source) or_return
    }

    if mats, has := cfg.get_field_table(schema.materials, "small_strain") or_return; has {
        cfg.load_section_plugins(mats, mesh, plugins^, &ctx, .Material) or_return
    }else{
        log.errorf("Atleast one material must be defined.")
        return {}, false
    }

    if ics, has := cfg.get_field_table(schema.initial_conditions, "small_strain") or_return; has {
        cfg.load_section_plugins(ics, mesh, plugins^, &ctx, .IC) or_return
    }

    von_mises_data := new(Von_Mises_Data, out_h.allocator)
	von_mises_data^ = Von_Mises_Data{
		system       = cd.sys,
		var          = cd.handle,
		displacement = cd.solution,
		materials    = cd.problem.materials,
	}

	cfg.output_handler_add_field(out_h, von_mises_output(von_mises_data))


    return cd, true
}


drive :: proc(cd: Small_Strain_Driver, mesh: fem.Mesh, prt: ^infra.Parallel_Runtime, allocator := context.allocator) -> bool{
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


    cfg.output_handler_step(cd.out_h, mesh, fem.timestepper_initial_step(cd.ts))

    for step in fem.timestepper_step(&cd.ts, cd.sys, cd.solution, &ad.step_state) {
        infra.ca_check(&ca); defer infra.ca_rewind(&ca)
        apply_constraints(cd.sys_ctx, mesh, cd.ts.start, cd.solution, ad.cm)

        fem.timestepper_update_derivatives(&cd.ts, &ad.step_state, cd.sys, cd.solution)

        infra.parallel_for(prt, {0, len(mesh.elements)}, ad, assembly_proc)

		fem.system_flush_orphans(ad.partitions, ad.tangent, ad.residual)
		fem.system_finalize_constraints(cd.sys, ad.tangent, ad.residual, ad.cm)

        fem.sparse_solve(ad.tangent, ad.update, ad.residual, cd.solver_options)

        fem.axpy(ad.update, cd.solution, 1.0)
        slice.zero(ad.residual)
        slice.zero(ad.tangent.values)

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
	bd        := fem.system_var_bd(sys_ctx.sys, sys_ctx.handle)
	quad_rule := fem.infer_quadrature(bd.order)

	quad  := fem.map_quadrature(element, {.Interior, quad_rule, 0})
	space := fem.basis_grad_space(quad, bd, .Vector)

	interior := bulk_response(sys_ctx, step_state, quad)
	assemble_bulk(ls, sys_ctx.handle, space, quad, step_state.step, interior)

	for facet in element.boundary_facets {
		quad  := fem.map_quadrature(element, {.Surface, quad_rule, facet})
		space := fem.basis_grad_space(quad, bd, .Vector)
		id    := element.boundary_ids[facet]

		if id not_in sys_ctx.problem.variational_bcs {continue}

		boundary := boundary_response(sys_ctx, step_state, quad, id)
		assemble_boundary(ls, sys_ctx.handle, space, quad, boundary)
	}

	bulk_response :: proc(sys_ctx: System_Context, time_ctx: fem.Timestep_State, mapped: fem.Mapped_Element) -> Bulk_Response {
	    br: Bulk_Response

	    space := fem.basis_grad_space(mapped, fem.system_var_bd(sys_ctx.sys, sys_ctx.handle), .Vector)

		br.u_coeffs = fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.handle, mapped.element.id, sys_ctx.solution)
		br.u_dot_coeffs = fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.handle, mapped.element.id, time_ctx.ddt)
		br.u_ddot_coeffs = fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.handle, mapped.element.id, time_ctx.d2dt2)
		br.material = empty_material(fem.space_points(space))
		br.source  = empty_source(fem.space_points(space))

		current_strain := fem.evaluate_var_symmetric_gradient(space, mapped, br.u_coeffs)
    	for qp in 0 ..< fem.space_points(space) {
    		current_strain[qp] = to_engineering_shear(current_strain[qp])
    	}


		mat := sys_ctx.problem.materials[mapped.element.section]
		mat.procedure(mapped, time_ctx.step.time, current_strain, mat.data, br.material)


		for source_int in sys_ctx.problem.sources[mapped.element.section] {
			source_int.procedure(mapped, time_ctx.step.time, source_int.data, br.source)
		}

		return br
	}

	boundary_response :: proc(sys_ctx: System_Context, time_ctx: fem.Timestep_State, mapped: fem.Mapped_Element, id: fem.Boundary_ID ) -> Boundary_Response {
	    br: Boundary_Response

	    space := fem.basis_grad_space(mapped, fem.system_var_bd(sys_ctx.sys, sys_ctx.handle), .Vector)

		u_coeffs := fem.system_gather_var_coeffs(sys_ctx.sys, sys_ctx.handle, mapped.element.id, sys_ctx.solution)

		br.displacement = fem.evaluate_var(space, u_coeffs)

		br.variational = empty_variational(fem.space_points(space))
		for var_int in sys_ctx.problem.variational_bcs[id] {
		   var_int.procedure(mapped, time_ctx.step.time, var_int.data, br.variational)
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
		iso := sys_ctx.problem.fixed_bcs[id] or_continue
		fem.system_apply_boundary_functional(sys_ctx.sys, sys_ctx.handle, mesh, id, time, vec, cm, iso.procedure, iso.data)
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
	for section_id, ic in sys_ctx.problem.ics {
		fem.system_project_dofs(sys_ctx.sys, sys_ctx.handle, mesh, section_id, start_time, vec, ic.procedure, ic.data)
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
