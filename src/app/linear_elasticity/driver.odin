package linear_elasticity

import "core:fmt"
import "core:log"
import "core:mem"
import "core:mem/virtual"
import "core:os"
import "core:slice"

import "../cfg"

import "../../fem"
import "../../fem/infra"
import "../../fem/serialization"


run_simulation :: proc(
	gcfg: cfg.General_Params,
	outputter: ^cfg.Outputter,
	model_cfg: Model_Config,
	arena: ^virtual.Arena,
	prt: ^infra.Parallel_Runtime,
) -> bool {
	arena_allocator := virtual.arena_allocator(arena)
	context.allocator = arena_allocator

	log.info("Setting up state...")

	sys_desc := fem.System_Description{}

	handle := fem.description_add_variable(&sys_desc, {bd = model_cfg.bd, components = 3})

	fem.description_couple(&sys_desc, handle, handle)

	system := fem.system_from_description(gcfg.mesh, sys_desc, context.allocator)

	cm := fem.system_constraint_mask(system)
	ics := fem.system_vector(system)

	transient_cfg, is_transient := gcfg.transient.?

	apply_ics(model_cfg.ics, system, handle, transient_cfg.start, gcfg.mesh, ics)
	apply_constraints(model_cfg.params, system, handle, gcfg.mesh, transient_cfg.start, ics, cm)

	ni_state := fem.nli_create_state(system, ics)

	ts: fem.Timestepper; ts_state: fem.Timestep_State
	if is_transient {
		ts, ts_state = fem.timestepper_create(transient_cfg.start, transient_cfg.end, transient_cfg.timestep, ics)
		fem.timestepper_set_scheme(&ts, &ts_state, handle, .Newmark)
	} else {
		ts, ts_state = fem.timestepper_create_steady(ics, system)
	}

	thread_partitions := fem.system_create_thread_partitions(system, prt, context.allocator)

	parallel_data: Parallel_Data
	parallel_data.mesh = gcfg.mesh
	parallel_data.partitions = thread_partitions
	parallel_data.handle = handle
	parallel_data.ni_state = ni_state
	parallel_data.system = system
	parallel_data.cm = cm
	parallel_data.params = model_cfg.params
	parallel_data.allocators = make([]infra.Checkpoint_Allocator, prt.total_threads)
	for &alloc in parallel_data.allocators {infra.ca_init(&alloc)}
	defer {
		for &alloc in parallel_data.allocators {infra.ca_deinit(&alloc)}
	}

	out_data: serialization.Output_Variable_Data = {system, handle, ni_state.solution}
	out_field := serialization.output_field_from_system_variable(fem.Grad_Space(.Vector), &out_data, "displacement")

	von_mises_data := Von_Mises_Data {
		system       = system,
		var          = handle,
		displacement = ni_state.solution,
		materials    = model_cfg.params.materials,
	}

	von_field := von_mises_output(&von_mises_data)

	solve_opt := fem.Solver_Options {
		kind       = .CG_SA,
		tol        = gcfg.solver.tolerance,
		max_iters  = gcfg.solver.max_linear_iters,
		block_size = 3,
	}

	log.info("Simulation starting...")

	ca: infra.Checkpoint_Allocator
	infra.ca_init(&ca)
	context.allocator = infra.ca_allocator(&ca)
	defer infra.ca_deinit(&ca)

	cfg.output_step(outputter, gcfg.mesh, fem.timestepper_initial_step(ts), {out_field, von_field}, true)

	for step in fem.timestepper_step(&ts, ts_state, ni_state.solution, system) {
		infra.ca_check(&ca); defer infra.ca_rewind(&ca)
		apply_constraints(model_cfg.params, system, handle, gcfg.mesh, step.time, ni_state.solution, cm)

		parallel_data.step = step

		slice.zero(ni_state.residual)
		slice.zero(ni_state.tangent.values)

		u_dot, u_ddot := fem.timestepper_derivatives(&ts, ts_state, ni_state.solution, system)
		parallel_data.u_dot = u_dot
		parallel_data.u_ddot = u_ddot

		infra.parallel_for(prt, {0, len(gcfg.mesh.elements)}, parallel_data, assembly_proc)
		fem.system_flush_orphans(parallel_data.partitions, ni_state.tangent, ni_state.residual)
		fem.system_finalize_constraints(system, ni_state.tangent, ni_state.residual, cm)

		if res, ok := fem.sparse_solve(ni_state.tangent, ni_state.update, ni_state.residual, solve_opt); !ok {
			log.errorf(
				"Linear solver %v failed to converge in %d iterations with final residual %f",
				solve_opt.kind,
				res.iters,
				res.residual,
			)
			return false
		}

		fem.axpy(ni_state.update, ni_state.solution, 1.0)

		cfg.output_step(outputter, gcfg.mesh, step, {out_field, von_field})
	}

	cfg.output_flush(outputter)

	log.info("Simulation complete.")
	return true
}

Parallel_Data :: struct {
	allocators: []infra.Checkpoint_Allocator,
	partitions: fem.Thread_Partitions,
	mesh:       fem.Mesh,
	params:     Model_Parameters,
	ni_state:   fem.Nonlinear_Iter_State,
	handle:     fem.Var_Handle,
	u_dot:      fem.Vector,
	u_ddot:     fem.Vector,
	step:       fem.Timestep,
	cm:         []bool,
	system:     fem.System,
}

assembly_proc :: proc(data: Parallel_Data, range: infra.Range) {
	ca := &data.allocators[infra.tid]
	infra.ca_check(ca); defer infra.ca_rewind(ca)
	context.allocator = infra.ca_allocator(ca)

	for i in range.start ..< range.end {
		element_id := fem.Entity_ID(i)

		ls := fem.system_local_problem(data.system, element_id)

		weak_form(
			data.params,
			data.system,
			ls,
			data.handle,
			data.ni_state.solution,
			data.u_dot,
			data.u_ddot,
			data.step,
			data.mesh.elements[element_id],
		)

		fem.system_thread_scatter(
			data.system,
			element_id,
			data.ni_state.tangent,
			data.ni_state.residual,
			data.cm,
			ls,
			data.partitions.partitions[infra.tid],
		)
	}

}
