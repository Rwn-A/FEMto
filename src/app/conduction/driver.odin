package conduction

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

	handle := fem.description_add_variable(&sys_desc, {bd = model_cfg.bd, components = 1})

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
		fem.timestepper_set_scheme(&ts, &ts_state, handle, model_cfg.time_scheme)
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
	out_field := serialization.output_field_from_system_variable(fem.Grad_Space(.Scalar), &out_data, "temperature")

	solve_opt := fem.Solver_Options {
		kind       = .CG_SA,
		max_iters  = gcfg.solver.max_linear_iters,
		block_size = 1,
	}

	log.info("Simulation starting...")

	ca: infra.Checkpoint_Allocator
	infra.ca_init(&ca)
	context.allocator = infra.ca_allocator(&ca)
	defer infra.ca_deinit(&ca)

	cfg.output_step(outputter, gcfg.mesh, fem.timestepper_initial_step(ts), {out_field}, true)


	for step in fem.timestepper_step(&ts, ts_state, ni_state.solution, system) {
		infra.ca_check(&ca); defer infra.ca_rewind(&ca)
		apply_constraints(model_cfg.params, system, handle, gcfg.mesh, step.time, ni_state.solution, cm)

		parallel_data.step = step

		ni := fem.nli_create(gcfg.solver.tolerance, gcfg.solver.max_nonlinear_iters)
		solve_opt.tol = fem.nli_linear_tolerance(ni)

		for iter, should_solve in fem.nli_step(&ni, ni_state) {
			infra.ca_check(&ca); defer infra.ca_rewind(&ca)

			if should_solve {
				if res, ok := fem.sparse_solve(
					ni_state.tangent,
					ni_state.update,
					ni_state.residual,
					solve_opt,
				); !ok {
					log.errorf(
						"Linear solver %v failed to converge in %d iterations with final residual %f",
						solve_opt.kind,
						res.iters,
						res.residual,
					)
					return false
				}
				fem.nli_update(&ni, ni_state)
			}

			u_dot, _ := fem.timestepper_derivatives(&ts, ts_state, ni_state.solution, system)
			parallel_data.u_dot = u_dot

			infra.parallel_for(prt, {0, len(gcfg.mesh.elements)}, parallel_data, assembly_proc)

			fem.system_flush_orphans(parallel_data.partitions, ni_state.tangent, ni_state.residual)

			fem.system_finalize_constraints(system, ni_state.tangent, ni_state.residual, cm)
		}

		if fem.nli_reached_max(&ni) {
			log.errorf("Non linear system did not converge to a result in %d iterations.", ni.current_iter)
			if is_transient {
				log.infof("Current timestep: %d", ts.current_step)
			}
			return false
		}

		cfg.output_step(outputter, gcfg.mesh, step, {out_field})
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
