package conduction

import "core:log"
import "core:mem/virtual"
import "core:slice"
import "core:fmt"
import "core:os"
import "core:mem"

import "../config"

import "../../fem"
import "../../fem/infra"
import "../../fem/serialization"


run_simulation :: proc(
	cfg: config.General_Config,
	model_schema: config.Model_Schema,
	arena: ^virtual.Arena,
	prt: ^infra.Parallel_Runtime,
) -> bool {
    arena_allocator := virtual.arena_allocator(arena)
	context.allocator = arena_allocator

	log.info("Loading parameters for conduction model.")

	model_params, bd, time_scheme := parse_conduction_schema(cfg.mesh, model_schema) or_return

	// any solver will work fine for conduction, CG is fastest and will likely work for most problems
	linsolve_kind := config.assign_default(model_schema.linear_solver, fem.Solver_Kind.CG_SA)

	log.info("Setting up state...")

	pvd_paths := make([dynamic]string)
	pvd_times := make([dynamic]f64)

	sys_desc := fem.System_Description{}

	handle := fem.description_add_variable(&sys_desc, {bd = bd, components = 1})

	fem.description_couple(&sys_desc, handle, handle)

	system := fem.system_from_description(cfg.mesh, sys_desc, context.allocator)

	cm := fem.system_constraint_mask(system)
	ics := fem.system_vector(system)

	load_ics(model_params, cfg.mesh, system, handle, ics)
	apply_constraints(model_params, system, handle, cfg.mesh, 0, ics, cm)

	ni_state := fem.nli_create_state(system, ics)

	ts: fem.Timestepper; ts_state: fem.Timestep_State
	transient_cfg, is_transient := cfg.transient.?
	if is_transient {
		ts, ts_state = fem.timestepper_create(transient_cfg.start, transient_cfg.end, transient_cfg.timestep, ics)
		fem.timestepper_set_scheme(&ts, &ts_state, handle, time_scheme)
	} else {
		ts, ts_state = fem.timestepper_create_steady(ics, system)
	}

	thread_partitions := fem.system_create_thread_partitions(system, prt, context.allocator)

	parallel_data: Parallel_Data
	parallel_data.mesh = cfg.mesh
	parallel_data.handle = handle
	parallel_data.ni_state = ni_state
	parallel_data.system = system
	parallel_data.cm = cm
	parallel_data.params = model_params
	parallel_data.allocators = make([]infra.Checkpoint_Allocator, prt.total_threads)
	for &alloc in parallel_data.allocators {infra.ca_init(&alloc)}

	out_data: serialization.Output_Variable_Data = {system, handle, ni_state.solution}
	out_field := serialization.output_field_from_system_variable(fem.Grad_Space(.Scalar), &out_data, "temperature")

	log.info("Simulation starting...")

	ca: infra.Checkpoint_Allocator
	infra.ca_init(&ca)
	context.allocator = infra.ca_allocator(&ca)

    ic_output_path := output_path(.VTU, cfg.output.directory, cfg.sim_name, 0, arena_allocator)
	serialization.write_vtu(ic_output_path, cfg.mesh, cfg.viz_mesh, {out_field})

	if is_transient {
	   append(&pvd_paths, ic_output_path)
	   append(&pvd_times, 0)
	}


	for step in fem.timestepper_step(&ts, ts_state, ni_state.solution, system) {
		infra.ca_check(&ca); defer infra.ca_rewind(&ca)
		apply_constraints(model_params, system, handle, cfg.mesh, step.time, ni_state.solution, cm)

		parallel_data.step = step

		ni := fem.nli_create(cfg.solver.tolerance, cfg.solver.max_nonlinear_iters)
		linear_tol := fem.nli_linear_tolerance(ni)

		for iter, should_solve in fem.nli_step(&ni, ni_state) {
			infra.ca_check(&ca); defer infra.ca_rewind(&ca)

			if should_solve {
				if !fem.sparse_solve(
					ni_state.tangent,
					ni_state.update,
					ni_state.residual,
					tol = linear_tol,
					kind = linsolve_kind,
				) {
				    log.info("Simulation failed.")
				    return false
				}
				fem.nli_update(&ni, ni_state)
			}

			u_dot, _ := fem.timestepper_derivatives(&ts, ts_state, ni_state.solution, system)
			parallel_data.u_dot = u_dot

			infra.parallel_for(prt, {0, len(cfg.mesh.elements)}, parallel_data, assembly_proc)


			fem.system_finalize_constraints(system, ni_state.tangent, ni_state.residual, cm)
		}

        path := output_path(.VTU, cfg.output.directory, cfg.sim_name, step.step + 1, arena_allocator)
		serialization.write_vtu(path, cfg.mesh, cfg.viz_mesh, {out_field})

        if is_transient {
            append(&pvd_paths, path)
             //current_time is the time after this iter is done so its what we want
            append(&pvd_times, ts.current_time)
        }
	}

	if is_transient {
	    path := output_path(.PVD, cfg.output.directory, cfg.sim_name, allocator = arena_allocator)
		serialization.write_pvd(path, pvd_paths[:], pvd_times[:])
	}

    log.info("Simulation complete.")
    return true
}

Parallel_Data :: struct {
	allocators: []infra.Checkpoint_Allocator,
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

		fem.system_scatter(data.system, element_id, data.ni_state.tangent, data.ni_state.residual, data.cm, ls)
	}

}

parse_conduction_schema :: proc(
    mesh: fem.Mesh,
	schema: config.Model_Schema,
) -> (
	params: Model_Parameters,
	bd: fem.Basis_Descriptor,
	scheme: fem.Time_Scheme,
	ok: bool,
) {
    config_field, exists := schema.fields["temperature"]

    if !exists {
		log.error("Conduction model expect a field `temperature` to be defined.");
		return {}, {}, {}, false
	}

	scheme = config.assign_default(config_field.time_scheme, fem.Time_Scheme.BDF2)
	bd.order = config.assign_default(config_field.order, fem.Basis_Order.Linear)
	bd.family = config.assign_default(config_field.basis_family, fem.Basis_Family.Lagrange)

	if scheme == .Newmark{
	   log.error("Selected time scheme is not valid for temperature field.")
	   return {}, {}, {}, false
	}

	if bd.family != .Lagrange{
	   log.error("Selected basis family is not valid for temperature field.")
	   return {}, {}, {}, false
	}

    params.isothermal_bcs = make(map[fem.Boundary_ID]Isothermal_Int)
	params.materials = make(map[fem.Section_ID]Material_Int)
	params.variational_bcs = make(map[fem.Boundary_ID][]Variational_Int)
	params.sources = make(map[fem.Section_ID][]Source_Int)
	params.ics = make(map[fem.Section_ID]f64)

    // initial conditions
    for section_name, ic_config in config_field.initial_conditions {
		id, exists := mesh.section_names[section_name]
		if !exists {
			log.errorf("section %s is not declared on the mesh.", section_name);
			return {}, {}, {}, false
		}

		type := config.property_get(string, ic_config, "type") or_return
		switch type {
		case "constant":
			params.ics[id] = config.property_get(f64, ic_config, "temperature") or_return
		case:
			log.errorf("%s is not a recognized initial condition type", type);
			return {}, {}, {}, false
		}
	}

    // boundary conditions
	for bnd_name, bnd_configs in config_field.boundaries {
		id, exists := mesh.boundary_names[bnd_name]
		if !exists {
			log.errorf("boundary %s is not declared on the mesh.", bnd_name);
			return {}, {}, {}, false
		}

		variational := make([dynamic]Variational_Int)
		defer params.variational_bcs[id] = variational[:]


		// these BC's are not additive, would be wrong to define a isothermal, insulated, convective boundary.
		had_iso: bool
		had_adiabatic: bool

		for bnd_config in bnd_configs {
			bnd_type := config.property_get(string, bnd_config, "type") or_return
			switch bnd_type {
			case "isothermal":
				had_iso = true
				params.isothermal_bcs[id] = {
				    data = transmute(rawptr)config.property_get(f64, bnd_config, "temperature") or_return,
				    procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr) -> f64 {
				        return transmute(f64)data
				    }
				   }

			case "adiabatic":
				had_adiabatic = true
			case:
				log.errorf("%s is not a recognized boundary for temperature.", bnd_type)
			}
		}
		if (had_iso && had_adiabatic) || (had_iso && (len(variational) > 0)) {
			log.errorf("Boundary %s was declared isothermal no other boundaries may be applied.", bnd_name)
			return {}, {}, {}, false
		}

		if had_adiabatic && len(variational) > 0 {
			log.errorf("Boundary %s was declared adiabatic no other boundaries may be applied.", bnd_name)
			return {}, {}, {}, false
		}

	}

    for section_name, section_config in schema.sections {
		id, exists := mesh.section_names[section_name]
		if !exists {
			log.errorf("section %s is not declared on the mesh.", section_name);
			return {}, {}, {}, false
		}

		mat := section_config.material

		mat_type := config.property_get(string, mat, "type") or_return

		switch mat_type {
		case "constant":
			data := new(Constant_Material_Data)
			data.k = config.property_get(f64, mat, "conductivity") or_return
			data.rho = config.property_get(f64, mat, "density") or_return
			data.cp = config.property_get(f64, mat, "specific_heat") or_return
			params.materials[id] = {
				data = data,
				procedure = proc(
					mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Material
				) {
					info := cast(^Constant_Material_Data)data
					slice.fill(out.k, info.k); slice.fill(out.rho, info.rho); slice.fill(out.cp, info.cp)
				},
			}
        }
	}

	return params, bd, scheme, true

	Constant_Material_Data :: struct {
		k, rho, cp: f64,
	}

}

Output_Path_Format :: enum {
    PVD,
    VTU,
}

output_path :: proc(format: Output_Path_Format, dir, sim_name: string, step: int = 0, allocator: mem.Allocator) -> string {
    context.allocator = allocator

    filename: string
    switch format{
        case .PVD:
            filename = fmt.aprintf("%s.pvd", sim_name)
        case .VTU:
            filename = fmt.aprintf("%s_%d.vtu", sim_name, step)

    }

	path, _ := os.join_path({dir, filename}, allocator)
	return path
}