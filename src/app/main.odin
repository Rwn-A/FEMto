package main

import "core:encoding/json"
import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:os/os2"
import "core:path/filepath"

import fem "../fe_core"
import "../la"
import models "../models"
import fio "../serialization"


main :: proc() {
	context.logger = log.create_console_logger()

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {log.panic("Allocation failure")}
	defer virtual.arena_destroy(&arena)
	context.allocator = virtual.arena_allocator(&arena)

	if !run() {
		os2.exit(1)
	}

}

run :: proc() -> bool {
	fem.setup_quadrature_rules()
	fem.setup_subcell_rules()

	if len(os2.args) < 2 {
		log.error("Expected a path to a config file."); return false
	}

	config_path := os2.args[1]

	cs := read_config(config_path) or_return

	mesh, warn, mesh_err := fio.gmsh_parse(cs.mesh)

	if mesh_err != nil {
		log.errorf("Failed to read mesh because %v", mesh_err); return false
	}

	if warn != nil {
		log.warn("Mesh warning: %v", warn)
	}

	sc := models.Solver_Config{}

	sc.output = {
		output_dir  = assign_default(cs.output.directory, "./"),
		frequency   = assign_default(cs.output.frequency, 1),
		prefix      = cs.name,
		output_rule = .Split_Linear if cs.output.order == .Linear else .Split_Quadratic,
	}

	sc.tc = {
		is_transient = (cs.transient != {}),
		timestep     = cs.transient.timestep,
		current_time = cs.transient.start,
		end_time     = cs.transient.end,
	}

	sc.linsolve_rtol = assign_default(cs.linear_solver.rtol, 1e-7)
	sc.non_linear_rtol = assign_default(cs.non_linear_solver.rtol, 1e-7)
	sc.linsolve_max_iter = assign_default(cs.linear_solver.max_iterations, 1000)
	sc.non_linear_max_iter = assign_default(cs.non_linear_solver.max_iterations, 100)

	model: models.Model
	switch cs.model {
	case .Conduction:
		model = configure_conduction(cs, mesh) or_return
	}

	res := models.solve_model(sc, &model, mesh)

	if res.converged == false {
		log.error(res)
	} else {
		log.info("Simulation complete!")
	}
	return true
}

Field_Config :: struct {
	order:      fem.Order,
	boundaries: map[string]Property_Config,
}

Property_Value :: union {
	f64,
	string,
	int,
}

Property_Config :: map[string]Property_Value

Config_Schema :: struct {
	name:              string,
	mesh:              string,
	model:             enum {
		Conduction,
	},
	output:            struct {
		directory: string,
		frequency: int,
		order:     fem.Order,
	},
	transient:         struct {
		timestep, start, end: f64,
	},
	linear_solver:     struct {
		max_iterations: int,
		rtol:           f64,
	},
	non_linear_solver: struct {
		max_iterations: int,
		rtol:           f64,
	},
	fields:            map[string]Field_Config,
}

assign_default :: proc(a: $T, default: T) -> T {
	if a != {} {return a}
	return default
}

read_config :: proc(config_path: string) -> (Config_Schema, bool) {
	cs := Config_Schema{}

	json_data, err := os2.read_entire_file_from_path(config_path, context.allocator)

	if err != nil {
		log.errorf(os2.error_string(err))
	}

	config_dir := filepath.dir(config_path)
	os2.chdir(config_dir)

	json.unmarshal(json_data, &cs, .MJSON)

	return cs, true

}

configure_conduction :: proc(cs: Config_Schema, mesh: fem.Mesh) -> (m: models.Model, ok: bool) {
	model_params := new(models.Conduction_Params)

	config_field, exists := cs.fields["temperature"]

	if !exists {
		log.error("Conduction model expect a field \"temperature\" to be defined."); return {}, false
	}

	model_params.soln_order = config_field.order
	model_params.isothermal_bnds = make(map[fem.Boundary_ID]f64)

	for bnd_name, bnd_config in config_field.boundaries {
		id, exists := mesh.boundary_names[bnd_name]
		if !exists {
			log.errorf("boundary %s is not declared on the mesh.", bnd_name); return {}, false
		}
		bnd_type := property_get(string, bnd_config, "type") or_return
		switch bnd_type {
		case "isothermal":
			model_params.isothermal_bnds[id] = property_get(f64, bnd_config, "temperature") or_return
		case "adiabatic": // no op
		case:
			log.errorf("%s is not a recognized boundary for temperature.", bnd_type)
		}


	}


	model := models.model_conduction(model_params)

	return model, true
}

property_get :: proc($T: typeid, pc: Property_Config, key: string) -> (T, bool) {
	val, exists := pc[key]

	if !exists {
		log.errorf("Expected a property %s of type %v", key, type_info_of(T)); return {}, false
	}

	unwrapped, ok := val.(T)

	if !ok {
		log.errorf("Property %s was expected to be of type %v", type_info_of(T)); return {}, false
	}

	return unwrapped, true


}
