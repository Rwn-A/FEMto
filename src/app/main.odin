// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package main

import "core:encoding/json"
import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:slice"

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
		os.exit(1)
	}

}

run :: proc() -> bool {
	fem.setup_quadrature_rules()
	fem.setup_subcell_rules()

	if len(os.args) < 2 {
		log.error("Expected a path to a config file."); return false
	}

	config_path := os.args[1]

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
	case .Elasticity:
		model = configure_elasticity(cs, mesh) or_return
	case:
		log.error("Unknown model"); return false
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
	order:              fem.Order,
	boundaries:         map[string][]Property_Config,
	initial_conditions: map[string]Property_Config,
}

Section_Config :: struct {
	material: Property_Config,
	sources:  []Property_Config,
}

Property_Value :: union {
	f64,
	string,
	int,
	fem.Vec3,
}

Property_Config :: map[string]Property_Value

Config_Schema :: struct {
	name:              string,
	mesh:              string,
	model:             enum {
		Conduction,
		Elasticity,
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
	sections:          map[string]Section_Config,
}

assign_default :: proc(a: $T, default: T) -> T {
	if a != {} {return a}
	return default
}

read_config :: proc(config_path: string) -> (Config_Schema, bool) {
	cs := Config_Schema{}

	json_data, err := os.read_entire_file_from_path(config_path, context.allocator)

	if err != nil {
		log.errorf(os.error_string(err))
		return {}, false
	}

	config_dir := filepath.dir(config_path)
	os.chdir(config_dir)

	if err := json.unmarshal(json_data, &cs, .MJSON); err != nil {
		log.info(err)
		return {}, false
	}


	return cs, true

}

property_get :: proc($T: typeid, pc: Property_Config, key: string, optional := false) -> (T, bool) {
	val, exists := pc[key]

	if !exists {
		if !optional {log.errorf("Expected a property %s of type %v", key, type_info_of(T))}
		return {}, false
	}

	unwrapped, ok := val.(T)

	if !ok {
		if !optional {log.errorf("Property %s was expected to be of type %v", key, type_info_of(T))}
		return {}, false
	}

	return unwrapped, true
}
