package config

import "core:encoding/json"
import "core:fmt"
import "core:log"
import "core:mem"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"

import "../../../vendor/toml"

import "../../fem"
import fio "../../fem/serialization"

Solver_Schema :: struct {
	max_linear_iters:    int,
	max_nonlinear_iters: int,
	tolerance:           f64,
}

Transient_Schema :: struct {
	timestep, start, end: f64 `toml:"required"`,
}

Output_Schema :: struct {
	directory: string,
	frequency: int,
	order:     fem.Basis_Order,
}

Mesh_Schema :: struct {
	path: string `toml:"required"`,
}


Schema :: struct {
	sim_name:     string,
	mesh:         Mesh_Schema `toml:"required"`,
	model:        string `toml:"required"`,
	solver:       Solver_Schema,
	output:       Output_Schema,
	time_control: Transient_Schema,

	//lazy, schema depends on chosen model
	materials:    ^toml.Table `toml:"required"`,
	boundaries:   ^toml.Table `toml:"required"`,
	sources:      ^toml.List,
	fields:       ^toml.Table,
}

DEFAULT_SCHEMA :: Schema {
	sim_name = "My_Sim",
	solver = {max_linear_iters = 1000, max_nonlinear_iters = 25, tolerance = 1e-8},
	output = {directory = "results", frequency = 1, order = .Linear},
}


// not model specific, sensible defaults have been filled.
General_Config :: struct {
	mesh:      fem.Mesh,
	viz_mesh:  fio.VTK_Raw_Mesh,
	sim_name:  string,
	output:    Output_Schema,
	transient: Maybe(Transient_Schema),
	solver:    Solver_Schema,
}

load_general_config :: proc(schema: Schema, arena: ^virtual.Arena) -> (cfg: General_Config, ok: bool) {
	context.allocator = virtual.arena_allocator(arena)

	log.infof("Loading mesh %s", schema.mesh.path)
	mesh, warn, err := fio.gmsh_parse(schema.mesh.path)
	if err != nil {log.error(err); return {}, false}
	if warn != nil {log.warn(warn)}

	cfg.mesh = mesh
	cfg.sim_name = schema.sim_name
	cfg.solver = schema.solver

	cfg.output = schema.output

	if cfg.output.frequency <= 0 {
		log.error("Output frequency must be positive and non-zero.")
		return {}, false
	}

	if schema.time_control == {} {
		cfg.transient = nil
	} else {
		if schema.time_control.end <= schema.time_control.start {
			log.error("config: Start time must be strcitly less than end time.")
			return {}, false
		}

		if schema.time_control.timestep <= 0 {
			log.error("config: timestep must be positive and non-zero.")
			return {}, false
		}
		cfg.transient = schema.time_control
	}

	cfg.viz_mesh = fio.vtk_create_visualization_mesh(mesh, cfg.output.order)

	if !os.exists(cfg.output.directory) {
		err := os.make_directory(cfg.output.directory)
		if err != nil {
			log.error(os.error_string(err))
			return {}, false
		}
	}


	return cfg, true
}

Output_Path_Format :: enum {
	PVD,
	VTU,
}

format_output_path :: proc(
	format: Output_Path_Format,
	dir, sim_name: string,
	step: int = 0,
	allocator: mem.Allocator,
) -> string {
	context.allocator = allocator

	filename: string
	switch format {
	case .PVD:
		filename = fmt.aprintf("%s.pvd", sim_name)
	case .VTU:
		filename = fmt.aprintf("%s_%d.vtu", sim_name, step)

	}

	path, _ := os.join_path({dir, filename}, allocator)
	return path
}

valid_option :: proc(allowed: bit_set[$T], got: T) -> bool {
	if got not_in allowed {
		log.errorf("config: %s is not a valid option, expected one of %v", got, allowed)
		return false
	}
	return true
}


// Property_Value :: union {
// 	f64,
// 	string,
// 	int,
// 	fem.Vec3,
// }

// Property_Config :: map[string]Property_Value

// Field_Config :: struct {
// 	order:              fem.Basis_Order,
// 	basis_family:       fem.Basis_Family,
// 	time_scheme:        fem.Time_Scheme,
// 	boundaries:         map[string][]Property_Config,
// 	initial_conditions: map[string]Property_Config,
// }

// Section_Config :: struct {
// 	material: Property_Config,
// 	sources:  []Property_Config,
// }


// Config_Schema :: struct {
// 	sim_name:      string,
// 	mesh_path:     string,
// 	model:         string,
// 	time_control:  Transient_Schema,
// 	solver:        Solver_Schema,
// 	output:        Output_Schema,
// 	linear_solver: fem.Solver_Kind,
// 	fields:        map[string]Field_Config,
// 	sections:      map[string]Section_Config,
// }

// Model_Schema :: struct {
// 	linear_solver: fem.Solver_Kind,
// 	fields:        map[string]Field_Config,
// 	sections:      map[string]Section_Config,
// }


// load_general_config :: proc(schema: Config_Schema, arena: ^virtual.Arena) -> (cfg: General_Config, ok: bool) {
// 	context.allocator = virtual.arena_allocator(arena)

// 	log.infof("Loading mesh %s", schema.mesh_path)
// 	mesh, warn, err := fio.gmsh_parse(schema.mesh_path)
// 	if err != nil {log.error(err); return {}, false}
// 	if warn != nil {log.warn(warn)}

// 	cfg.mesh = mesh

// 	cfg.sim_name = assign_default(schema.sim_name, "my_simulation")

// 	cfg.output.directory = assign_default(schema.output.directory, "output")
// 	cfg.output.frequency = assign_default(schema.output.frequency, 1)
// 	cfg.output.order = assign_default(schema.output.order, fem.Basis_Order.Linear)

// 	cfg.viz_mesh = fio.vtk_create_visualization_mesh(mesh, cfg.output.order)

// 	if !os.exists(cfg.output.directory) {
// 		os.make_directory(cfg.output.directory)
// 	}

// 	if schema.time_control == {} {
// 		cfg.transient = nil
// 	} else {
// 		if schema.time_control.timestep <= 0 {
// 			log.errorf("Timestep must be a positive non-zero number.")
// 			return {}, false
// 		}

// 		if schema.time_control.end <= schema.time_control.start {
// 			log.errorf("Start time must occur before end time.")
// 			return {}, false
// 		}

// 		cfg.transient = schema.time_control
// 	}

// 	cfg.solver.max_linear_iters = assign_default(schema.solver.max_linear_iters, 500)
// 	cfg.solver.max_nonlinear_iters = assign_default(schema.solver.max_nonlinear_iters, 10)
// 	cfg.solver.tolerance = assign_default(schema.solver.tolerance, 1e-6)

// 	return cfg, true
// }


// read_config :: proc(config_path: string, arena: ^virtual.Arena) -> (Config_Schema, bool) {
// 	context.allocator = virtual.arena_allocator(arena)
// 	cs := Config_Schema{}

// 	json_data, err := os.read_entire_file_from_path(config_path, context.allocator)

// 	if err != nil {
// 		log.error(os.error_string(err))
// 		return {}, false
// 	}

// 	config_dir := filepath.dir(config_path)
// 	os.chdir(config_dir)

// 	if err := json.unmarshal(json_data, &cs, .MJSON); err != nil {
// 		log.error(err)
// 		return {}, false
// 	}


// 	return cs, true

// }

// assign_default :: proc(a: $T, default: T) -> T {
// 	if a != {} {return a}
// 	return default
// }


// property_get :: proc($T: typeid, pc: Property_Config, key: string, optional := false) -> (T, bool) {
// 	val, exists := pc[key]

// 	if !exists {
// 		if !optional {log.errorf("Expected a property %s of type %v", key, type_info_of(T))}
// 		return {}, false
// 	}

// 	unwrapped, ok := val.(T)

// 	if !ok {
// 		if !optional {log.errorf("Property %s was expected to be of type %v", key, type_info_of(T))}
// 		return {}, false
// 	}

// 	return unwrapped, true
// }
