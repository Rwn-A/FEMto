package config

import "core:encoding/json"
import "core:fmt"
import "core:log"
import "core:mem"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:strings"

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
	cfg.sim_name = strings.clone(schema.sim_name)
	cfg.solver = schema.solver

	cfg.output = schema.output
	cfg.output.directory = strings.clone(schema.output.directory)

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

active_region_ids :: proc(cfg: General_Config, tbl: ^toml.Table, region_key: string, allocator := context.temp_allocator) -> ([]fem.Section_ID, bool) {
    regions, has_regions := table_get_opt(^toml.List, tbl, region_key)
    ids := make([dynamic]fem.Section_ID, allocator)
    for section_name, id in cfg.mesh.section_names {
        if !has_regions {
            append(&ids, id)
            continue
        }
        for region in regions {
            region_str, is_str := region.(string)
            if !is_str {
                log.error("region must be a string")
                return nil, false
            }
            if region_str not_in cfg.mesh.section_names {
                log.errorf("region %s not defined on mesh", region_str)
                return nil, false
            }
            if region_str == section_name do append(&ids, id)
        }
    }
    return ids[:], true
}