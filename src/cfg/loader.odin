/*
 Defines the top level schema for the configuration file and helpers for models to parse their specific configuration.
*/
package cfg

import "core:os"
import "core:log"
import "core:path/filepath"
import "core:dynlib"
import "core:slice"

import "../fem"
import fio"../fem/serialization"

MAX_PLUGINS :: 16

Mesh_Schema :: struct {
    path: string `validate:"required"`
}

Output_Format :: enum {
    CSV,
    VTU,
}

Output_Schema :: struct {
    directory: string,
    frequency: int,
    order: fem.Basis_Order,
    formats: []Output_Format,
    post_process: ^Table, // model specific hook
}

Solver_Schema :: struct {
    max_linear_iters:    int `validate:"gte=1"`,
	max_nonlinear_iters: int `validate:"gte=1"`,
	tolerance:           f64 `validate:"gt=0"`,
	override:           ^Table, // model specific hook
}

Transient_Schema :: struct {
	timestep: f64 `validate:"required,gt=0"`,
	end:      f64 `validate:"required"`,
	start:    f64,
	scheme:  ^Table, // model specific hook
}

Root_Schema :: struct {
	sim_name:          string,
	model:             string `validate:"required"`,
	mesh:              Mesh_Schema `validate:"required"`,
	solver:            Solver_Schema,
	output:            Output_Schema,
	time_control:      Transient_Schema, // an emtpy struct indicates steady-state
	plugins:           []string,
	boundaries:          ^Table, // model specific hook
	sources:            ^Table, // model specific hook
	materials:          ^Table, // model specific hook
	initial_conditions: ^Table, // model specific hook
	discretization:    ^Table, // model specific hook
}

@(rodata, private)
DEFAULT_CONFIG := Root_Schema {
	sim_name = "My_Sim",
	solver = {max_linear_iters = 1000, max_nonlinear_iters = 25, tolerance = 1e-7},
	output = {directory = "results", frequency = 1, order = .Linear, formats = {.VTU}},
}

Plugin :: #type proc(tbl: ^Table, ctx: Plugin_Context) -> bool
Plugin_Registry :: map[string]Plugin

// this is the type of procedure the load plugin looks for in the dll
Register_Proc :: #type proc "c" (reg: ^Plugin_Registry)

// not really used by the plugin itself but needed for the model to verify
// what the plugin is allowed to do.
Plugin_Context :: struct {
	current_section:                        fem.Section_ID,
	current_boundary:                       fem.Boundary_ID,
	current_section_name, current_bnd_name: string,
	allowed_kind:                           Plugin_Kind, // plugin was called in a `source` block for example, so should only make a source
	problem_data:                       any,
}

Plugin_Kind :: enum {
	Source,
	BC,
	Material,
	IC,
}

load_root :: proc(path: string, allocator := context.allocator) -> (schema: Root_Schema, ok: bool) {
    context.allocator = allocator

	data, err := os.read_entire_file_from_path(path, context.allocator)

	if err != nil {
		log.fatalf(os.error_string(err)); return {}, false
	}

	os.chdir(filepath.dir(path))

	root, toml_err := parse(string(data), path)

	if toml_err.type != .None {
		log.fatal(format_error(&toml_err)); return {}, false
	}

	schema = DEFAULT_CONFIG
	unmarshal(root, &schema) or_return

	return schema, true
}

load_mesh :: proc(mesh_schema: Mesh_Schema, allocator := context.allocator) -> (fem.Mesh, bool) {
    log.infof("Loading mesh %s", mesh_schema.path)

	mesh, warn, err := fio.gmsh_parse(mesh_schema.path, allocator)
	if err != nil {log.error(err); return {}, false}
	if warn != nil {log.warn(warn)}

	return mesh, true
}

load_plugins :: proc(plugins: []string, allocator := context.allocator) -> (libs: [dynamic; MAX_PLUGINS]dynlib.Library, reg: Plugin_Registry, ok: bool) {
    context.allocator = allocator

	reg = make(Plugin_Registry)

	for path in plugins {
		lib, ok := dynlib.load_library(path)

		if !ok {
			log.errorf("config: Failed to load plugin at %s", path)
			return {}, {}, false
		}

		if append(&libs, lib) == 0 {
		  log.errorf("Maximum plugin limit of %d reached", MAX_PLUGINS)
		  return {}, {}, false
		}

		register_sym, found := dynlib.symbol_address(lib, "register")

		if !found {
			log.errorf("config: Plugin at %s missing register symbol", path)
			return {}, {}, false
		}

		register := cast(Register_Proc)register_sym

		register(&reg)
	}
    return libs, reg, true
}

unload_plugins :: proc(libs: [dynamic; MAX_PLUGINS]dynlib.Library) {
    for lib in libs{ dynlib.unload_library(lib) }
}

// helpers for models

@(private = "file")
Base_Region_Schema :: struct {
	kind:    string `validate:"required"`,
	regions: []string,
}

@(private = "file")
Base_Boundary_Schema :: struct {
	kind:       string `validate:"required"`,
	boundaries: []string `validate:"required"`,
}

load_bc_plugins :: proc(bcs: []^Table, mesh: fem.Mesh, registry: Plugin_Registry, ctx: ^Plugin_Context) -> bool {
	if bcs == nil {return true}
	ctx.allowed_kind = .BC

	for bc in bcs {
		base: Base_Boundary_Schema
		unmarshal(bc, &base) or_return
		plugin, exists := registry[base.kind]
		if !exists {
			log.errorf("config: BC kind %s could not be found.", base.kind); return false
		}
		if len(base.boundaries) == 0 {
			log.error("config: Must provide at least one boundary name for a boundary condition.")
			return false
		}
		for boundary in base.boundaries {
			id, exists := mesh.boundary_names[boundary]
			if !exists {
				log.errorf("config: Boundary %s does not exist on provided mesh.", boundary)
				names, _ := slice.map_keys(mesh.boundary_names)
				log.infof("config: Available boundaries: %s", names)
				return false
			}
			ctx.current_boundary = id
			ctx.current_bnd_name = boundary
			plugin(bc, ctx^) or_return
		}
	}
	return true
}

load_section_plugins :: proc(
	entries: []^Table,
	mesh: fem.Mesh,
	registry: Plugin_Registry,
	ctx: ^Plugin_Context,
	kind: Plugin_Kind,
) -> bool {
	if entries == nil { return true }
	ctx.allowed_kind = kind

	for entry in entries {
		base: Base_Region_Schema
		unmarshal(entry, &base) or_return
		plugin, exists := registry[base.kind]
		if !exists {
			log.errorf("config: kind %s could not be found.", base.kind)
			return false
		}
		if len(base.regions) == 0 {
			for name, id in mesh.section_names {
				ctx.current_section = id
				ctx.current_section_name = name
				plugin(entry, ctx^) or_return
			}
		} else {
			for region in base.regions {
				id, exists := mesh.section_names[region]
				if !exists {
					log.errorf("config: Region %s does not exist on provided mesh.", region)
					names, _ := slice.map_keys(mesh.section_names)
					log.infof("config: Available regions: %s", names)
					return false
				}
				ctx.current_section = id
				ctx.current_section_name = region
				plugin(entry, ctx^) or_return
			}
		}
	}
	return true
}

// usually in time_control.scheme or a sub table of it
load_time_scheme :: proc(allowed: bit_set[fem.Time_Scheme], default: fem.Time_Scheme, scheme_tbl: ^Table) -> (scheme: fem.Time_Scheme, ok: bool) {
    res: struct {scheme: fem.Time_Scheme} = {default}
    if scheme_tbl == nil {return default, true}
    unmarshal(scheme_tbl, &res)
    valid_enum_option(allowed, res.scheme) or_return
    return res.scheme, true
}

// solver.override or a subtable
load_solver_kind :: proc(allowed: bit_set[fem.Solver_Kind], default: fem.Solver_Kind, table : ^Table) -> (kind: fem.Solver_Kind, ok: bool) {
    res: struct {solver: fem.Solver_Kind} = {default}
    if table == nil {return default, true}
    unmarshal(table, &res)
    valid_enum_option(allowed, res.solver) or_return
    return res.solver, true
}

// discretization or a subtable
load_discretization :: proc(allowed: bit_set[fem.Basis_Family], default: fem.Basis_Family, table: ^Table) -> (bd: fem.Basis_Descriptor, ok: bool) {
    res: struct {family: fem.Basis_Family, order: fem.Basis_Order} = {default, .Linear}
    if table == nil {return {default, .Linear}, true}
    unmarshal(table, &res)
    valid_enum_option(allowed, res.family) or_return
    return {res.family, res.order}, true
}


// for bcs, mat props, sources
get_field_table :: proc(base: ^Table, field_name: string) -> ([]^Table, bool, bool) {
	if base == nil {return {}, false, true}
    if field_name not_in base {return {}, false, true}

    x: []^Table
    field, ok := table_get(^List, base, field_name)
    if !ok {return {}, false, false}
    if !unmarshal_value(x, field, field_name, {}) {return {}, false, false}

    return x, true, true
}

valid_enum_option :: proc(allowed: bit_set[$T], got: T) -> bool {
	if got not_in allowed {
		log.errorf("config: %s is not a valid option, expected one of %v", got, allowed)
		return false
	}
	return true
}

