/*
 Defines the schema for the simulations config, and provides utilities for models to load their config.
*/
package cfg

import "core:log"
import "core:mem/virtual"
import "core:os"
import "core:path/filepath"
import "core:slice"
import "core:dynlib"

import "../../../vendor/toml"
import "../../fem"
import fio "../../fem/serialization"

Output_Format :: enum {
	CSV,
	VTU,
}

Solver_Schema :: struct {
	max_linear_iters:    int `validate:"gte=1"`,
	max_nonlinear_iters: int `validate:"gte=1"`,
	tolerance:           f64 `validate:"gt=0"`,
}

Transient_Schema :: struct {
	timestep:   f64 `validate:"required,gt=0"`,
	start, end: f64 `validate:"required"`,
}

Output_Schema :: struct {
	directory: string,
	frequency: int `validate:"gte=1"`,
	order:     fem.Basis_Order,
	formats:   []Output_Format,
}

Mesh_Schema :: struct {
	path: string `validate:"required"`,
}

Schema :: struct {
	sim_name:     string,
	model:        string `validate:"required"`,
	plugins: 	  []string,
	mesh:         Mesh_Schema `validate:"required"`,
	solver:       Solver_Schema,
	output:       Output_Schema,
	time_control: Transient_Schema,
	//lazy, depends on chosen model how these are parsed.
	fields:       ^Table,
	material:    []^Table,
	source:      []^Table,
	boundary_condition:          []^Table,
}

@(rodata)
DEFAULT_CONFIG := Schema {
	sim_name = "My_Sim",
	solver = {max_linear_iters = 1000, max_nonlinear_iters = 25, tolerance = 1e-8},
	output = {directory = "results", frequency = 1, order = .Linear, formats = {.VTU}},
}

General_Params :: struct {
	mesh:      fem.Mesh,
	transient: Maybe(Transient_Schema),
	solver:    Solver_Schema,
}

load_schema :: proc(path: string, allocator := context.allocator) -> (config: Schema, ok: bool) {
	context.allocator = allocator

	data, err := os.read_entire_file_from_path(path, context.allocator)

	if err != nil {
		log.fatalf(os.error_string(err))
		return {}, false
	}

	os.chdir(filepath.dir(path))

	root, toml_err := toml.parse(string(data), path)

	if toml_err.type != .None {
		log.fatal(toml.format_error(&toml_err))
		return {}, false
	}

	config = DEFAULT_CONFIG

	unmarshal(root, &config) or_return

	return config, true
}

load_general_params :: proc(schema: Schema, allocator := context.allocator) -> (params: General_Params, output: Outputter, ok: bool) {
	context.allocator = allocator

	params.solver = schema.solver

	log.infof("Loading mesh %s", schema.mesh.path)
	mesh, warn, err := fio.gmsh_parse(schema.mesh.path)
	if err != nil {log.error(err); return {}, {}, false}
	if warn != nil {log.warn(warn)}

	params.mesh = mesh

	if schema.time_control == {} {
		params.transient = nil
	} else {
		if schema.time_control.end <= schema.time_control.start {
			log.error("config: Start time must be strictly less than end time.")
			return {}, {}, false
		}
		params.transient = schema.time_control
	}



	output = outputter_create(schema, mesh) or_return


	return params, output, true
}

valid_enum_option :: proc(allowed: bit_set[$T], got: T) -> bool {
	if got not_in allowed {
		log.errorf("config: %s is not a valid option, expected one of %v", got, allowed)
		return false
	}
	return true
}

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


Plugin_Context :: struct($Material, $Source, $BC, $IC: typeid) {
    materials: map[string]proc(^Table) -> (Material, bool),
    sources:   map[string]proc(^Table) -> (Source, bool),
    bcs:       map[string]proc(^Table) -> (BC, bool),
    ics:       map[string]proc(^Table) -> (IC, bool),
}

init_plugin_context :: proc(
	ctx: ^Plugin_Context($Material, $Source, $BC, $IC),
	 allocator := context.allocator
) -> Plugin_Context(Material, Source, BC, IC) {
	context.allocator = allocator
	return {
		materials = make(map[string]proc(^Table) -> (Material, bool)),
		sources = make(map[string]proc(^Table) -> (Source, bool)),
		bcs =     make(map[string]proc(^Table) -> (BC, bool)),
		ics = make(map[string]proc(^Table) -> (IC, bool)),

	}
}

load_plugin :: proc(path: string, ctx: ^Plugin_Context($Material, $Source, $BC, $IC)) -> (success:bool) {
    lib, ok := dynlib.load_library(path)

    defer if !dynlib.unload_library(lib){
		log.error(dynlib.last_error())
		success = false
    }

    if !ok {
        log.errorf("config: Failed to load plugin at %s", path)
        return false
    }

    register_sym, found := dynlib.symbol_address(lib, "register")

    if !found {
        log.errorf("config: Plugin at %s missing register symbol", path)
        return false
    }
    register := cast(proc "c" (^Plugin_Context(Material, Source, BC, IC)))register_sym
    register(ctx)
    return true
}


load_material :: proc(
	$T: typeid,
	schema: Schema,
	mesh: fem.Mesh,
	loaded_materials: map[string]proc(^Table) -> (T, bool),
	allocator := context.allocator,
) -> (
	res: map[fem.Section_ID]T,
	ok: bool,
) {
	res = make(map[fem.Section_ID]T)

	if len(schema.material) == 0 {
		log.error("config: Must define atleast one material.")
		return {}, false
	}

	for material_schema in schema.material {
		base: Base_Region_Schema
		unmarshal(material_schema, &base) or_return

		material_proc, exists := loaded_materials[base.kind]
		if !exists {
			log.errorf("config: Material kind %s could not be found.", base.kind)
			return res, false
		}
		material := material_proc(material_schema) or_return

		ids := resolve_region_ids(base.regions, mesh) or_return
		for id in ids {
			if id in res {
				log.errorf("config: Cannot define two materials over the same region.")
				return res, false
			}
			res[id] = material
		}
	}
	return res, true
}

load_sources :: proc(
	$T: typeid,
	schema: Schema,
	mesh: fem.Mesh,
	loaded_sources: map[string]proc(^Table) -> (T, bool),
	allocator := context.allocator,
) -> (
	res: map[fem.Section_ID][]T,
	ok: bool,
) {
	res = make(map[fem.Section_ID][]T)

	res_builder := make(map[fem.Section_ID][dynamic]T)
	for _, id in mesh.section_names {res_builder[id] = make([dynamic]T)}

	for source_schema in schema.source {
		base: Base_Region_Schema
		unmarshal(source_schema, &base) or_return

		source_proc, exists := loaded_sources[base.kind]
		if !exists {
			log.errorf("config: Source kind %s could not be found.", base.kind)
			return res, false
		}
		source := source_proc(source_schema) or_return

		ids := resolve_region_ids(base.regions, mesh) or_return
		for id in ids {append(&res_builder[id], source)}
	}
	for region, &sources in res_builder {res[region] = sources[:]}
	return res, true
}

load_bcs :: proc($BC, $C, $V: typeid, schema: Schema, mesh: fem.Mesh, loaded_bcs: map[string]proc(^Table) -> (BC, bool), allocator := context.allocator) -> (variational: map[fem.Boundary_ID][]V, constraint: map[fem.Boundary_ID]C, ok: bool) {
	variational = make(map[fem.Boundary_ID][]V)
	constraint = make(map[fem.Boundary_ID]C)
	var_builder := make(map[fem.Boundary_ID][dynamic]V)

	for _, id in mesh.boundary_names {var_builder[id] = make([dynamic]V)}

	if len(schema.boundary_condition) == 0 {
		log.error("config: Must provide atleast one boundary condition")
		return {}, {}, false
	}

	for bc_schema in schema.boundary_condition {
		base: Base_Boundary_Schema
		unmarshal(bc_schema, &base) or_return

		bc_proc, exists := loaded_bcs[base.kind]
		if !exists {
			log.errorf("config: BC kind %s could not be found.", base.kind)
			return variational, constraint, false
		}
		bc := bc_proc(bc_schema) or_return

		ids := resolve_boundary_ids(base.boundaries, mesh) or_return
		switch v in bc {
		case V:
			for id in ids {
				if id in constraint {
					log.errorf("config: Boundary already has a constraint BC, cannot apply variational BC.")
					return variational, constraint, false
				}
				append(&var_builder[id], v)
			}
		case C:
			for id in ids {
				if id in constraint {
					log.errorf("config: Cannot define two constraint BCs over the same boundary.")
					return variational, constraint, false
				}
				if len(var_builder[id]) > 0 {
					log.errorf("config: Boundary already has a variational BC, cannot apply constraint BC.")
					return variational, constraint, false
				}
				constraint[id] = v
			}
		}
	}

	for id, &bcs in var_builder {variational[id] = bcs[:]}
	return variational, constraint, true
}


load_ics :: proc(
	$T: typeid,
	ics: []^Table,
	mesh: fem.Mesh,
	loaded_ics: map[string]proc(^Table) -> (T, bool),
	allocator := context.allocator,
) -> (
	res: map[fem.Section_ID]T,
	ok: bool,
) {
	res = make(map[fem.Section_ID]T)

	for ic_schema in ics {
		base: Base_Region_Schema
		unmarshal(ic_schema, &base) or_return

		ic_proc, exists := loaded_ics[base.kind]
		if !exists {
			log.errorf("config: IC kind %s could not be found.", base.kind)
			return res, false
		}
		ic := ic_proc(ic_schema) or_return

		ids := resolve_region_ids(base.regions, mesh) or_return
		for id in ids {
			if id in res {
				log.errorf("config: Cannot define two initial conditions over the same region.")
				return res, false
			}
			res[id] = ic
		}
	}
	return res, true
}

@(private = "file")
resolve_region_ids :: proc(regions: []string, mesh: fem.Mesh) -> (ids: []fem.Section_ID, ok: bool) {
	result: [dynamic; 32]fem.Section_ID
	if len(regions) == 0 {
		for _, id in mesh.section_names {append(&result, id)}
	} else {
		for region in regions {
			id, exists := mesh.section_names[region]
			if !exists {
				log.errorf("config: Region %s does not exist on provided mesh.", region)
				return nil, false
			}
			if append(&result, id) == 0 {
				log.errorf("config: Maximum limit of boundary tags reached.")
				return {}, false
			}
		}
	}
	return result[:], true
}

@(private = "file")
resolve_boundary_ids :: proc(boundaries: []string, mesh: fem.Mesh) -> (ids: []fem.Boundary_ID, ok: bool) {
	result: [dynamic; 32]fem.Boundary_ID
	if len(boundaries) == 0 {
		log.error("Must provide atleast one boundary name for a boundary condition.")
		return {}, false
	} else {
		for boundary in boundaries {
			id, exists := mesh.boundary_names[boundary]
			if !exists {
				log.errorf("config: Boundary %s does not exist on provided mesh.", boundary)
				return nil, false
			}
			if append(&result, id) == 0 {
				log.errorf("config: Maximum limit of boundary tags reached.")
				return {}, false
			}
		}
	}
	return result[:], true
}
