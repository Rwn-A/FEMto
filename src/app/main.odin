package main

import "core:encoding/json"
import "core:log"
import "core:math/linalg"
import "core:mem/virtual"
import "core:os/os2"
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
	boundaries: map[string][]Property_Config,
	initial_conditions: map[string]Property_Config,
}

Section_Config :: struct {
	material: Property_Config,
	sources: []Property_Config,
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
	sections: 		   map[string]Section_Config,
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
	model_params.materials = make(map[fem.Section_ID]models.Conduction_Material_Int)
	model_params.variational_bcs = make(map[fem.Boundary_ID][]models.Conduction_BC_Int)
	model_params.sources = make(map[fem.Section_ID][]models.Conduction_Source_Int)
	model_params.ics = make(map[fem.Section_ID]f64)

	for bnd_name, bnd_configs in config_field.boundaries {
		id, exists := mesh.boundary_names[bnd_name]
		if !exists {
			log.errorf("boundary %s is not declared on the mesh.", bnd_name); return {}, false
		}
		variational := make([dynamic]models.Conduction_BC_Int)
		defer model_params.variational_bcs[id] = variational[:]


		// these BC's are not additive, would be wrong to define a isothermal, insulated, convective boundary.
		had_iso: bool
		had_adiabatic: bool

		for bnd_config in bnd_configs {
			bnd_type := property_get(string, bnd_config, "type") or_return
			switch bnd_type {
			case "isothermal":
				had_iso = true
				model_params.isothermal_bnds[id] = property_get(f64, bnd_config, "temperature") or_return
			case "adiabatic":
				had_adiabatic = true
			case "flux":
				data := new(Fixed_Flux_Data)
				data.flux =  property_get(f64, bnd_config, "prescribed_flux") or_return
				append(&variational, models.Conduction_BC_Int{
					data = data,
					procedure = proc(ctx: fem.Element_Context, start, end: int, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_BC) {
						info := cast(^Fixed_Flux_Data)data
						slice.fill(out.q, info.flux)
					}
				}
				)
			case "convective":
				data := new(Convective_Data)
				data.h =  property_get(f64, bnd_config, "convective_coefficient") or_return
				data.ambient = property_get(f64, bnd_config, "ambient_temp") or_return
				append(&variational, models.Conduction_BC_Int{
					data = data,
					procedure = proc(ctx: fem.Element_Context, start, end: int, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_BC) {
						info := cast(^Convective_Data)data
						slice.fill(out.h, info.h)
						slice.fill(out.T_amb, info.ambient)
					}
				})
			case "radiative":
			    data := new(Radiative_Data)
			    data.emissivity = property_get(f64, bnd_config, "emissivity") or_return
			    data.ambient    = property_get(f64, bnd_config, "ambient_temp") or_return
			    data.stefan_boltzmann = property_get(f64, bnd_config, "stefan_boltzmann") or_return
			    append(&variational, models.Conduction_BC_Int{
			    	 data = data,
			        procedure = proc(ctx: fem.Element_Context, start, end: int, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_BC) {
			            info := cast(^Radiative_Data)data
			            sigma := info.stefan_boltzmann
			            for i in 0..<(end-start) {
			                T     := current_soln[start + i]
			                T_amb := info.ambient
			                out.q[i]   = info.emissivity * sigma * (T_amb*T_amb*T_amb*T_amb - T*T*T*T)
			                out.d_q[i] = -4 * info.emissivity * sigma * T*T*T
			            }
			        }
			    })
			case:
				log.errorf("%s is not a recognized boundary for temperature.", bnd_type)
			}
		}
		if (had_iso && had_adiabatic) || (had_iso && (len(variational) > 0)) {
			log.errorf("Boundary %s was declared isothermal no other boundaries may be applied.", bnd_name)
			return {}, false
		}

		if had_adiabatic && len(variational) > 0 {
			log.errorf("Boundary %s was declared adiabatic no other boundaries may be applied.", bnd_name)
		}

	}

	for section_name, ic_config in config_field.initial_conditions {
		id, exists := mesh.section_names[section_name]
		if !exists {
			log.errorf("section %s is not declared on the mesh.", section_name); return {}, false
		}

		type := property_get(string, ic_config, "type") or_return
		switch type{
			case "constant":
				model_params.ics[id] = property_get(f64, ic_config, "temperature") or_return
			case:
				log.errorf("%s is not a recognized initial condition type", type); return {}, false
		}
	}

	for section_name, section_config in cs.sections {
		id, exists := mesh.section_names[section_name]
		if !exists {
			log.errorf("section %s is not declared on the mesh.", section_name); return {}, false
		}

		mat := section_config.material

		mat_type := property_get(string, mat, "type") or_return

		switch mat_type{
			case "constant":
				data := new(Constant_Material_Data)
				data.k = property_get(f64, mat, "conductivity") or_return
				data.rho = property_get(f64, mat, "density") or_return
				data.cp = property_get(f64, mat, "specific_heat") or_return
				model_params.materials[id] = {
					data = data,
					procedure =  proc(ctx: fem.Element_Context, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_Material) {
						info := cast(^Constant_Material_Data)data
						slice.fill(out.k, info.k); slice.fill(out.rho, info.rho); slice.fill(out.cp, info.cp)
					}
				}
			case "linear":
				data := new(Linear_Material_Data)
				data.k0 = property_get(f64, mat, "reference_conductivity") or_return
				data.k1 = property_get(f64, mat, "conductivity_coefficient", true) or_else 0

				data.cp0 = property_get(f64, mat, "reference_specific_heat") or_return
				data.cp1 = property_get(f64, mat, "specific_heat_coefficient", true) or_else 0

				data.rho0 = property_get(f64, mat, "reference_density") or_return
				data.rho1 = property_get(f64, mat, "density_coefficient", true) or_else 0

				data.reference_t = property_get(f64, mat, "reference_temperature") or_return

				model_params.materials[id] = {
					data = data,
					procedure =  proc(ctx: fem.Element_Context, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_Material) {
						info := cast(^Linear_Material_Data)data
						slice.fill(out.d_k, info.k1)
						slice.fill(out.d_cp, info.cp1)
						slice.fill(out.d_rho, info.rho1)
						for temperature, i in current_soln {
							out.k[i] = info.k0 + info.k1 * (temperature - info.reference_t)
							out.cp[i] = info.cp0 + info.cp1 * (temperature - info.reference_t)
							out.rho[i] = info.rho0 + info.rho1 * (temperature - info.reference_t)
						}
					}
				}
			case:
				log.errorf("%s is not a recognized material type.", mat_type); return {}, false
		}
		model_params.sources[id] = make([]models.Conduction_Source_Int, len(section_config.sources))
		for source, i in section_config.sources {
			type := property_get(string, source, "type") or_return
			switch type {
			case "constant":
				data := new(Constant_Source_Data)
				data.Q = property_get(f64, source, "heat_source") or_return
				model_params.sources[id][i] = {
					data = data,
					procedure =  proc(ctx: fem.Element_Context, current_time: f64, current_soln: []f64, data: rawptr, out: models.Conduction_Source) {
						info := cast(^Constant_Source_Data)data
						for &entry in out.Q{entry += info.Q}
					}
				}
			case:
				log.errorf("%s is not a recognized source term", type); return {}, false
			}
		}

	}


	model := models.model_conduction(model_params)

	Constant_Material_Data :: struct {
		k, rho, cp: f64
	}

	// P = P0 + P1 * T
	Linear_Material_Data :: struct {
		k0, k1: f64,
		rho0, rho1: f64,
		cp0, cp1: f64,
		reference_t: f64,
	}

	Constant_Source_Data :: struct {
		Q: f64,
	}

	Fixed_Flux_Data :: struct {
		flux: f64,
	}

	Convective_Data :: struct {
		h: f64,
		ambient: f64,
	}

	Radiative_Data :: struct {
    	emissivity: f64,
    	ambient:    f64,
    	stefan_boltzmann: f64,
	}


	return model, true
}

property_get :: proc($T: typeid, pc: Property_Config, key: string, optional := false) -> (T, bool) {
	val, exists := pc[key]

	if !exists {
		if !optional {log.errorf("Expected a property %s of type %v", key, type_info_of(T))}
		return {}, false
	}

	unwrapped, ok := val.(T)

	if !ok {
		if !optional {log.errorf("Property %s was expected to be of type %v", type_info_of(T))}
		 return {}, false
	}

	return unwrapped, true
}
