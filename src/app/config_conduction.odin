// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package main

import fem "../fe_core"
import "../la"
import models "../models"
import fio "../serialization"
import "core:log"
import "core:slice"

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
				data.flux = property_get(f64, bnd_config, "prescribed_flux") or_return
				append(&variational, models.Conduction_BC_Int {
					data = data,
					procedure = proc(
						ctx: fem.Element_Context,
						start, end: int,
						current_time: f64,
						current_soln: []f64,
						data: rawptr,
						out: models.Conduction_BC,
					) {
						info := cast(^Fixed_Flux_Data)data
						slice.fill(out.q, info.flux)
					},
				})
			case "convective":
				data := new(Convective_Data)
				data.h = property_get(f64, bnd_config, "convective_coefficient") or_return
				data.ambient = property_get(f64, bnd_config, "ambient_temp") or_return
				append(&variational, models.Conduction_BC_Int {
					data = data,
					procedure = proc(
						ctx: fem.Element_Context,
						start, end: int,
						current_time: f64,
						current_soln: []f64,
						data: rawptr,
						out: models.Conduction_BC,
					) {
						info := cast(^Convective_Data)data
						slice.fill(out.h, info.h)
						slice.fill(out.T_amb, info.ambient)
					},
				})
			case "radiative":
				data := new(Radiative_Data)
				data.emissivity = property_get(f64, bnd_config, "emissivity") or_return
				data.ambient = property_get(f64, bnd_config, "ambient_temp") or_return
				data.stefan_boltzmann = property_get(f64, bnd_config, "stefan_boltzmann") or_return
				append(&variational, models.Conduction_BC_Int {
					data = data,
					procedure = proc(
						ctx: fem.Element_Context,
						start, end: int,
						current_time: f64,
						current_soln: []f64,
						data: rawptr,
						out: models.Conduction_BC,
					) {
						info := cast(^Radiative_Data)data
						sigma := info.stefan_boltzmann
						for i in 0 ..< (end - start) {
							T := current_soln[start + i]
							T_amb := info.ambient
							out.q[i] = info.emissivity * sigma * (T_amb * T_amb * T_amb * T_amb - T * T * T * T)
							out.d_q[i] = -4 * info.emissivity * sigma * T * T * T
						}
					},
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
		switch type {
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

		switch mat_type {
		case "constant":
			data := new(Constant_Material_Data)
			data.k = property_get(f64, mat, "conductivity") or_return
			data.rho = property_get(f64, mat, "density") or_return
			data.cp = property_get(f64, mat, "specific_heat") or_return
			model_params.materials[id] = {
				data = data,
				procedure = proc(
					ctx: fem.Element_Context,
					current_time: f64,
					current_soln: []f64,
					data: rawptr,
					out: models.Conduction_Material,
				) {
					info := cast(^Constant_Material_Data)data
					slice.fill(out.k, info.k); slice.fill(out.rho, info.rho); slice.fill(out.cp, info.cp)
				},
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
				procedure = proc(
					ctx: fem.Element_Context,
					current_time: f64,
					current_soln: []f64,
					data: rawptr,
					out: models.Conduction_Material,
				) {
					info := cast(^Linear_Material_Data)data
					slice.fill(out.d_k, info.k1)
					slice.fill(out.d_cp, info.cp1)
					slice.fill(out.d_rho, info.rho1)
					for temperature, i in current_soln {
						out.k[i] = info.k0 + info.k1 * (temperature - info.reference_t)
						out.cp[i] = info.cp0 + info.cp1 * (temperature - info.reference_t)
						out.rho[i] = info.rho0 + info.rho1 * (temperature - info.reference_t)
					}
				},
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
					procedure = proc(
						ctx: fem.Element_Context,
						current_time: f64,
						current_soln: []f64,
						data: rawptr,
						out: models.Conduction_Source,
					) {
						info := cast(^Constant_Source_Data)data
						for &entry in out.Q {entry += info.Q}
					},
				}
			case:
				log.errorf("%s is not a recognized source term", type); return {}, false
			}
		}

	}


	model := models.model_conduction(model_params)

	Constant_Material_Data :: struct {
		k, rho, cp: f64,
	}

	// P = P0 + P1 * T
	Linear_Material_Data :: struct {
		k0, k1:      f64,
		rho0, rho1:  f64,
		cp0, cp1:    f64,
		reference_t: f64,
	}

	Constant_Source_Data :: struct {
		Q: f64,
	}

	Fixed_Flux_Data :: struct {
		flux: f64,
	}

	Convective_Data :: struct {
		h:       f64,
		ambient: f64,
	}

	Radiative_Data :: struct {
		emissivity:       f64,
		ambient:          f64,
		stefan_boltzmann: f64,
	}


	return model, true
}
