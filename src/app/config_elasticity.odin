// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package main

import fem "../fe_core"
import "../la"
import models "../models"
import fio "../serialization"
import "core:log"

import "core:slice"

configure_elasticity :: proc(cs: Config_Schema, mesh: fem.Mesh) -> (m: models.Model, ok: bool) {
	assert(
		mesh.dim == .D3,
		"Elasticity is only supported on 3D meshes, lower dimensional simulations require a different model.",
	)

	model_params := new(models.Elasticity_Params)

	config_field, exists := cs.fields["displacement"]

	if !exists {
		log.error("Elasticity model expect a field \"displacement\" to be defined."); return {}, false
	}

	model_params.soln_order = config_field.order

	model_params.dirichlet_displacements = make(map[fem.Boundary_ID]fem.Vec3)
	model_params.materials = make(map[fem.Section_ID]models.Elastic_Material_Int)
	model_params.variational_bcs = make(map[fem.Boundary_ID][]models.Elastic_BC_Int)
	model_params.sources = make(map[fem.Section_ID][]models.Elastic_Source_Int)
	model_params.ics = make(map[fem.Section_ID]fem.Vec3)

	for bnd_name, bnd_configs in config_field.boundaries {
		id, exists := mesh.boundary_names[bnd_name]
		if !exists {
			log.errorf("boundary %s is not declared on the mesh.", bnd_name); return {}, false
		}

		variational := make([dynamic]models.Elastic_BC_Int)
		defer if len(variational) != 0 {model_params.variational_bcs[id] = variational[:]}


		// these BC's are not additive, would be wrong to define a isothermal, insulated, convective boundary.
		had_fixed: bool
		had_free: bool

		for bnd_config in bnd_configs {
			bnd_type := property_get(string, bnd_config, "type") or_return
			switch bnd_type {
			case "fixed":
				had_fixed = true
				model_params.dirichlet_displacements[id] = property_get(fem.Vec3, bnd_config, "displacement") or_return
			case "free":
				had_free = true
			case "traction":
				data := new(Traction_Data)
				data.t = property_get(fem.Vec3, bnd_config, "traction") or_return
				append(&variational, models.Elastic_BC_Int {
					data = data,
					procedure = proc(
						ctx: fem.Element_Context,
						point_start, point_end: int,
						current_time: f64,
						current_soln: []fem.Vec3,
						data: rawptr,
						out: models.Elastic_BC,
					) {
						info := cast(^Traction_Data)data
						slice.fill(out.t, info.t)
					},
				})
			case:
				log.errorf("%s is not a recognized boundary for displacement.", bnd_type)
			}
		}
		if (had_fixed && had_free) || (had_fixed && (len(variational) > 0)) {
			log.errorf("Boundary %s was declared fixed no other boundaries may be applied.", bnd_name)
			return {}, false
		}

		if had_free && len(variational) > 0 {
			log.errorf("Boundary %s was declared free no other boundaries may be applied.", bnd_name)
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

			E := property_get(f64, mat, "elastic_modulus") or_return
			nu := property_get(f64, mat, "poissons_ratio") or_return

			data.C = isotropic_constitutive_tensor(E, nu)

			model_params.materials[id] = {
				data = data,
				procedure = proc(
					ctx: fem.Element_Context,
					current_strain: []fem.Voigt6,
					data: rawptr,
					out: models.Elastic_Material,
				) {
					info := cast(^Constant_Material_Data)data
					slice.fill(out.constitutive_tensor, info.C)
					for qp in 0 ..< len(ctx.points) {
						fem.voigt_gemv(info.C, current_strain[qp], &out.stress[qp], 1, 0)
					}

				},
			}

		case:
			log.errorf("%s is not a recognized material type.", mat_type); return {}, false
		}

		model_params.sources[id] = make([]models.Elastic_Source_Int, len(section_config.sources))
		for source, i in section_config.sources {
			type := property_get(string, source, "type") or_return
			switch type {
			case "constant":
				data := new(Constant_Source_Data)
				data.F = property_get(fem.Vec3, source, "load") or_return
				model_params.sources[id][i] = {
					data = data,
					procedure = proc(
						ctx: fem.Element_Context,
						current_time: f64,
						current_soln: []fem.Vec3,
						data: rawptr,
						out: models.Elastic_Source,
					) {
						info := cast(^Constant_Source_Data)data
						for &entry in out.F {entry += info.F}
					},
				}
			case:
				log.errorf("%s is not a recognized source term", type); return {}, false
			}
		}
	}

	isotropic_constitutive_tensor :: proc(E, nu: f64) -> fem.Voigt6x6 {
		lambda := E * nu / ((1 + nu) * (1 - 2 * nu))
		mu := E / (2 * (1 + nu))

		C: fem.Voigt6x6

		C[fem.voigt6x6_idx(0, 0)] = lambda + 2 * mu
		C[fem.voigt6x6_idx(1, 1)] = lambda + 2 * mu
		C[fem.voigt6x6_idx(2, 2)] = lambda + 2 * mu
		C[fem.voigt6x6_idx(0, 1)] = lambda; C[fem.voigt6x6_idx(1, 0)] = lambda
		C[fem.voigt6x6_idx(0, 2)] = lambda; C[fem.voigt6x6_idx(2, 0)] = lambda
		C[fem.voigt6x6_idx(1, 2)] = lambda; C[fem.voigt6x6_idx(2, 1)] = lambda
		C[fem.voigt6x6_idx(3, 3)] = mu
		C[fem.voigt6x6_idx(4, 4)] = mu
		C[fem.voigt6x6_idx(5, 5)] = mu

		return C
	}

	Constant_Material_Data :: struct {
		C: fem.Voigt6x6,
	}

	Constant_Source_Data :: struct {
		F: fem.Vec3,
	}

	Traction_Data :: struct {
		t: fem.Vec3,
	}


	return models.model_elasticity(model_params), true
}
