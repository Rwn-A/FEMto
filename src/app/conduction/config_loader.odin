package conduction

import "core:log"
import "core:slice"


import "../../fem"
import "../config"

import "../../../vendor/toml"

Field_Schema :: struct {
	order:              fem.Basis_Order,
	basis_family:       fem.Basis_Family,
	time_scheme:        fem.Time_Scheme,
	initial_conditions: ^toml.Table,
}

Model_Config :: struct {
	params:      Model_Parameters,
	bd:          fem.Basis_Descriptor,
	time_scheme: fem.Time_Scheme,
	ics:         map[fem.Section_ID]Initial_Condition_Int,
}

SUPPORTED_FAMILIES := bit_set[fem.Basis_Family]{.Lagrange}
SUPPORTED_TIME_SCHEMES := bit_set[fem.Time_Scheme]{.BE, .BDF2}


load_model_config :: proc(
	schema: config.Schema,
	cfg: config.General_Config,
	allocator := context.allocator,
) -> (
	model_cfg: Model_Config,
	ok: bool,
) {
	log.info("Loading parameters for conduction model.")
	context.allocator = allocator

	//field

	field_cfg := Field_Schema {
		order        = .Linear,
		basis_family = .Lagrange,
		time_scheme  = .BE,
	}

	model_cfg.ics = make(map[fem.Section_ID]Initial_Condition_Int)
	model_cfg.params.materials = make(map[fem.Section_ID]Material_Int)
	model_cfg.params.isothermal_bcs = make(map[fem.Boundary_ID]Isothermal_Int)


	if temperature_field, exists := config.table_get_opt(^toml.Table, schema.fields, "temperature"); exists {
		config.unmarshal(temperature_field, &field_cfg) or_return

		for section_name, id in cfg.mesh.section_names {
			ic_tbl := config.table_get_opt(^toml.Table, field_cfg.initial_conditions, section_name) or_continue

			form := config.table_get(string, ic_tbl, "form") or_return

			switch form {
			case "constant":
				model_cfg.ics[id] = constant_ic(ic_tbl) or_return
			case:
				log.errorf("config: unknown form %s for initial condition.", form)
				return {}, false
			}
		}
	}

	config.valid_option(SUPPORTED_FAMILIES, field_cfg.basis_family) or_return
	config.valid_option(SUPPORTED_TIME_SCHEMES, field_cfg.time_scheme) or_return

	model_cfg.bd = {
		family = field_cfg.basis_family,
		order  = field_cfg.order,
	}
	model_cfg.time_scheme = field_cfg.time_scheme

	// material

	for section_name, id in cfg.mesh.section_names {
		mat_tbl := config.table_get(^toml.Table, schema.materials, section_name) or_return

		form := config.table_get(string, mat_tbl, "form") or_return

		switch form {
		case "constant":
			model_cfg.params.materials[id] = constant_material(mat_tbl) or_return
		case:
			log.errorf("config: unknown form %s for material.", form)
			return {}, false
		}
	}

	// sources

	for source in schema.sources {
		source_tbl, is_tbl := source.(^toml.Table)
		if !is_tbl {
			log.error("Each source entry was expected to be a table")
			return {}, false
		}

		form := config.table_get(string, source_tbl, "form") or_return

		source: Source_Int
		switch form {
		case "constant":
			source = constant_source(source_tbl) or_return
		case:
			log.errorf("config: unknown form %s for source", form)
			return {}, false
		}

		regions, given := config.table_get_opt(^toml.List, source_tbl, "region")
		for section_name, id in cfg.mesh.section_names {
			source_builder := make([dynamic]Source_Int)
			defer model_cfg.params.sources[id] = source_builder[:]

			if !given {
				append(&source_builder, source)
			}else{
				for region in regions {
					region_str, is_str := region.(string)
					if !is_str{
						log.error("region in source `region` must be a string")
						return {}, false
					}
					if region_str not_in cfg.mesh.section_names {
						log.errorf("region %s was not defined on the mesh", region_str)
						return {}, false
					}
					if region_str == section_name {
						append(&source_builder, source)
					}
				}
			}


		}


	}

	// boundaries

	for boundary_name, id in cfg.mesh.boundary_names {
		bc_arr := config.table_get_opt(^toml.List, schema.boundaries, boundary_name) or_continue

		variational_builder := make([dynamic]Variational_Int)
		defer model_cfg.params.variational_bcs[id] = variational_builder[:]

		// some boundaries cant be stacked (isothermal, adiabatic)
		had_exclusive: bool

		for bc in bc_arr {
			if had_exclusive {
				log.errorf(
					"Boundary %s can not have multiple conditions because one exclusive condition has already been applied.",
					boundary_name,
				)
				return {}, false
			}
			bc_tbl, is_tbl := bc.(^toml.Table)
			if !is_tbl {
				log.errorf("Boundary condition on boundary %s was expected to be a table", boundary_name)
				return {}, false
			}

			kind := config.table_get(string, bc_tbl, "kind") or_return

			switch kind {
			case "isothermal":
				model_cfg.params.isothermal_bcs[id] = isothermal(bc_tbl) or_return
				had_exclusive = true
			case "convective":
				append(&variational_builder, convective(bc_tbl) or_return)
			case "fixed_flux":
				append(&variational_builder, fixed_flux(bc_tbl) or_return)
			case "radiative":
				append(&variational_builder, radiative(bc_tbl) or_return)
			case "adiabatic":
				had_exclusive = true
			case:
				log.errorf("config: unknown kind %s for boundary condition on boundary %s.", kind, boundary_name)
				return {}, false

			}
		}
	}


	return model_cfg, true
}

constant_material :: proc(tbl: ^toml.Table) -> (mat: Material_Int, ok: bool) {
	Constant_Material_Data :: struct {
		k, rho, cp: f64,
	}
	data := new(Constant_Material_Data)

	data.k = config.table_get(f64, tbl, "conductivity") or_return
	data.rho = config.table_get(f64, tbl, "density") or_return
	data.cp = config.table_get(f64, tbl, "specific_heat_capacity") or_return

	mat.data = data
	mat.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Material) {
		info := cast(^Constant_Material_Data)data
		slice.fill(out.k, info.k); slice.fill(out.rho, info.rho); slice.fill(out.cp, info.cp)
	}

	return mat, true
}

constant_ic :: proc(tbl: ^toml.Table) -> (ic: Initial_Condition_Int, ok: bool) {
	temp := config.table_get(f64, tbl, "temperature") or_return

	ic.data = transmute(rawptr)temp
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
		slice.fill(out, transmute(f64)data)
	}

	return ic, true
}

constant_source :: proc(tbl: ^toml.Table) -> (source: Source_Int, ok: bool) {
	heat := config.table_get(f64, tbl, "heat") or_return
	source.data = transmute(rawptr)heat
	source.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Source) {
		for &Q in out.Q {Q += transmute(f64)data}
	}

	return source, true
}

isothermal :: proc(tbl: ^toml.Table) -> (bc: Isothermal_Int, ok: bool) {
	temp := config.table_get(f64, tbl, "temperature") or_return

	bc.data = transmute(rawptr)temp
	bc.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
		slice.fill(out, transmute(f64)data)
	}

	return bc, true
}

convective :: proc(tbl: ^toml.Table) -> (bc: Variational_Int, ok: bool) {
	Convective_Data :: struct {
		h:       f64,
		ambient: f64,
	}

	data := new(Convective_Data)
	data.h = config.table_get(f64, tbl, "convective_coefficient") or_return
	data.ambient = config.table_get(f64, tbl, "ambient_temp") or_return

	return {
			data = data,
			procedure = proc(
				mapped: fem.Mapped_Element,
				time: f64,
				data: rawptr,
				temperature: []f64,
				out: Variational_BC,
			) {
				info := cast(^Convective_Data)data
				for &h in out.h {h += info.h}
				for &T_amb in out.h {T_amb += info.ambient}
			},
		},
		true
}

radiative :: proc(tbl: ^toml.Table) -> (bc: Variational_Int, ok: bool) {
	Radiative_Data :: struct {
		emissivity:       f64,
		ambient:          f64,
		stefan_boltzmann: f64,
	}

	data := new(Radiative_Data)
	data.emissivity = config.table_get(f64, tbl, "emissivity") or_return
	data.ambient = config.table_get(f64, tbl, "ambient_temp") or_return
	data.stefan_boltzmann = config.table_get(f64, tbl, "stefan_boltzmann") or_return

	return {
			data = data,
			procedure = proc(
				mapped: fem.Mapped_Element,
				time: f64,
				data: rawptr,
				temperature: []f64,
				out: Variational_BC,
			) {
				info := cast(^Radiative_Data)data
				sigma := info.stefan_boltzmann
				for i in 0 ..<len(temperature) {
					T := temperature[i]
					T_amb := info.ambient
					out.q[i] = info.emissivity * sigma * (T_amb * T_amb * T_amb * T_amb - T * T * T * T)
					out.d_q[i] = -4 * info.emissivity * sigma * T * T * T
				}
			},
		},
		true

}

fixed_flux :: proc(tbl: ^toml.Table) -> (bc: Variational_Int, ok: bool) {

	flux := config.table_get(f64, tbl, "flux") or_return

	return {
			data = transmute(rawptr)flux,
			procedure = proc(
				mapped: fem.Mapped_Element,
				time: f64,
				data: rawptr,
				temperature: []f64,
				out: Variational_BC,
			) {
				flux := transmute(f64)data
				for &q in out.q {q += flux}
			},
		},
		true
}
