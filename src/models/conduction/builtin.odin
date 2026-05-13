package conduction

import fem "../../fem"
import cfg "../../cfg"
import "core:slice"
import "core:log"


builtin_register :: proc(reg: ^cfg.Plugin_Registry) {
	reg["builtin:constant_material"] = constant_material
	reg["builtin:constant_source"] = constant_source
	reg["builtin:constant_ic"] = constant_ic
	reg["builtin:isothermal"] = isothermal
	reg["builtin:convective"] = convective
	reg["builtin:radiative"] = radiative
	reg["builtin:fixed_flux"] = fixed_flux
}

// this is so multiple models that are loaded can have conflicitng names
// once a model is done loading its bc's etc it unloads it builtins.
// user plugins arent tied to one model so they are always loaded.
builtin_deregister :: proc(reg: ^cfg.Plugin_Registry) {
	delete_key(reg, "builtin:constant_material")
	delete_key(reg, "builtin:constant_source")
	delete_key(reg, "builtin:constant_ic")
	delete_key(reg, "builtin:isothermal")
	delete_key(reg, "builtin:convective")
	delete_key(reg, "builtin:radiative")
	delete_key(reg, "builtin:fixed_flux")
}

add_variational :: proc(ctx: cfg.Plugin_Context, v: Variational_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")
	if append(&p.variational_bcs[ctx.current_boundary], v) == 0 {
		log.errorf("config: Too many boundaries declared on boundary %s.", ctx.current_bnd_name)
		return false
	}
	return true
}

add_constraint :: proc(ctx: cfg.Plugin_Context, c: Isothermal_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_boundary in p.isothermal_bcs {
		log.errorf("Cannot apply multiple constraints to boundary %s", ctx.current_bnd_name)
		return false
	}
	p.isothermal_bcs[ctx.current_boundary] = c
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")
	return true
}

add_source :: proc(ctx: cfg.Plugin_Context, s: Source_Int) -> bool {
	if ctx.allowed_kind != .Source {
		log.errorf("Cannot use a source term in a non `source` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if append(&p.sources[ctx.current_section], s) == 0 {
		log.errorf("config: Too many sources declared on section %s.", ctx.current_section_name)
		return false
	}
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")
	return true
}

add_material :: proc(ctx: cfg.Plugin_Context, m: Material_Int) -> bool {
	if ctx.allowed_kind != .Material {
		log.errorf("Cannot use a material in a non `material` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_section in p.materials {
		log.errorf("Cannot apply multiple materials to one region %s", ctx.current_section_name)
		return false
	}
	p.materials[ctx.current_section] = m
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")
	return true
}

add_initial_condition :: proc(ctx: cfg.Plugin_Context, ic: Initial_Condition_Int) -> bool {
	if ctx.allowed_kind != .IC {
		log.errorf("Cannot use a initial condition in a non `initial_condition block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_section in p.ics {
		log.errorf("Cannot apply multiple initial conditions to one region %s", ctx.current_section_name)
		return false
	}
	p.ics[ctx.current_section] = ic
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")

	return true
}



constant_material :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Constant_Material_Data :: struct {
		k, rho, cp: f64,
	}

	Schema :: struct {
		conductivity:           f64 `validate:"required,gt=0"`,
		density:                f64 `validate:"required,gt=0"`,
		specific_heat_capacity: f64 `validate:"required,gt=0"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(Constant_Material_Data)
	data.k = s.conductivity
	data.rho = s.density
	data.cp = s.specific_heat_capacity

	mat: Material_Int
	mat.data = data
	mat.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Material) {
		info := cast(^Constant_Material_Data)data
		slice.fill(out.k, info.k)
		slice.fill(out.rho, info.rho)
		slice.fill(out.cp, info.cp)
	}

	add_material(ctx, mat) or_return
	return true
}

constant_ic :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		temperature: f64 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	ic: Initial_Condition_Int
	ic.data = transmute(rawptr)s.temperature
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
		slice.fill(out, transmute(f64)data)
	}

	add_initial_condition(ctx, ic) or_return

	return true
}

constant_source :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		heat: f64 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	source: Source_Int
	source.data = transmute(rawptr)s.heat
	source.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Source) {
		for &Q in out.Q {Q += transmute(f64)data}
	}

	add_source(ctx, source) or_return
	return true
}

isothermal :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		temperature: f64 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	bc: Isothermal_Int

	bc.data = transmute(rawptr)s.temperature
	bc.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
		slice.fill(out, transmute(f64)data)
	}
	add_constraint(ctx, bc)
	return true
}

convective :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Convective_Data :: struct {
		h:       f64,
		ambient: f64,
	}

	Schema :: struct {
		convective_coefficient: f64 `validate:"required,gt=0"`,
		ambient_temp:           f64 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(Convective_Data)
	data.h = s.convective_coefficient
	data.ambient = s.ambient_temp

	v := Variational_Int {
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
			for &T_amb in out.T_amb {T_amb += info.ambient}
		},
	}
	add_variational(ctx, v) or_return
	return true
}

radiative :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Radiative_Data :: struct {
		emissivity:       f64,
		ambient:          f64,
		stefan_boltzmann: f64,
	}

	Schema :: struct {
		emissivity:       f64 `validate:"required,gte=0,lte=1"`,
		ambient_temp:     f64 `validate:"required"`,
		stefan_boltzmann: f64 `validate:"required,gt=0"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(Radiative_Data)
	data.emissivity = s.emissivity
	data.ambient = s.ambient_temp
	data.stefan_boltzmann = s.stefan_boltzmann

	v := Variational_Int {
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
			for i in 0 ..< len(temperature) {
				T := temperature[i]
				T_amb := info.ambient
				out.q[i] += info.emissivity * sigma * (T_amb * T_amb * T_amb * T_amb - T * T * T * T)
				out.d_q[i] += -4 * info.emissivity * sigma * T * T * T
			}
		},
	}
	add_variational(ctx, v) or_return
	return true
}

fixed_flux :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		flux: f64 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	v: Variational_Int = {
		data = transmute(rawptr)s.flux,
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
	}
	add_variational(ctx, v) or_return

	return true

}
