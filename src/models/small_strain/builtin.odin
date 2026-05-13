package small_strain

import fem "../../fem"
import cfg "../../cfg"
import "core:log"
import "core:slice"

add_variational :: proc(ctx: cfg.Plugin_Context, v: Variational_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	assert(ok, "Attempted to registry a elasticity property where plugin context was not for elasticity.")
	if append(&p.variational_bcs[ctx.current_boundary], v) == 0 {
		log.errorf("config: Too many boundaries declared on boundary %s.", ctx.current_bnd_name)
		return false
	}
	return true
}

add_constraint :: proc(ctx: cfg.Plugin_Context, c: Fixed_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_boundary in p.fixed_bcs {
		log.errorf("Cannot apply multiple constraints to boundary %s", ctx.current_bnd_name)
		return false
	}
	p.fixed_bcs[ctx.current_boundary] = c
	assert(ok, "Attempted to registry an elastic property where plugin context was not for elasticity.")
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
	assert(ok, "Attempted to registry an elastic property where plugin context was not for elasticity.")
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
	assert(ok, "Attempted to registry an elastic property where plugin context was not for elasticity.")
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
	assert(ok, "Attempted to registry a elastic property where plugin context was not for elasticity.")

	return true
}



isotropic_constitutive_tensor :: proc(E, nu: f64) -> fem.Voigt6x6 {
	lambda := E * nu / ((1 + nu) * (1 - 2 * nu))
	mu := E / (2 * (1 + nu))
	C: fem.Voigt6x6
	C[0, 0] = lambda + 2 * mu
	C[1, 1] = lambda + 2 * mu
	C[2, 2] = lambda + 2 * mu
	C[0, 1] = lambda; C[1, 0] = lambda
	C[0, 2] = lambda; C[2, 0] = lambda
	C[1, 2] = lambda; C[2, 1] = lambda
	C[3, 3] = mu
	C[4, 4] = mu
	C[5, 5] = mu
	return C
}

builtin_register :: proc(reg: ^cfg.Plugin_Registry) {
	reg["builtin:constant_material"] = constant_material
	reg["builtin:constant_source"] = constant_source
	reg["builtin:constant_ic"] = constant_ic
	reg["builtin:fixed"] = fixed
	reg["builtin:constant_traction"] = constant_traction
}

builtin_deregister :: proc(reg: ^cfg.Plugin_Registry) {
	delete_key(reg, "builtin:constant_material")
	delete_key(reg, "builtin:constant_source")
	delete_key(reg, "builtin:constant_ic")
	delete_key(reg, "builtin:fixed")
	delete_key(reg, "builtin:constant_traction")
}

constant_material :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		elastic_modulus: f64 `validate:"required,gt=0"`,
		poissons_ratio:  f64 `validate:"required,gte=0,lt=0.5"`,
		density:         f64 `validate:"required,gt=0"`,
	}
	Elastic_Constant_Material_Data :: struct {
    	C:   fem.Voigt6x6,
    	rho: f64,
    }

	s: Schema
	cfg.unmarshal(tbl, &s) or_return
	data := new(Elastic_Constant_Material_Data)
	data.C = isotropic_constitutive_tensor(s.elastic_modulus, s.poissons_ratio)
	data.rho = s.density
	mat: Material_Int
	mat.data = data
	mat.procedure = proc(
		mapped: fem.Mapped_Element,
		time: f64,
		current_strain: []fem.Voigt6,
		data: rawptr,
		out: Material,
	) {
		info := cast(^Elastic_Constant_Material_Data)data
		slice.fill(out.constitutive_tensor, info.C)
		slice.fill(out.density, info.rho)
		for qp in 0 ..< len(out.stress) {
			out.stress[qp] = info.C * current_strain[qp]
		}
	}
	add_material(ctx, mat) or_return
	return true
}

constant_ic :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		displacement: fem.Vec3 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return
	data := new(fem.Vec3)
	data^ = s.displacement
	ic: Initial_Condition_Int
	ic.data = data
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []fem.Vec3) {
		slice.fill(out, (cast(^fem.Vec3)data)^)
	}
	add_initial_condition(ctx, ic) or_return
	return true
}

constant_source :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		body_force: fem.Vec3 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return
	data := new(fem.Vec3)
	data^ = s.body_force
	src: Source_Int
	src.data = data
	src.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Source) {
		for &v in out.body {v += (cast(^fem.Vec3)data)^}
	}
	add_source(ctx, src) or_return
	return true
}

fixed :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		displacement: fem.Vec3 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return
	data := new(fem.Vec3)
	data^ = s.displacement
	f: Fixed_Int
	f.data = data
	f.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []fem.Vec3) {
		slice.fill(out, (cast(^fem.Vec3)data)^)
	}
	add_constraint(ctx, f) or_return
	return true
}

constant_traction :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Schema :: struct {
		traction: fem.Vec3 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return
	data := new(fem.Vec3)
	data^ = s.traction
	v: Variational_Int
	v.data = data
	v.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Variational) {
		for &t in out.t {t += (cast(^fem.Vec3)data)^}
	}
	add_variational(ctx, v) or_return
	return true
}