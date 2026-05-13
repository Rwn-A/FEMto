package inc_flow

import fem "../../fem"
import cfg "../../cfg"
import "core:slice"
import "core:log"

builtin_register :: proc(reg: ^cfg.Plugin_Registry) {
    reg["builtin:no_slip"] = no_slip
    reg["builtin:inflow"] = inflow
    reg["builtin:outflow"] = outflow
    reg["builtin:constant_material"] = constant_material
}

builtin_deregister :: proc(reg: ^cfg.Plugin_Registry) {
	delete_key(reg, "builtin:no_slip")
	delete_key(reg, "builtin:inflow")
	delete_key(reg, "builtin:outflow")
	delete_key(reg, "builtin:constant_material")
}


no_slip :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	ic: V_Fixed_Int
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []fem.Vec3) {
		slice.zero(out)
	}

	add_constraint_v(ctx, ic) or_return

	return true
}

inflow :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
    Schema :: struct {
        velocity: fem.Vec3 `validate:"required"`,
    }

    data := new(Schema)
    cfg.unmarshal(tbl, data) or_return

	ic: V_Fixed_Int
	ic.data = data
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []fem.Vec3) {
	    info := cast(^Schema)data
		slice.fill(out, info.velocity)
	}

	add_constraint_v(ctx, ic) or_return

	return true
}


outflow :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
    ic: P_Fixed_Int
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
	    slice.zero(out)
	}

	add_constraint_p(ctx, ic) or_return

	return true
}

constant_material :: proc(tbl: ^cfg.Table, ctx: cfg.Plugin_Context) -> bool {
	Constant_Material_Data :: struct {
		rho, mu: f64
	}

	Schema :: struct {
		viscosity:           f64 `validate:"required,gt=0"`,
		density:                f64 `validate:"required,gt=0"`,
	}

	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(Constant_Material_Data)
	data.rho = s.density
	data.mu = s.viscosity

	mat: Material_Int
	mat.data = data
	mat.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Material) {
		info := cast(^Constant_Material_Data)data
		slice.fill(out.density, info.rho)
		slice.fill(out.dynamic_viscosity, info.mu)
	}

	add_material(ctx, mat) or_return
	return true
}


add_constraint_v :: proc(ctx: cfg.Plugin_Context, c: V_Fixed_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_boundary in p.v_fixed_bcs {
		log.errorf("Cannot apply multiple constraints to boundary %s", ctx.current_bnd_name)
		return false
	}
	p.v_fixed_bcs[ctx.current_boundary] = c
	assert(ok, "Attemped to registry a conduction property where plugin context was not for conduction.")
	return true
}

add_constraint_p :: proc(ctx: cfg.Plugin_Context, c: P_Fixed_Int) -> bool {
	if ctx.allowed_kind != .BC {
		log.errorf("Cannot use a boundary condition in a non `boundary` block.")
		return false
	}
	p, ok := ctx.problem_data.(^Problem_Data)
	if ctx.current_boundary in p.p_fixed_bcs {
		log.errorf("Cannot apply multiple constraints to boundary %s", ctx.current_bnd_name)
		return false
	}
	p.p_fixed_bcs[ctx.current_boundary] = c
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