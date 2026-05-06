package conduction

import "core:math"
import "core:slice"
import cfg "../cfg"
import fem "../../fem"

Constant_Material_Data :: struct {
    k, rho, cp: f64,
}

Ramped_ISO :: struct {
    i_temp, f_temp: f64,
    slope:          f64,
}

Pulsed_ISO :: struct {
    angular_frequency: f64,
    amplitude:         f64,
    offset:            f64,
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


builtin_register :: proc (ctx: ^Plug_Context) {
    ctx.materials["builtin:constant_material"] = constant_material
    ctx.sources["builtin:constant_source"]     = constant_source
    ctx.ics["builtin:constant_ic"]             = constant_ic
    ctx.bcs["builtin:isothermal"]              = isothermal
    ctx.bcs["builtin:ramped_isothermal"]       = ramped_iso
    ctx.bcs["builtin:pulsed_isothermal"]       = pulsed_iso
    ctx.bcs["builtin:convective"]              = convective
    ctx.bcs["builtin:radiative"]               = radiative
    ctx.bcs["builtin:fixed_flux"]              = fixed_flux
}

constant_material :: proc(tbl: ^cfg.Table) -> (mat: Material_Int, ok: bool) {
    Schema :: struct {
        conductivity:           f64 `validate:"required,gt=0"`,
        density:                f64 `validate:"required,gt=0"`,
        specific_heat_capacity: f64 `validate:"required,gt=0"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    data := new(Constant_Material_Data)
    data.k   = s.conductivity
    data.rho  = s.density
    data.cp   = s.specific_heat_capacity

    mat.data = data
    mat.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Material) {
        info := cast(^Constant_Material_Data)data
        slice.fill(out.k, info.k)
        slice.fill(out.rho, info.rho)
        slice.fill(out.cp, info.cp)
    }
    return mat, true
}

constant_ic :: proc(tbl: ^cfg.Table) -> (ic: Initial_Condition_Int, ok: bool) {
    Schema :: struct {
        temperature: f64 `validate:"required"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    ic.data = transmute(rawptr)s.temperature
    ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
        slice.fill(out, transmute(f64)data)
    }
    return ic, true
}

constant_source :: proc(tbl: ^cfg.Table) -> (source: Source_Int, ok: bool) {
    Schema :: struct {
        heat: f64 `validate:"required"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    source.data = transmute(rawptr)s.heat
    source.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Source) {
        for &Q in out.Q { Q += transmute(f64)data }
    }
    return source, true
}

isothermal :: proc(tbl: ^cfg.Table) -> (bc_u: BC_Int, ok: bool) {
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
    return bc, true
}

ramped_iso :: proc(tbl: ^cfg.Table) -> (bc_u: BC_Int, ok: bool) {
    Schema :: struct {
        initial_temperature: f64 `validate:"required"`,
        final_temperature:   f64 `validate:"required"`,
        duration:            f64 `validate:"required,gt=0"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    data := new(Ramped_ISO)
    data.i_temp = s.initial_temperature
    data.f_temp = s.final_temperature
    data.slope  = (s.final_temperature - s.initial_temperature) / s.duration

    bc: Isothermal_Int

    bc.data = data
    bc.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
        info := cast(^Ramped_ISO)data
        temp := info.slope * time + info.i_temp
        slice.fill(out, temp if temp <= info.f_temp else info.f_temp)
    }
    return bc, true
}

pulsed_iso :: proc(tbl: ^cfg.Table) -> (bc_u: BC_Int, ok: bool) {
    Schema :: struct {
        frequency: f64 `validate:"required,gt=0"`,
        amplitude: f64 `validate:"required"`,
        offset:    f64 `validate:"required"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    data := new(Pulsed_ISO)
    data.angular_frequency = s.frequency * 2 * math.PI
    data.amplitude         = s.amplitude
    data.offset            = s.offset

    bc: Isothermal_Int

    bc.data = data
    bc.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []f64) {
        info := cast(^Pulsed_ISO)data
        slice.fill(out, info.offset + info.amplitude * math.sin(info.angular_frequency * time))
    }
    return bc, true
}

convective :: proc(tbl: ^cfg.Table) -> (bc: BC_Int, ok: bool) {
    Schema :: struct {
        convective_coefficient: f64 `validate:"required,gt=0"`,
        ambient_temp:           f64 `validate:"required"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    data := new(Convective_Data)
    data.h       = s.convective_coefficient
    data.ambient = s.ambient_temp


    return Variational_Int{
        data = data,
        procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Variational_BC) {
            info := cast(^Convective_Data)data
            for &h in out.h { h += info.h }
            for &T_amb in out.T_amb { T_amb += info.ambient }
        },
    }, true
}

radiative :: proc(tbl: ^cfg.Table) -> (bc: BC_Int, ok: bool) {
    Schema :: struct {
        emissivity:       f64 `validate:"required,gte=0,lte=1"`,
        ambient_temp:     f64 `validate:"required"`,
        stefan_boltzmann: f64 `validate:"required,gt=0"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    data := new(Radiative_Data)
    data.emissivity       = s.emissivity
    data.ambient          = s.ambient_temp
    data.stefan_boltzmann = s.stefan_boltzmann

    return Variational_Int{
        data = data,
        procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Variational_BC) {
            info := cast(^Radiative_Data)data
            sigma := info.stefan_boltzmann
            for i in 0 ..< len(temperature) {
                T     := temperature[i]
                T_amb := info.ambient
                out.q[i]   += info.emissivity * sigma * (T_amb*T_amb*T_amb*T_amb - T*T*T*T)
                out.d_q[i] += -4 * info.emissivity * sigma * T*T*T
            }
        },
    }, true
}

fixed_flux :: proc(tbl: ^cfg.Table) -> (bc: BC_Int, ok: bool) {
    Schema :: struct {
        flux: f64 `validate:"required"`,
    }
    s: Schema
    cfg.unmarshal(tbl, &s) or_return

    return Variational_Int{
        data = transmute(rawptr)s.flux,
        procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: Variational_BC) {
            flux := transmute(f64)data
            for &q in out.q { q += flux }
        },
    }, true
}