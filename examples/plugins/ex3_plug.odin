//  example plugin
// built with: odin build .\examples\plugins\ex3_plug.odin -file -out:./.build/my_plug.dll -build-mode:shared
package main

import "core:math"

import "core:fmt"

import "../../src/app/cfg"
import "../../src/fem"
import "../../src/app/conduction"
import "../../src/app/linear_elasticity"

@(export)
register :: proc "c" (reg: ^cfg.Plugin_Registry) {
    reg["my_plug:custom_source"] = custom_source
}


custom_source :: proc(tbl: ^cfg.Table,  ctx: cfg.Plugin_Context) -> (ok: bool) {
    Plug_Schema :: struct {
        frequency: f64 `validate:"required,gt=0"`,
        amplitude: f64 `validate:"required"`,
    }

    fmt.println("Hello from my plugin!")


    Plug_Params :: struct{
        amplitude: f64,
        angular_frequency: f64
    }

    s: Plug_Schema
    cfg.unmarshal(tbl, &s) or_return

    source: conduction.Source_Int

    data := new(Plug_Params)
    data.angular_frequency = s.frequency * 2 * math.PI
    data.amplitude = s.amplitude
    source.data = data

    source.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, temperature: []f64, out: conduction.Source) {
        info := cast(^Plug_Params)data
        // always += because multiple sources might be active.
        for &Q in out.Q { Q += info.amplitude * math.sin(info.angular_frequency * time) }
    }

    conduction.add_source(ctx, source)

    return true
}
