package linear_elasticity

import fem "../../fem"
import cfg "../cfg"
import "core:log"
import "core:slice"

Elastic_Constant_Material_Data :: struct {
	C:   fem.Voigt6x6,
	rho: f64,
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

builtin_register :: proc(ctx: ^Plug_Context) {
	ctx.materials["builtin:constant_material"] = constant_material
	ctx.sources["builtin:constant_source"] = constant_source
	ctx.ics["builtin:constant_ic"] = constant_ic
	ctx.bcs["builtin:fixed"] = fixed
	ctx.bcs["builtin:constant_traction"] = constant_traction
}

constant_material :: proc(tbl: ^cfg.Table) -> (mat: Material_Int, ok: bool) {
	Schema :: struct {
		elastic_modulus: f64 `validate:"required,gt=0"`,
		poissons_ratio:  f64 `validate:"required,gte=0,lt=0.5"`,
		density:         f64 `validate:"required,gt=0"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(Elastic_Constant_Material_Data)
	data.C = isotropic_constitutive_tensor(s.elastic_modulus, s.poissons_ratio)
	data.rho = s.density

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
	return mat, true
}

constant_ic :: proc(tbl: ^cfg.Table) -> (ic: Initial_Condition_Int, ok: bool) {
	Schema :: struct {
		displacement: fem.Vec3 `validate:"required"`,
	}

	data := new(fem.Vec3)

	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data^ = s.displacement

	ic.data = data
	ic.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: []fem.Vec3) {
		d := cast(^fem.Vec3)data
		slice.fill(out, d^)
	}
	return ic, true
}

constant_source :: proc(tbl: ^cfg.Table) -> (src: Source_Int, ok: bool) {
	Schema :: struct {
		body_force: fem.Vec3 `validate:"required"`,
	}
	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(fem.Vec3)
	data^ = s.body_force

	src.data = data
	src.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Source) {
		for &v in out.body {
			v += (cast(^fem.Vec3)data)^
		}
	}
	return src, true
}

fixed :: proc(tbl: ^cfg.Table) -> (bc: BC_Int, ok: bool) {
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
		d := cast(^fem.Vec3)data
		slice.fill(out, d^)
	}
	return f, true
}

constant_traction :: proc(tbl: ^cfg.Table) -> (bc: BC_Int, ok: bool) {
	Schema :: struct {
		traction: fem.Vec3 `validate:"required"`,
	}

	s: Schema
	cfg.unmarshal(tbl, &s) or_return

	data := new(fem.Vec3)
	data^ = s.traction

	b: Variational_Int
	b.data = data
	b.procedure = proc(mapped: fem.Mapped_Element, time: f64, data: rawptr, out: Variational_BC){
		for &t in out.t {t += (cast(^fem.Vec3)data)^}
	}

	return b, true

}