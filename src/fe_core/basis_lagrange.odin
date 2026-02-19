// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package fem

Basis_LS :: struct {
	using info: Basis_Info,
}
Basis_LV :: distinct Basis_LS

ls_create :: proc(element: Element, order: Order) -> Basis_LS {
	raw := RAW_LAGRANGE_BASIS[element.type][order]
	return {order = order, _support = raw.support, arity = raw.arity, geometry_required = {.Jac_Inv}, components = 1}
}

lv_create :: proc(element: Element, order: Order) -> Basis_LV {
	ls := ls_create(element, order)
	ls.components = AMBIENT_DIM
	ls.arity = ls.arity * ls.components
	return cast(Basis_LV)ls
}

ls_value :: proc(ls: Basis_LS, ctx: Element_Context, point, basis: int) -> f64 {
	return ctx.basis_lagrange_table[ls.order].values[point][basis]
}

ls_gradient :: proc(ls: Basis_LS, ctx: Element_Context, point, basis: int) -> Vec3 {
	return ctx.basis_lagrange_table[ls.order].gradients[point][basis] * ctx_inverse_jacobian(ctx, point)
}

lv_value :: proc(lv: Basis_LV, ctx: Element_Context, point, basis: int) -> (r: Vec3) {
	sb, cmpnt := basis_decompose_component(lv, basis)
	r[cmpnt] = ctx.basis_lagrange_table[lv.order].values[point][sb]
	return r
}

lv_gradient :: proc(lv: Basis_LV, ctx: Element_Context, point, basis: int) -> (m: Mat3) {
	sb, cmpnt := basis_decompose_component(lv, basis)

	g := ls_gradient(cast(Basis_LS)lv, ctx, point, basis)

	m[0, cmpnt] = g[0]
	m[1, cmpnt] = g[1]
	m[2, cmpnt] = g[2]

	return m
}

// precomputed in ref space available to element context.
// indexed via point, basis index.
Lagrange_Context_Table :: struct {
	values:    [][]f64,
	gradients: [][]Vec3,
}

@(private)
populate_lagrange_context_table :: proc(
	element_type: Element_Type,
	points: []Vec3,
	allocator := context.allocator,
) -> (
	table: [Order]Lagrange_Context_Table,
) {
	for order in Order {
		basis := RAW_LAGRANGE_BASIS[element_type][order]
		table[order].values = make([][]f64, len(points))
		table[order].gradients = make([][]Vec3, len(points))

		for point, point_idx in points {
			table[order].values[point_idx] = make([]f64, basis.arity)
			table[order].gradients[point_idx] = make([]Vec3, basis.arity)
			for basis_index in 0 ..< basis.arity {
				val, grad := basis.reference_data(basis_index, point)
				table[order].values[point_idx][basis_index] = val
				table[order].gradients[point_idx][basis_index] = grad
			}
		}
	}
	return table
}


Raw_Lagrange_Info :: struct {
	support:        []Basis_Support,
	arity:          int,
	reference_data: proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3),
}

RAW_LAGRANGE_BASIS := [Element_Type][Order]Raw_Lagrange_Info {
	.Point         = LAGRANGE_POINT_BASIS,
	.Line          = LAGRANGE_LINE_BASIS,
	.Quadrilateral = LAGRANGE_QUAD_BASIS,
	.Triangle      = LAGRANGE_TRI_BASIS,
	.Tetrahedron   = LAGRANGE_TET_BASIS,
}

raw_lagrange_value :: proc(type: Element_Type, order: Order, idx: int, ref_point: Vec3) -> f64 {
	v, _ := RAW_LAGRANGE_BASIS[type][order].reference_data(idx, ref_point)
	return v
}

raw_lagrange_gradient :: proc(type: Element_Type, order: Order, idx: int, ref_point: Vec3) -> Vec3 {
	_, g := RAW_LAGRANGE_BASIS[type][order].reference_data(idx, ref_point)
	return g
}

@(rodata)
LAGRANGE_POINT_BASIS := [Order]Raw_Lagrange_Info {
	.Linear = {
		support = {{.D0, 0, 0, 1}},
		arity = 1,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {return 1, 0},
	},
	.Quadratic = {
		support = {{.D0, 0, 0, 1}},
		arity = 1,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {return 1, 0},
	},
}

@(rodata)
LAGRANGE_LINE_BASIS := [Order]Raw_Lagrange_Info {
	.Linear = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}},
		arity = 2,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			switch idx {
			case 0:
				return (1.0 - r.x) / 2.0, Vec3{-0.5, 0, 0}
			case 1:
				return (1.0 + r.x) / 2.0, Vec3{0.5, 0, 0}
			case:
				unreachable()
			}
		},
	},
	.Quadratic = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D1, 0, 0, 1}},
		arity = 3,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			switch idx {
			case 0:
				return -r.x * (1.0 - r.x) / 2.0, Vec3{-(1.0 - 2.0 * r.x) / 2.0, 0, 0}
			case 1:
				return r.x * (1.0 + r.x) / 2.0, Vec3{(1.0 + 2.0 * r.x) / 2.0, 0, 0}
			case 2:
				return 1.0 - r.x * r.x, Vec3{-2.0 * r.x, 0, 0}
			case:
				unreachable()
			}
		},
	},
}

@(rodata)
LAGRANGE_QUAD_BASIS := [Order]Raw_Lagrange_Info {
	.Linear = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}, {.D0, 3, 0, 1}},
		arity = 4,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y
			switch idx {
			case 0:
				return ((1 - x) * (1 - y)) * 0.25, Vec3{-(1 - y) * 0.25, -(1 - x) * 0.25, 0}
			case 1:
				return ((1 + x) * (1 - y)) * 0.25, Vec3{(1 - y) * 0.25, -(1 + x) * 0.25, 0}
			case 2:
				return ((1 + x) * (1 + y)) * 0.25, Vec3{(1 + y) * 0.25, (1 + x) * 0.25, 0}
			case 3:
				return ((1 - x) * (1 + y)) * 0.25, Vec3{-(1 + y) * 0.25, (1 - x) * 0.25, 0}
			case:
				unreachable()
			}
		},
	},
	.Quadratic = {
		support = {
			{.D0, 0, 0, 1},
			{.D0, 1, 0, 1},
			{.D0, 2, 0, 1},
			{.D0, 3, 0, 1},
			{.D1, 0, 0, 1},
			{.D1, 1, 0, 1},
			{.D1, 2, 0, 1},
			{.D1, 3, 0, 1},
			{.D2, 0, 0, 1},
		},
		arity = 9,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y

			Lx := [3]f64{0.5 * x * (x - 1.0), 0.5 * x * (x + 1.0), 1.0 - x * x}
			Ly := [3]f64{0.5 * y * (y - 1.0), 0.5 * y * (y + 1.0), 1.0 - y * y}
			dLx := [3]f64{x - 0.5, x + 0.5, -2.0 * x}
			dLy := [3]f64{y - 0.5, y + 0.5, -2.0 * y}

			switch idx {
			case 0:
				return Lx[0] * Ly[0], Vec3{dLx[0] * Ly[0], Lx[0] * dLy[0], 0}
			case 1:
				return Lx[1] * Ly[0], Vec3{dLx[1] * Ly[0], Lx[1] * dLy[0], 0}
			case 2:
				return Lx[1] * Ly[1], Vec3{dLx[1] * Ly[1], Lx[1] * dLy[1], 0}
			case 3:
				return Lx[0] * Ly[1], Vec3{dLx[0] * Ly[1], Lx[0] * dLy[1], 0}
			case 4:
				return Lx[2] * Ly[0], Vec3{dLx[2] * Ly[0], Lx[2] * dLy[0], 0}
			case 5:
				return Lx[1] * Ly[2], Vec3{dLx[1] * Ly[2], Lx[1] * dLy[2], 0}
			case 6:
				return Lx[2] * Ly[1], Vec3{dLx[2] * Ly[1], Lx[2] * dLy[1], 0}
			case 7:
				return Lx[0] * Ly[2], Vec3{dLx[0] * Ly[2], Lx[0] * dLy[2], 0}
			case 8:
				return Lx[2] * Ly[2], Vec3{dLx[2] * Ly[2], Lx[2] * dLy[2], 0}
			case:
				unreachable()
			}
		},
	},
}

@(rodata)
LAGRANGE_TRI_BASIS := [Order]Raw_Lagrange_Info {
	.Linear = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}},
		arity = 3,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y
			switch idx {
			case 0:
				return 1.0 - x - y, Vec3{-1, -1, 0}
			case 1:
				return x, Vec3{1, 0, 0}
			case 2:
				return y, Vec3{0, 1, 0}
			case:
				unreachable()
			}
		},
	},
	.Quadratic = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}, {.D1, 0, 0, 1}, {.D1, 1, 0, 1}, {.D1, 2, 0, 1}},
		arity = 6,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y

			l0 := 1.0 - x - y
			l1 := x
			l2 := y

			dl0 := Vec3{-1, -1, 0}
			dl1 := Vec3{1, 0, 0}
			dl2 := Vec3{0, 1, 0}

			switch idx {
			case 0:
				return l0 * (2 * l0 - 1), (4 * l0 - 1) * dl0
			case 1:
				return l1 * (2 * l1 - 1), (4 * l1 - 1) * dl1
			case 2:
				return l2 * (2 * l2 - 1), (4 * l2 - 1) * dl2
			case 3:
				return 4 * l0 * l1, 4 * (l0 * dl1 + l1 * dl0)
			case 4:
				return 4 * l1 * l2, 4 * (l1 * dl2 + l2 * dl1)
			case 5:
				return 4 * l2 * l0, 4 * (l2 * dl0 + l0 * dl2)
			case:
				unreachable()
			}
		},
	},
}

@(rodata)
LAGRANGE_TET_BASIS := [Order]Raw_Lagrange_Info {
	.Linear = {
		support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}, {.D0, 3, 0, 1}},
		arity = 4,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y
			z := r.z
			switch idx {
			case 0:
				return 1.0 - x - y - z, Vec3{-1, -1, -1}
			case 1:
				return x, Vec3{1, 0, 0}
			case 2:
				return y, Vec3{0, 1, 0}
			case 3:
				return z, Vec3{0, 0, 1}
			case:
				unreachable()
			}
		},
	},
	.Quadratic = {
		support = {
			{.D0, 0, 0, 1},
			{.D0, 1, 0, 1},
			{.D0, 2, 0, 1},
			{.D0, 3, 0, 1},
			{.D1, 0, 0, 1},
			{.D1, 1, 0, 1},
			{.D1, 2, 0, 1},
			{.D1, 3, 0, 1},
			{.D1, 4, 0, 1},
			{.D1, 5, 0, 1},
		},
		arity = 10,
		reference_data = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
			x := r.x
			y := r.y
			z := r.z

			l0 := 1.0 - x - y - z
			l1 := x
			l2 := y
			l3 := z

			dl0 := Vec3{-1, -1, -1}
			dl1 := Vec3{1, 0, 0}
			dl2 := Vec3{0, 1, 0}
			dl3 := Vec3{0, 0, 1}

			switch idx {
			case 0:
				return l0 * (2 * l0 - 1), (4 * l0 - 1) * dl0
			case 1:
				return l1 * (2 * l1 - 1), (4 * l1 - 1) * dl1
			case 2:
				return l2 * (2 * l2 - 1), (4 * l2 - 1) * dl2
			case 3:
				return l3 * (2 * l3 - 1), (4 * l3 - 1) * dl3
			case 4:
				return 4 * l0 * l1, 4 * (l0 * dl1 + l1 * dl0)
			case 5:
				return 4 * l0 * l2, 4 * (l0 * dl2 + l2 * dl0)
			case 6:
				return 4 * l0 * l3, 4 * (l0 * dl3 + l3 * dl0)
			case 7:
				return 4 * l1 * l2, 4 * (l1 * dl2 + l2 * dl1)
			case 8:
				return 4 * l1 * l3, 4 * (l1 * dl3 + l3 * dl1)
			case 9:
				return 4 * l2 * l3, 4 * (l2 * dl3 + l3 * dl2)
			case:
				unreachable()
			}
		},
	},
}
