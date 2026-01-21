// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
Basis function definitions and interface.

Compile-time dispatch for now, possible changes down the road.
*/
package fem


Basis_Support :: struct {
	entity_dim:                                           Dimension,
	entity_index, entity_basis_index, entity_basis_arity: int,
}

Lagrange_Scalar :: struct {
	support:   []Basis_Support,
	reference: [Sample_Group]Lagrange_Reference_Table,
	arity:     int,
}

Lagrange_Vector :: distinct Lagrange_Scalar

basis_create :: proc($T: typeid, element: Element, order: Order) -> (basis: T) {
	when T == Lagrange_Scalar || T == Lagrange_Vector {
		raw := RAW_LAGRANGE_BASIS[element.type][order]
		return {support = raw.support, arity = raw.arity, reference = LAGRANGE_REFERENCE_TABLES[element.type][order]}
	} else {
		#panic("Type T is not a valid basis type.")
	}
}


basis_ls_gradient :: proc(ls: Lagrange_Scalar, ctx: ^Sample_Context, basis: int) -> Vec3 {
	return ls.reference[ctx.s.group].gradients[ctx.s.index][basis] * sample_pullback(ctx)
}

basis_lv_gradient :: proc(ls: Lagrange_Scalar, ctx: ^Sample_Context, basis: int) -> (m: Mat3) {
	scalar_basis := basis / AMBIENT_DIM
    comp := basis % AMBIENT_DIM

    g := basis_ls_gradient(ls, ctx, scalar_basis)

    m[0, comp] = g[0]
    m[1, comp] = g[1]
    m[2, comp] = g[2]

    return m
}

basis_gradient :: proc {
	basis_ls_gradient,
}

basis_ls_value :: proc(ls: Lagrange_Scalar, ctx: ^Sample_Context, basis: int) -> f64 {
	return ls.reference[ctx.s.group].values[ctx.s.index][basis]
}

basis_lv_value :: proc(lv: Lagrange_Vector, ctx: ^Sample_Context, basis: int) -> (r: Vec3) {
	r[basis % AMBIENT_DIM] = lv.reference[ctx.s.group].values[ctx.s.index][basis / AMBIENT_DIM]
	return r
}

basis_value :: proc {
	basis_ls_value,
	basis_lv_value,
}


basis_ls_support :: proc(ls: Lagrange_Scalar, basis: int) -> Basis_Support {
	return ls.support[basis]
}

basis_lv_support :: proc(lv: Lagrange_Vector, basis: int) -> Basis_Support {
	scalar_support := lv.support[basis / AMBIENT_DIM]
	scalar_support.entity_basis_index = (scalar_support.entity_basis_index * AMBIENT_DIM) + (basis % AMBIENT_DIM)
	scalar_support.entity_basis_arity *= AMBIENT_DIM
	return scalar_support
}

basis_support :: proc {
	basis_ls_support,
	basis_lv_support,
}


Lagrange_Reference_Table :: struct {
	values:    [][]f64,
	gradients: [][]Vec3,
}

// filled at startup.
LAGRANGE_REFERENCE_TABLES := [Element_Type][Order][Sample_Group]Lagrange_Reference_Table{}


populate_basis_tables :: proc(allocator := context.allocator) {
	context.allocator = allocator
	for type in Element_Type {
		for order in Order {
			basis := RAW_LAGRANGE_BASIS[type][order]
			for ps, ps_id in SAMPLE_POINTS[type] {
				table := &LAGRANGE_REFERENCE_TABLES[type][order][ps_id]
				table.values = make([][]f64, len(ps))
				table.gradients = make([][]Vec3, len(ps))
				for point, point_idx in ps {
					table.values[point_idx] = make([]f64, basis.arity)
					table.gradients[point_idx] = make([]Vec3, basis.arity)
					for basis_index in 0 ..< basis.arity {
						val, grad := basis.reference_data(basis_index, point)
						table.values[point_idx][basis_index] = val
						table.gradients[point_idx][basis_index] = grad
					}
				}
			}
		}
	}
}

// raw basis

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

// Use for evalauting scalar lagrange basis at non-sample points.
raw_shape_value :: proc(type: Element_Type, order: Order, idx: int, ref_point: Vec3) -> f64 {
	v, _ := RAW_LAGRANGE_BASIS[type][order].reference_data(idx, ref_point)
	return v
}

// Use for evalauting scalar lagrange basis at non-sample points.
raw_shape_gradient :: proc(type: Element_Type, order: Order, idx: int, ref_point: Vec3) -> Vec3 {
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
