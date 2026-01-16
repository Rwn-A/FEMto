// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package fem

import "core:math"

Quadrature_Groups :: bit_set[Sample_Group]{.Integration_1, .Integration_3, .Facet_Integration_1, .Facet_Integration_3}

dV :: #force_inline proc(ctx: ^Sample_Context) -> f64 {
	return QUADRATURE_WEIGHTS[ctx.element.type][ctx.s.group][ctx.s.index] * sample_scale(ctx)
}

dS :: #force_inline proc(ctx: ^Sample_Context) -> f64 {
	return QUADRATURE_WEIGHTS[ctx.element.type][ctx.fs.group][ctx.fs.index] * sample_facet_scale(ctx)
}

// TODO: unify tables into this structure
// Quadrature_Rule :: struct {
//     points: []Vec3,
//     weights: []f64
//     facet_indices: []int,
// }

REC_ROOT_3 :: 1 / math.SQRT_THREE

POINT_QUADRATURE :: #partial [Sample_Group][]Vec3 {
	.Integration_1       = {{0, 0, 0}},
	.Integration_3       = {{0, 0, 0}},
	.Facet_Integration_1 = {},
	.Facet_Integration_3 = {},
}

LINE_QUADRATURE :: #partial [Sample_Group][]Vec3 {
	.Integration_1       = {{0, 0, 0}},
	.Integration_3       = {{-REC_ROOT_3, 0, 0}, {REC_ROOT_3, 0, 0}},
	.Facet_Integration_1 = {{-1, 0, 0}, {1, 0, 0}},
	.Facet_Integration_3 = {{-1, 0, 0}, {1, 0, 0}},
}

TRIANGLE_QUADRATURE :: #partial [Sample_Group][]Vec3 {
	.Integration_1       = {{1.0 / 3.0, 1.0 / 3.0, 0}},
	.Integration_3       = {{1.0 / 6.0, 1.0 / 6.0, 0}, {2.0 / 3.0, 1.0 / 6.0, 0}, {1.0 / 6.0, 2.0 / 3.0, 0}},
	.Facet_Integration_1 = {{0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}},
	.Facet_Integration_3 = {
		{0.5 - REC_ROOT_3 / 6.0, 0, 0},
		{0.5 + REC_ROOT_3 / 6.0, 0, 0},
		{0.5 + REC_ROOT_3 / 6.0, 0.5 - REC_ROOT_3 / 6.0, 0},
		{0.5 - REC_ROOT_3 / 6.0, 0.5 + REC_ROOT_3 / 6.0, 0},
		{0, 0.5 - REC_ROOT_3 / 6.0, 0},
		{0, 0.5 + REC_ROOT_3 / 6.0, 0},
	},
}

QUADRILATERAL_QUADRATURE :: #partial [Sample_Group][]Vec3 {
	.Integration_1       = {{0, 0, 0}},
	.Integration_3       = {
		{-REC_ROOT_3, -REC_ROOT_3, 0},
		{REC_ROOT_3, -REC_ROOT_3, 0},
		{REC_ROOT_3, REC_ROOT_3, 0},
		{-REC_ROOT_3, REC_ROOT_3, 0},
	},
	.Facet_Integration_1 = {{0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}},
	.Facet_Integration_3 = {
		{-REC_ROOT_3, -1, 0},
		{REC_ROOT_3, -1, 0},
		{1, -REC_ROOT_3, 0},
		{1, REC_ROOT_3, 0},
		{-REC_ROOT_3, 1, 0},
		{REC_ROOT_3, 1, 0},
		{-1, -REC_ROOT_3, 0},
		{-1, REC_ROOT_3, 0},
	},
}

TETRAHEDRON_QUADRATURE :: #partial [Sample_Group][]Vec3 {
	.Integration_1       = {{0.25, 0.25, 0.25}},
	.Integration_3       = {
		{0.138196601125011, 0.138196601125011, 0.138196601125011},
		{0.585410196624968, 0.138196601125011, 0.138196601125011},
		{0.138196601125011, 0.585410196624968, 0.138196601125011},
		{0.138196601125011, 0.138196601125011, 0.585410196624968},
	},
	.Facet_Integration_1 = {
		{1.0 / 3.0, 1.0 / 3.0, 0},
		{1.0 / 3.0, 0, 1.0 / 3.0},
		{0, 1.0 / 3.0, 1.0 / 3.0},
		{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
	},
	.Facet_Integration_3 = {
		{1.0 / 6.0, 1.0 / 6.0, 0},
		{2.0 / 3.0, 1.0 / 6.0, 0},
		{1.0 / 6.0, 2.0 / 3.0, 0},
		{1.0 / 6.0, 0, 1.0 / 6.0},
		{2.0 / 3.0, 0, 1.0 / 6.0},
		{1.0 / 6.0, 0, 2.0 / 3.0},
		{0, 1.0 / 6.0, 1.0 / 6.0},
		{0, 2.0 / 3.0, 1.0 / 6.0},
		{0, 1.0 / 6.0, 2.0 / 3.0},
		{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
		{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
		{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
	},
}

@(rodata)
QUADRATURE_WEIGHTS := #partial [Element_Type][Sample_Group][]f64 {
	.Point = #partial{.Integration_1 = {1}, .Integration_3 = {1}},
	.Line = #partial{
		.Integration_1 = {2},
		.Integration_3 = {1, 1},
		.Facet_Integration_1 = {1, 1},
		.Facet_Integration_3 = {1, 1},
	},
	.Triangle = #partial{
		.Integration_1 = {0.5},
		.Integration_3 = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
		.Facet_Integration_1 = {1.0, 1.0, 1.0},
		.Facet_Integration_3 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
	},
	.Quadrilateral = #partial{
		.Integration_1 = {4},
		.Integration_3 = {1, 1, 1, 1},
		.Facet_Integration_1 = {2, 2, 2, 2},
		.Facet_Integration_3 = {1, 1, 1, 1, 1, 1, 1, 1},
	},
	.Tetrahedron = #partial{
		.Integration_1 = {1.0 / 6.0},
		.Integration_3 = {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0},
		.Facet_Integration_1 = {0.5, 0.5, 0.5, 0.5},
		.Facet_Integration_3 = {
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
			1.0 / 6.0,
		},
	},
}
