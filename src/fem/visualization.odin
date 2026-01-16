// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package fem

Visualization_Groups :: bit_set[Sample_Group]{.Visualization_Linear, .Visualization_Quadratic}

MAX_OUTPUT_FIELD_COMPONENTS :: 9

Visualization_Proc :: #type proc(ctx: ^Sample_Context, data: rawptr) -> [MAX_OUTPUT_FIELD_COMPONENTS]f64

Visualization_Field :: struct {
	friendly_name:  string,
	components:     int,
	data:           rawptr,
	value_provider: Visualization_Proc,
}

POINT_VISUALIZATION :: #partial [Sample_Group][]Vec3 {
	.Visualization_Linear    = {{0, 0, 0}},
	.Visualization_Quadratic = {{0, 0, 0}},
}

LINE_VISUALIZATION :: #partial [Sample_Group][]Vec3 {
	.Visualization_Linear    = {{-1, 0, 0}, {1, 0, 0}},
	.Visualization_Quadratic = {{-1, 0, 0}, {0, 0, 0}, {1, 0, 0}},
}

TRIANGLE_VISUALIZATION :: #partial [Sample_Group][]Vec3 {
	.Visualization_Linear    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}},
	.Visualization_Quadratic = {{0, 0, 0}, {0.5, 0, 0}, {1, 0, 0}, {0, 0.5, 0}, {0.5, 0.5, 0}, {0, 1, 0}},
}

QUADRILATERAL_VISUALIZATION :: #partial [Sample_Group][]Vec3 {
	.Visualization_Linear    = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}},
	.Visualization_Quadratic = {
		{-1, -1, 0},
		{0, -1, 0},
		{1, -1, 0},
		{-1, 0, 0},
		{0, 0, 0},
		{1, 0, 0},
		{-1, 1, 0},
		{0, 1, 0},
		{1, 1, 0},
	},
}

TETRAHEDRON_VISUALIZATION :: #partial [Sample_Group][]Vec3 {
	.Visualization_Linear    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
	.Visualization_Quadratic = {
		{0, 0, 0},
		{0.5, 0, 0},
		{1, 0, 0},
		{0, 0.5, 0},
		{0.5, 0.5, 0},
		{0, 1, 0},
		{0, 0, 0.5},
		{0.5, 0, 0.5},
		{0, 0.5, 0.5},
		{0, 0, 1},
	},
}

@(rodata)
VISUALIZATION_CONNECTIVITY := #partial [Element_Type][Sample_Group][][]int {
	.Line = #partial{.Visualization_Linear = {{0, 1}}, .Visualization_Quadratic = {{0, 1}, {1, 2}}},
	.Quadrilateral = #partial{
		.Visualization_Linear = {{0, 1, 2, 3}},
		.Visualization_Quadratic = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}},
	},
	.Triangle = #partial{
		.Visualization_Linear = {{0, 1, 2}},
		.Visualization_Quadratic = {{0, 3, 1}, {1, 4, 2}, {2, 5, 0}, {3, 4, 5}},
	},
	.Tetrahedron = #partial{
		.Visualization_Linear = {{0, 1, 2, 3}},
		.Visualization_Quadratic = {
			{0, 1, 3, 6},
			{2, 4, 1, 7},
			{5, 3, 4, 8},
			{9, 6, 8, 7},
			{1, 4, 3, 7},
			{3, 4, 7, 8},
			{1, 3, 6, 7},
			{3, 6, 7, 8},
		},
	},
}
