// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
Subcells split an element into smaller sub elements, this is used for visualization.

Split_Quadratic for example, will create enough sub-elements that a linear approximation within each element will
accurately represent a quadratically accurate solution over the whole element.
*/
package fem

Subcell_Rule :: enum {
	Split_Linear,
	Split_Quadratic,
}

SUBCELL_RULE_TO_REF_ID := [Element_Type][Subcell_Rule]Ref_Context_ID{}

setup_subcell_rules :: proc(allocator := context.allocator) {
	context.allocator = allocator
	for rule in Subcell_Rule {
		SUBCELL_RULE_TO_REF_ID[.Point][rule] = register_reference_context(.Point, POINT_SUBCELL[rule])
		SUBCELL_RULE_TO_REF_ID[.Line][rule] = register_reference_context(.Line, LINE_SUBCELL[rule])
		SUBCELL_RULE_TO_REF_ID[.Triangle][rule] = register_reference_context(.Triangle, TRIANGLE_SUBCELL[rule])
		SUBCELL_RULE_TO_REF_ID[.Quadrilateral][rule] = register_reference_context(.Quadrilateral, QUADRILATERAL_SUBCELL[rule])
		SUBCELL_RULE_TO_REF_ID[.Tetrahedron][rule] = register_reference_context(.Tetrahedron, TETRAHEDRON_SUBCELL[rule])

	}
}

subcell_for :: proc(
	element: Mesh_Element,
	rule: Subcell_Rule,
	required_geometry: Geometry_Options,
	allocator := context.allocator,
) -> (
	ctx: Element_Context,
	count: int,
) {
	assert(required_geometry & Facet_Options == {}, "Subcell rules cannot supply facet geometry.")
	ctx = element_create_context(element, SUBCELL_RULE_TO_REF_ID[element.type][rule], required_geometry, allocator)
	return ctx, len(ctx.points)
}

@(rodata)
POINT_SUBCELL := [Subcell_Rule][]Vec3 {
	.Split_Linear    = {{0, 0, 0}},
	.Split_Quadratic = {{0, 0, 0}},
}

@(rodata)
LINE_SUBCELL := [Subcell_Rule][]Vec3 {
	.Split_Linear    = {{-1, 0, 0}, {1, 0, 0}},
	.Split_Quadratic = {{-1, 0, 0}, {0, 0, 0}, {1, 0, 0}},
}

@(rodata)
TRIANGLE_SUBCELL := [Subcell_Rule][]Vec3 {
	.Split_Linear    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}},
	.Split_Quadratic = {{0, 0, 0}, {0.5, 0, 0}, {1, 0, 0}, {0, 0.5, 0}, {0.5, 0.5, 0}, {0, 1, 0}},
}

@(rodata)
QUADRILATERAL_SUBCELL := [Subcell_Rule][]Vec3 {
	.Split_Linear    = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}},
	.Split_Quadratic = {{-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}},
}

@(rodata)
TETRAHEDRON_SUBCELL := [Subcell_Rule][]Vec3 {
	.Split_Linear    = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
	.Split_Quadratic = {
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


SUBCELL_CONNECTIVITY := #partial [Element_Type][Subcell_Rule][][]int {
	.Line = #partial{.Split_Linear = {{0, 1}}, .Split_Quadratic = {{0, 1}, {1, 2}}},
	.Quadrilateral = #partial{
		.Split_Linear = {{0, 1, 2, 3}},
		.Split_Quadratic = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}},
	},
	.Triangle = #partial{.Split_Linear = {{0, 1, 2}}, .Split_Quadratic = {{0, 3, 1}, {1, 4, 2}, {2, 5, 0}, {3, 4, 5}}},
	.Tetrahedron = #partial{
		.Split_Linear = {{0, 1, 2, 3}},
		.Split_Quadratic = {
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
