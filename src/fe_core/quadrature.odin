// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package fem

import "core:math"

import "core:log"

Quadrature_Rule :: enum {
	Quad_1,
	Quad_3,
	Facet_Quad_1,
	Facet_Quad_3,
}

Facet_Rules :: bit_set[Quadrature_Rule]{.Facet_Quad_1, .Facet_Quad_3}

QUADRATURE_RULE_TO_CTX_ID := [Element_Type][Quadrature_Rule]Ref_Context_ID{}

Quadrature_Context :: struct {
	using ctx: Element_Context,
	weights:   []f64,
}

infer_quadrature_rule :: proc(element: Element, max_basis_order: Order, domain: enum {
		Interior,
		Facet,
	}) -> Quadrature_Rule {
	dimension := element_dim(element.type) if domain == .Interior else element_facet_dim(element.type)

	// +1 because enums start at 0
	basis_int_order := int(max_basis_order) + 1
	element_int_order := int(element.order) + 1
	dimension_int_order := int(dimension) + 1

	//0 is okay for linear element order, because linear elements dont modify the needed order
	required_order := (basis_int_order * 2) + (element_int_order - 1) * (dimension_int_order)

	switch {
	case required_order <= 1:
		return .Quad_1 if domain == .Interior else .Facet_Quad_1
	case required_order <= 3:
		return .Quad_3 if domain == .Interior else .Facet_Quad_3
	case:
		panic("Simulation does not support accurate enough quadrature for the requested integrand.")
	}
}

// Wrapper over `Element_Context` for quadrature specifically. required_geometry should be queried from your basis.
quadrature_for :: proc(
	element: Mesh_Element,
	rule: Quadrature_Rule,
	required_geometry: Geometry_Options,
	allocator := context.allocator,
) -> (
	qctx: Quadrature_Context,
	num_points: int,
) {
	context.allocator = allocator
	quad_geometry: Geometry_Options = {.Jac_Det} if rule not_in Facet_Rules else {.Facet_Jac_Det, .Facet_Normal}

	qctx.weights = QUADRATURE_RULES[element.type][rule].weights
	qctx.ctx = element_create_context(element, QUADRATURE_RULE_TO_CTX_ID[element.type][rule], quad_geometry | required_geometry)

	return qctx, len(qctx.weights)
}

// TODO: quadrature needs to ensure facet indices are ordered, and we can probably cache start, end.
facet_quad :: proc(qctx: Quadrature_Context, facet: int) -> (start: int, end: int) {
	start = -1
	for facet_index, i in qctx.facet_indices {
		if facet_index == facet {start = i}
		if start != -1 && facet_index != facet {end = i}
	}
	if end == 0 {end = len(qctx.facet_indices)}
	return start, end
}

dV :: #force_inline proc(qctx: Quadrature_Context, point: int) -> f64 {
	return qctx.weights[point] * ctx_jacobian_determinant(qctx, point)
}

dS :: #force_inline proc(qctx: Quadrature_Context, point: int) -> f64 {
	return qctx.weights[point] * ctx_facet_jacobian_determinant(qctx, point)
}

setup_quadrature_rules :: proc(allocator := context.allocator) {
	context.allocator = allocator

	//generate facet rules and register
	for element in Element_Type {
		for rule in Quadrature_Rule {
			if rule not_in Facet_Rules {
				id := register_reference_context(element, QUADRATURE_RULES[element][rule].points)
				QUADRATURE_RULE_TO_CTX_ID[element][rule] = id
				continue
			}

			// no facet rules for points.
			if element == .Point {continue}

			point_builder := make([dynamic]Vec3)
			weight_builder := make([dynamic]f64)
			indices_builder := make([dynamic]int)

			for facet_index in 0 ..< element_num_facets(element) {
				volume_rule: Quadrature_Rule
				#partial switch rule {
				case .Facet_Quad_1:
					volume_rule = .Quad_1
				case .Facet_Quad_3:
					volume_rule = .Quad_3
				case:
					unreachable()
				}
				facet_rule_table := QUADRATURE_RULES[element_facet_type(element, facet_index)][volume_rule]

				for qp, index in soa_zip(point = facet_rule_table.points, weight = facet_rule_table.weights) {
					append(&point_builder, facet_reference_to_volume_reference(element, qp.point, facet_index))
					append(&weight_builder, qp.weight)
					append(&indices_builder, facet_index)
				}
			}

			rule_table := Quadrature_Rule_Table {
				points  = point_builder[:],
				weights = weight_builder[:],
			}

			QUADRATURE_RULES[element][rule] = rule_table

			id := register_reference_context(element, rule_table.points, facet_indices = indices_builder[:])
			QUADRATURE_RULE_TO_CTX_ID[element][rule] = id
		}
	}
}

@(private)
Quadrature_Rule_Table :: struct {
	points:  []Vec3,
	weights: []f64,
}

QUADRATURE_RULES := #partial [Element_Type][Quadrature_Rule]Quadrature_Rule_Table {
	.Point         = POINT_QUADRATURE,
	.Line          = LINE_QUADRATURE,
	.Triangle      = TRIANGLE_QUADRATURE,
	.Quadrilateral = QUADRILATERAL_QUADRATURE,
	.Tetrahedron   = TETRAHEDRON_QUADRATURE,
}


REC_ROOT_3 :: 1 / math.SQRT_THREE


POINT_QUADRATURE :: #partial [Quadrature_Rule]Quadrature_Rule_Table {
	.Quad_1 = {points = {{0, 0, 0}}, weights = {1}},
	.Quad_3 = {points = {{0, 0, 0}}, weights = {1}},
}

LINE_QUADRATURE :: #partial [Quadrature_Rule]Quadrature_Rule_Table {
	.Quad_1 = {points = {{0, 0, 0}}, weights = {2}},
	.Quad_3 = {points = {{-REC_ROOT_3, 0, 0}, {REC_ROOT_3, 0, 0}}, weights = {1, 1}},
}

QUADRILATERAL_QUADRATURE :: #partial [Quadrature_Rule]Quadrature_Rule_Table {
	.Quad_1 = {points = {{0, 0, 0}}, weights = {4}},
	.Quad_3 = {
		points = {
			{-REC_ROOT_3, -REC_ROOT_3, 0},
			{REC_ROOT_3, -REC_ROOT_3, 0},
			{REC_ROOT_3, REC_ROOT_3, 0},
			{-REC_ROOT_3, REC_ROOT_3, 0},
		},
		weights = {1, 1, 1, 1},
	},
}

TRIANGLE_QUADRATURE :: #partial [Quadrature_Rule]Quadrature_Rule_Table {
	.Quad_1 = {points = {{1.0 / 3.0, 1.0 / 3.0, 0}}, weights = {0.5}},
	.Quad_3 = {
		points = {{1.0 / 6.0, 1.0 / 6.0, 0}, {2.0 / 3.0, 1.0 / 6.0, 0}, {1.0 / 6.0, 2.0 / 3.0, 0}},
		weights = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
	},
}

TETRAHEDRON_QUADRATURE :: #partial [Quadrature_Rule]Quadrature_Rule_Table {
	.Quad_1 = {points = {{0.25, 0.25, 0.25}}, weights = {1.0 / 6.0}},
	.Quad_3 = {
		points = {
			{0.138196601125011, 0.138196601125011, 0.138196601125011},
			{0.585410196624968, 0.138196601125011, 0.138196601125011},
			{0.138196601125011, 0.585410196624968, 0.138196601125011},
			{0.138196601125011, 0.138196601125011, 0.585410196624968},
		},
		weights = {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0},
	},
}
