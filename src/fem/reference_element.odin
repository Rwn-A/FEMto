// Defines reference data for all supported elements and wraps topology queries.
// Only topology is likely to be used by user code directly.
// The rest of the reference data is wrapped into different systems
// with the mapped element structure and basis interfaces.
package fem

import "core:math"

REC_ROOT_3 :: 1 / math.SQRT_THREE
ROOT_3_5 :: 0.7745966692414834
MAX_FACETS :: 6

Element_Type :: enum {
	Point,
	Line,
	Triangle,
	Quadrilateral,
	Tetrahedron,
	Hexahedron,
}

Dimension :: enum {
	D0, // vertex
	D1, // edge
	D2, // face
	D3, // cell
}


Element_Order :: enum {
	Linear,
	Quadratic,
}

// might make sense down the road to seperate these,
// as we might have expect far higher order basis then we support for geometry.
Basis_Order :: Element_Order

Reference_Topology :: struct {
	dimension:           Dimension,
	// (order, dimension, which sub entity of that dimension) -> nodes for that sub entity.
	sub_entity_nodes:    [Element_Order][Dimension][][]int,
	facet_element_types: []Element_Type,
	facet_tangents:      [][]Vec3,
}

Quadrature_Rule :: enum {
	Quad_1,
	Quad_3,
	Quad_5,
}

Reference_Quadrature :: struct {
	points:  []Vec3,
	weights: []f64,
}

Reference_Subcell :: struct {
	points:       []Vec3,
	connectivity: [][]int,
}

Basis_Support :: struct {
	entity_dim:         Dimension, //which entity supports the basis
	entity_index:       int, // which local entity index on the element is this (vertex 0, face 1, etc.)
	entity_basis_index: int, // which basis on the specific entity is this.
	entity_basis_arity: int, // how many total dofs exist on this entity.
}

Reference_Basis_Lagrange :: struct {
	support:          []Basis_Support,
	reference_values: proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3),
	arity:            int,
}

Element_Reference_Data :: struct {
	topology:   Reference_Topology,
	quadrature: [Quadrature_Rule]Reference_Quadrature,
	subcell:    [Basis_Order]Reference_Subcell,
	lagrange:   [Basis_Order]Reference_Basis_Lagrange,
}

// table

@(rodata)
REFERENCE_ELEMENTS := [Element_Type]Element_Reference_Data {
	.Point         = POINT_DATA,
	.Line          = LINE_DATA,
	.Quadrilateral = QUADRILATERAL_DATA,
	.Triangle      = TRIANGLE_DATA,
	.Tetrahedron   = TETRAHEDRON_DATA,
	.Hexahedron    = HEXAHEDRON_DATA,
}


// topology access functions

element_dim :: proc(type: Element_Type) -> Dimension {
	return REFERENCE_ELEMENTS[type].topology.dimension
}

element_facet_type :: proc(parent_type: Element_Type, facet: int, loc := #caller_location) -> Element_Type {
	assert(parent_type != .Point, loc = loc)
	return REFERENCE_ELEMENTS[parent_type].topology.facet_element_types[facet]
}

element_facet_count :: proc(parent_type: Element_Type, loc := #caller_location) -> int {
	assert(parent_type != .Point, loc = loc)
	return len(REFERENCE_ELEMENTS[parent_type].topology.facet_element_types)
}

element_facet_tangents :: proc(parent_type: Element_Type, facet: int, loc := #caller_location) -> []Vec3 {
	assert(parent_type != .Point && parent_type != .Line, loc = loc)
	return REFERENCE_ELEMENTS[parent_type].topology.facet_tangents[facet]
}

// returns a list of which nodes on the given element make up the given entity.
// ex: element_sub_entity_nodes(.Triangle, .Linear, .D1, 2)
// gets all the nodes that make up the third edge (index 2) of the triangle.
element_sub_entity_nodes :: proc(
	type: Element_Type,
	order: Element_Order,
	entity_dimension: Dimension,
	entity: int,
	loc := #caller_location,
) -> []int {
	assert(type != .Point, loc = loc)
	assert(entity_dimension < element_dim(type), loc = loc)
	return REFERENCE_ELEMENTS[type].topology.sub_entity_nodes[order][entity_dimension][entity]
}

// returns the sub entity nodes for all entities on the element of that dimension.
element_sub_entity_nodes_all :: proc(
	type: Element_Type,
	order: Element_Order,
	entity_dimension: Dimension,
	loc := #caller_location,
) -> [][]int {
	assert(type != .Point, loc = loc)
	assert(entity_dimension < element_dim(type), loc = loc)
	return REFERENCE_ELEMENTS[type].topology.sub_entity_nodes[order][entity_dimension]
}

// Returns true if two sub-entities on the same reference element share a node.
// Nodes can be acquired with `element_sub_entity_nodes`.
// Iterates over node indices rather than using a precomputed table.
// Queries should be off hot paths anyway.
element_sub_entities_incident :: proc(a_nodes, b_nodes: []int) -> bool {
	for an in a_nodes {
		for bn in b_nodes {
			if an == bn {return true}
		}
	}
	return false
}

support_on_facet :: proc(element: Element, s: Basis_Support, facet_index: int) -> bool {
	facet_dim := element_dim(element_facet_type(element.type, facet_index))
	facet_nodes := element_sub_entity_nodes(element.type, element.order, facet_dim, facet_index)
	dof_support_nodes := element_sub_entity_nodes(element.type, element.order, s.entity_dim, s.entity_index)
	return element_sub_entities_incident(facet_nodes, dof_support_nodes)
}

// Moves a point from the facets reference space, to the volume elements reference space.
facet_reference_to_volume_reference :: proc(element_type: Element_Type, facet_point: Vec3, facet_index: int) -> Vec3 {
	assert(element_type != .Point)
	switch element_type {
	case .Point:
		unreachable()
	case .Line:
		return facet_index == 0 ? Vec3{-1, 0, 0} : Vec3{1, 0, 0}
	case .Triangle:
		switch facet_index {
		case 0:
			return Vec3{facet_point.x, 0, 0}
		case 1:
			return Vec3{1 - facet_point.x, facet_point.x, 0}
		case 2:
			return Vec3{0, 1 - facet_point.x, 0}
		case:
			unreachable()
		}
	case .Quadrilateral:
		switch facet_index {
		case 0:
			return Vec3{facet_point.x, -1, 0}
		case 1:
			return Vec3{1, facet_point.x, 0}
		case 2:
			return Vec3{facet_point.x, 1, 0}
		case 3:
			return Vec3{-1, facet_point.x, 0}
		case:
			unreachable()
		}
	case .Tetrahedron:
		switch facet_index {
		case 0:
			return Vec3{facet_point.x, facet_point.y, 0}
		case 1:
			return Vec3{facet_point.x, 0, facet_point.y}
		case 2:
			return Vec3{0, facet_point.x, facet_point.y}
		case 3:
			return Vec3{1 - facet_point.x - facet_point.y, facet_point.x, facet_point.y}
		case:
			unreachable()
		}
	case .Hexahedron:
		switch facet_index {
		case 0:
			return Vec3{facet_point.x, facet_point.y, -1}
		case 1:
			return Vec3{facet_point.x, facet_point.y, +1}
		case 2:
			return Vec3{facet_point.x, -1, facet_point.y}
		case 3:
			return Vec3{facet_point.x, +1, facet_point.y}
		case 4:
			return Vec3{-1, facet_point.x, facet_point.y}
		case 5:
			return Vec3{+1, facet_point.x, facet_point.y}
		case:
			unreachable()
		}
	case:
		unreachable()
	}
}


POINT_DATA :: Element_Reference_Data {
	topology = {dimension = .D0},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
		.Quad_1 = {points = {{0, 0, 0}}, weights = {1}},
		.Quad_3 = {points = {{0, 0, 0}}, weights = {1}},
		.Quad_5 = {points = {{0, 0, 0}}, weights = {1}},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {{.D0, 0, 0, 1}},
			arity = 1,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {return 1, 0},
		},
		.Quadratic = {
			support = {{.D0, 0, 0, 1}},
			arity = 1,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {return 1, 0},
		},
	},
}


LINE_DATA :: Element_Reference_Data {
	topology = {
		dimension = .D1,
		facet_element_types = {.Point, .Point},
		sub_entity_nodes = {.Linear = #partial{.D0 = {{0}, {1}}}, .Quadratic = #partial{.D0 = {{0}, {1}}}},
	},
	subcell = [Basis_Order]Reference_Subcell {
		.Linear = {points = {{-1, 0, 0}, {1, 0, 0}}, connectivity = {{0, 1}}},
		.Quadratic = {points = {{-1, 0, 0}, {0, 0, 0}, {1, 0, 0}}, connectivity = {{0, 1}, {1, 2}}},
	},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
		.Quad_1 = {points = {{0, 0, 0}}, weights = {2}},
		.Quad_3 = {points = {{-REC_ROOT_3, 0, 0}, {REC_ROOT_3, 0, 0}}, weights = {1, 1}},
		.Quad_5 = {points = {{-ROOT_3_5, 0, 0}, {0, 0, 0}, {ROOT_3_5, 0, 0}}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9}},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}},
			arity = 2,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
	},
}

QUADRILATERAL_DATA :: Element_Reference_Data {
	topology = {
		dimension = .D2,
		facet_element_types = {.Line, .Line, .Line, .Line},
		facet_tangents = {{{1, 0, 0}}, {{0, 1, 0}}, {{-1, 0, 0}}, {{0, -1, 0}}},
		sub_entity_nodes = {
			.Linear = #partial{.D0 = {{0}, {1}, {2}, {3}}, .D1 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}}},
			.Quadratic = #partial{.D0 = {{0}, {1}, {2}, {3}}, .D1 = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}}},
		},
	},
	subcell = [Basis_Order]Reference_Subcell {
		.Linear = {points = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}}, connectivity = {{0, 1, 2, 3}}},
		.Quadratic = {
			points = {
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
			connectivity = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}},
		},
	},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
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
		.Quad_5 = {
			points = {
				{-ROOT_3_5, -ROOT_3_5, 0},
				{0., -ROOT_3_5, 0},
				{ROOT_3_5, -ROOT_3_5, 0},
				{-ROOT_3_5, 0., 0},
				{0., 0., 0},
				{ROOT_3_5, 0., 0},
				{-ROOT_3_5, ROOT_3_5, 0},
				{0., ROOT_3_5, 0},
				{ROOT_3_5, ROOT_3_5, 0},
			},
			weights = {
				25. / 81.,
				40. / 81.,
				25. / 81.,
				40. / 81.,
				64. / 81.,
				40. / 81.,
				25. / 81.,
				40. / 81.,
				25. / 81.,
			},
		},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}, {.D0, 3, 0, 1}},
			arity = 4,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
	},
}

TRIANGLE_DATA :: Element_Reference_Data {
	topology = {
		dimension = .D2,
		facet_element_types = {.Line, .Line, .Line},
		facet_tangents = {{{0.5, 0, 0}}, {{-0.5, 0.5, 0}}, {{0, -0.5, 0}}}, // because lines arent (-1, 1)
		sub_entity_nodes = {
			.Linear = #partial{.D0 = {{0}, {1}, {2}}, .D1 = {{0, 1}, {1, 2}, {2, 0}}},
			.Quadratic = #partial{.D0 = {{0}, {1}, {2}}, .D1 = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}}},
		},
	},
	subcell = [Basis_Order]Reference_Subcell {
		.Linear = {points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, connectivity = {{0, 1, 2}}},
		.Quadratic = {
			points = {{0, 0, 0}, {0.5, 0, 0}, {1, 0, 0}, {0, 0.5, 0}, {0.5, 0.5, 0}, {0, 1, 0}},
			connectivity = {{0, 3, 1}, {1, 4, 2}, {2, 5, 0}, {3, 4, 5}},
		},
	},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
		.Quad_1 = {points = {{1.0 / 3.0, 1.0 / 3.0, 0}}, weights = {0.5}},
		.Quad_3 = {
			points = {{1.0 / 6.0, 1.0 / 6.0, 0}, {2.0 / 3.0, 1.0 / 6.0, 0}, {1.0 / 6.0, 2.0 / 3.0, 0}},
			weights = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
		},
		.Quad_5 = {
			points = {
				{0.091576213509771, 0.091576213509771, 0},
				{0.816847572980459, 0.091576213509771, 0},
				{0.091576213509771, 0.816847572980459, 0},
				{0.445948490915965, 0.108103018168070, 0},
				{0.108103018168070, 0.445948490915965, 0},
				{0.445948490915965, 0.445948490915965, 0},
			},
			weights = {
				0.054975871827661,
				0.054975871827661,
				0.054975871827661,
				0.111690794839005,
				0.111690794839005,
				0.111690794839005,
			},
		},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}},
			arity = 3,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
	},
}

TETRAHEDRON_DATA :: Element_Reference_Data {
	topology = {
		dimension = .D3,
		facet_element_types = {.Triangle, .Triangle, .Triangle, .Triangle},
		facet_tangents = {
			{{0, 1, 0}, {1, 0, 0}},
			{{1, 0, 0}, {0, 0, 1}},
			{{-1, 1, 0}, {-1, 0, 1}},
			{{0, 0, 1}, {0, 1, 0}},
		},
		sub_entity_nodes = {
			.Linear = #partial{
				.D0 = {{0}, {1}, {2}, {3}},
				.D1 = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}},
				.D2 = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}},
			},
			.Quadratic = #partial{
				.D0 = {{0}, {1}, {2}, {3}},
				.D1 = {{0, 1, 4}, {1, 2, 5}, {2, 0, 6}, {0, 3, 7}, {1, 3, 8}, {2, 3, 9}},
				.D2 = {{0, 1, 2, 4, 5, 6}, {0, 1, 3, 4, 8, 7}, {1, 2, 3, 5, 9, 8}, {2, 0, 3, 6, 7, 9}},
			},
		},
	},
	subcell = [Basis_Order]Reference_Subcell {
		.Linear = {points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, connectivity = {{0, 1, 2, 3}}},
		.Quadratic = {
			points = {
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
			connectivity = {
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
	},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
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
		.Quad_5 = {
			points = {
				{0.0927352503108912, 0.0927352503108912, 0.0927352503108912},
				{0.7217942490673264, 0.0927352503108912, 0.0927352503108912},
				{0.0927352503108912, 0.7217942490673264, 0.0927352503108912},
				{0.0927352503108912, 0.0927352503108912, 0.7217942490673264},
				{0.3108859192633006, 0.3108859192633006, 0.3108859192633006},
				{0.0673422421100982, 0.3108859192633006, 0.3108859192633006},
				{0.3108859192633006, 0.0673422421100982, 0.3108859192633006},
				{0.3108859192633006, 0.3108859192633006, 0.0673422421100982},
				{0.0455037041256495, 0.0455037041256495, 0.4544962958743505},
				{0.0455037041256495, 0.4544962958743505, 0.0455037041256495},
				{0.4544962958743505, 0.0455037041256495, 0.0455037041256495},
				{0.4544962958743505, 0.4544962958743505, 0.0455037041256495},
				{0.4544962958743505, 0.0455037041256495, 0.4544962958743505},
				{0.0455037041256495, 0.4544962958743505, 0.4544962958743505},
			},
			weights = {
				0.01224884051939365,
				0.01224884051939365,
				0.01224884051939365,
				0.01224884051939365,
				0.01878132095300263,
				0.01878132095300263,
				0.01878132095300263,
				0.01878132095300263,
				0.00709100346284690,
				0.00709100346284690,
				0.00709100346284690,
				0.00709100346284690,
				0.00709100346284690,
				0.00709100346284690,
			},
		},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {{.D0, 0, 0, 1}, {.D0, 1, 0, 1}, {.D0, 2, 0, 1}, {.D0, 3, 0, 1}},
			arity = 4,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
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
					return 4 * l1 * l2, 4 * (l1 * dl2 + l2 * dl1)
				case 6:
					return 4 * l0 * l2, 4 * (l0 * dl2 + l2 * dl0)
				case 7:
					return 4 * l0 * l3, 4 * (l0 * dl3 + l3 * dl0)
				case 8:
					return 4 * l1 * l3, 4 * (l1 * dl3 + l3 * dl1)
				case 9:
					return 4 * l2 * l3, 4 * (l2 * dl3 + l3 * dl2)
				case:
					unreachable()
				}
			},
		},
	},
}

HEXAHEDRON_DATA :: Element_Reference_Data {
	topology = {
		dimension = .D3,
		facet_element_types = {
			.Quadrilateral,
			.Quadrilateral,
			.Quadrilateral,
			.Quadrilateral,
			.Quadrilateral,
			.Quadrilateral,
		},
		facet_tangents = {
			{{0, 1, 0}, {1, 0, 0}},
			{{1, 0, 0}, {0, 1, 0}},
			{{1, 0, 0}, {0, 0, 1}},
			{{0, 0, 1}, {1, 0, 0}},
			{{0, 0, 1}, {0, 1, 0}},
			{{0, 1, 0}, {0, 0, 1}},
		},
		sub_entity_nodes = {
			.Linear = #partial{
				.D0 = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}},
				.D1 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}},
				.D2 = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}},
			},
			.Quadratic = #partial{
				.D0 = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}},
				.D1 = {
					{0, 1, 8},
					{1, 2, 9},
					{2, 3, 10},
					{3, 0, 11},
					{4, 5, 12},
					{5, 6, 13},
					{6, 7, 14},
					{7, 4, 15},
					{0, 4, 16},
					{1, 5, 17},
					{2, 6, 18},
					{3, 7, 19},
				},
				.D2 = {
					{0, 1, 2, 3, 8, 9, 10, 11, 20},
					{4, 5, 6, 7, 12, 13, 14, 15, 21},
					{0, 1, 5, 4, 8, 17, 12, 16, 22},
					{3, 2, 6, 7, 10, 18, 14, 19, 23},
					{0, 3, 7, 4, 11, 19, 15, 16, 24},
					{1, 2, 6, 5, 9, 18, 13, 17, 25},
				},
			},
		},
	},
	subcell = [Basis_Order]Reference_Subcell {
		.Linear = {
			points = {
				{-1, -1, -1},
				{1, -1, -1},
				{1, 1, -1},
				{-1, 1, -1},
				{-1, -1, 1},
				{1, -1, 1},
				{1, 1, 1},
				{-1, 1, 1},
			},
			connectivity = {{0, 1, 2, 3, 4, 5, 6, 7}},
		},
		.Quadratic = {
			points = {
				{-1, -1, -1},
				{1, -1, -1},
				{1, 1, -1},
				{-1, 1, -1},
				{-1, -1, 1},
				{1, -1, 1},
				{1, 1, 1},
				{-1, 1, 1},
				{0, -1, -1},
				{1, 0, -1},
				{0, 1, -1},
				{-1, 0, -1},
				{0, -1, 1},
				{1, 0, 1},
				{0, 1, 1},
				{-1, 0, 1},
				{-1, -1, 0},
				{1, -1, 0},
				{1, 1, 0},
				{-1, 1, 0},
				{0, 0, -1},
				{0, 0, 1},
				{0, -1, 0},
				{0, 1, 0},
				{-1, 0, 0},
				{1, 0, 0},
				{0, 0, 0},
			},
			connectivity = {
				{0, 8, 20, 11, 16, 22, 26, 24},
				{8, 1, 9, 20, 22, 17, 25, 26},
				{11, 20, 10, 3, 24, 26, 23, 19},
				{20, 9, 2, 10, 26, 25, 18, 23},
				{16, 22, 26, 24, 4, 12, 21, 15},
				{22, 17, 25, 26, 12, 5, 13, 21},
				{24, 26, 23, 19, 15, 21, 14, 7},
				{26, 25, 18, 23, 21, 13, 6, 14},
			},
		},
	},
	quadrature = [Quadrature_Rule]Reference_Quadrature {
		.Quad_1 = {points = {{0, 0, 0}}, weights = {8}},
		.Quad_3 = {
			points = {
				{-REC_ROOT_3, -REC_ROOT_3, -REC_ROOT_3},
				{REC_ROOT_3, -REC_ROOT_3, -REC_ROOT_3},
				{REC_ROOT_3, REC_ROOT_3, -REC_ROOT_3},
				{-REC_ROOT_3, REC_ROOT_3, -REC_ROOT_3},
				{-REC_ROOT_3, -REC_ROOT_3, REC_ROOT_3},
				{REC_ROOT_3, -REC_ROOT_3, REC_ROOT_3},
				{REC_ROOT_3, REC_ROOT_3, REC_ROOT_3},
				{-REC_ROOT_3, REC_ROOT_3, REC_ROOT_3},
			},
			weights = {1, 1, 1, 1, 1, 1, 1, 1},
		},
		.Quad_5 = {
			points = {
				{-ROOT_3_5, -ROOT_3_5, -ROOT_3_5},
				{0., -ROOT_3_5, -ROOT_3_5},
				{ROOT_3_5, -ROOT_3_5, -ROOT_3_5},
				{-ROOT_3_5, 0., -ROOT_3_5},
				{0., 0., -ROOT_3_5},
				{ROOT_3_5, 0., -ROOT_3_5},
				{-ROOT_3_5, ROOT_3_5, -ROOT_3_5},
				{0., ROOT_3_5, -ROOT_3_5},
				{ROOT_3_5, ROOT_3_5, -ROOT_3_5},
				{-ROOT_3_5, -ROOT_3_5, 0.},
				{0., -ROOT_3_5, 0.},
				{ROOT_3_5, -ROOT_3_5, 0.},
				{-ROOT_3_5, 0., 0.},
				{0., 0., 0.},
				{ROOT_3_5, 0., 0.},
				{-ROOT_3_5, ROOT_3_5, 0.},
				{0., ROOT_3_5, 0.},
				{ROOT_3_5, ROOT_3_5, 0.},
				{-ROOT_3_5, -ROOT_3_5, ROOT_3_5},
				{0., -ROOT_3_5, ROOT_3_5},
				{ROOT_3_5, -ROOT_3_5, ROOT_3_5},
				{-ROOT_3_5, 0., ROOT_3_5},
				{0., 0., ROOT_3_5},
				{ROOT_3_5, 0., ROOT_3_5},
				{-ROOT_3_5, ROOT_3_5, ROOT_3_5},
				{0., ROOT_3_5, ROOT_3_5},
				{ROOT_3_5, ROOT_3_5, ROOT_3_5},
			},
			weights = {
				0.17146776,
				0.27469136,
				0.17146776,
				0.27469136,
				0.43209877,
				0.27469136,
				0.17146776,
				0.27469136,
				0.17146776,
				0.27469136,
				0.43209877,
				0.27469136,
				0.43209877,
				0.68659221,
				0.43209877,
				0.27469136,
				0.43209877,
				0.27469136,
				0.17146776,
				0.27469136,
				0.17146776,
				0.27469136,
				0.43209877,
				0.27469136,
				0.17146776,
				0.27469136,
				0.17146776,
			},
		},
	},
	lagrange = [Basis_Order]Reference_Basis_Lagrange {
		.Linear = {
			support = {
				{.D0, 0, 0, 1},
				{.D0, 1, 0, 1},
				{.D0, 2, 0, 1},
				{.D0, 3, 0, 1},
				{.D0, 4, 0, 1},
				{.D0, 5, 0, 1},
				{.D0, 6, 0, 1},
				{.D0, 7, 0, 1},
			},
			arity = 8,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
				sx := [8]f64{-1, 1, 1, -1, -1, 1, 1, -1}
				sy := [8]f64{-1, -1, 1, 1, -1, -1, 1, 1}
				sz := [8]f64{-1, -1, -1, -1, 1, 1, 1, 1}

				xi := sx[idx]
				eta := sy[idx]
				zeta := sz[idx]

				fx := 1 + xi * r.x
				fy := 1 + eta * r.y
				fz := 1 + zeta * r.z

				val = fx * fy * fz * 0.125
				grad = Vec3{xi * fy * fz * 0.125, eta * fx * fz * 0.125, zeta * fx * fy * 0.125}
				return
			},
		},
		.Quadratic = {
			support = {
				{.D0, 0, 0, 1},
				{.D0, 1, 0, 1},
				{.D0, 2, 0, 1},
				{.D0, 3, 0, 1},
				{.D0, 4, 0, 1},
				{.D0, 5, 0, 1},
				{.D0, 6, 0, 1},
				{.D0, 7, 0, 1},
				{.D1, 0, 0, 1},
				{.D1, 1, 0, 1},
				{.D1, 2, 0, 1},
				{.D1, 3, 0, 1},
				{.D1, 4, 0, 1},
				{.D1, 5, 0, 1},
				{.D1, 6, 0, 1},
				{.D1, 7, 0, 1},
				{.D1, 8, 0, 1},
				{.D1, 9, 0, 1},
				{.D1, 10, 0, 1},
				{.D1, 11, 0, 1},
				{.D2, 0, 0, 1},
				{.D2, 1, 0, 1},
				{.D2, 2, 0, 1},
				{.D2, 3, 0, 1},
				{.D2, 4, 0, 1},
				{.D2, 5, 0, 1},
				{.D3, 0, 0, 1},
			},
			arity = 27,
			reference_values = proc(idx: int, r: Vec3) -> (val: f64, grad: Vec3) {
				node_to_tensor := [27]int {
					0,
					1,
					4,
					3,
					9,
					10,
					13,
					12,
					2,
					7,
					5,
					6,
					11,
					16,
					14,
					15,
					18,
					19,
					22,
					21,
					8,
					17,
					20,
					23,
					24,
					25,
					26,
				}
				t := node_to_tensor[idx]
				ix := t % 3
				iy := (t / 3) % 3
				iz := (t / 9) % 3

				Lx := [3]f64{0.5 * r.x * (r.x - 1), 0.5 * r.x * (r.x + 1), 1 - r.x * r.x}
				dLx := [3]f64{r.x - 0.5, r.x + 0.5, -2 * r.x}
				Ly := [3]f64{0.5 * r.y * (r.y - 1), 0.5 * r.y * (r.y + 1), 1 - r.y * r.y}
				dLy := [3]f64{r.y - 0.5, r.y + 0.5, -2 * r.y}
				Lz := [3]f64{0.5 * r.z * (r.z - 1), 0.5 * r.z * (r.z + 1), 1 - r.z * r.z}
				dLz := [3]f64{r.z - 0.5, r.z + 0.5, -2 * r.z}

				val = Lx[ix] * Ly[iy] * Lz[iz]
				grad = Vec3{dLx[ix] * Ly[iy] * Lz[iz], Lx[ix] * dLy[iy] * Lz[iz], Lx[ix] * Ly[iy] * dLz[iz]}
				return
			},
		},
	},
}
