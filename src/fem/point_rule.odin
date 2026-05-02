// A point rule is a set of points with precomputed basis values.
// An example of a point rule is a quadrature rule.
// The points are reference space values, geometry transforms are handled via the mapped element.
package fem

Basis_Table :: union {
	Grad_Space_Table,
}

Point_Rule_Flag :: enum {
	Facet_Supported,
	Weighted,
}

Point_Rule :: struct {
	element_type: Element_Type,
	points:       []Vec3,
	flags:        bit_set[Point_Rule_Flag],
	// possibly sparse
	basis_tables: [Basis_Family][Basis_Order]Basis_Table,
	//optional: depending on flags
	weights:      []f64,
	facet_index:  int,
}

INTERIOR_QUADRATURE_POINT_RULES := [Element_Type][Quadrature_Rule]Point_Rule{}

SURFACE_QUADRATURE_POINT_RULES := [Element_Type][Quadrature_Rule][]Point_Rule{}

SUBCELL_POINT_RULES := [Element_Type][Basis_Order]Point_Rule{}


// Use an allocator with free all support, no seperate destroy proc.
setup_default_rules :: proc(allocator := context.allocator) {
	context.allocator = allocator

	for element in Element_Type {
		for rule in Quadrature_Rule {
			pr: Point_Rule
			pr.element_type = element
			pr.flags = {.Weighted}
			pr.weights = REFERENCE_ELEMENTS[element].quadrature[rule].weights
			pr.points = REFERENCE_ELEMENTS[element].quadrature[rule].points

			tabulate_all_basis(&pr)

			INTERIOR_QUADRATURE_POINT_RULES[element][rule] = pr
		}

		for rule in Basis_Order {
			pr: Point_Rule
			pr.element_type = element
			pr.points = REFERENCE_ELEMENTS[element].subcell[rule].points

			tabulate_all_basis(&pr)

			SUBCELL_POINT_RULES[element][rule] = pr
		}
	}

	for element in Element_Type {
		if element == .Point {continue}
		for rule in Quadrature_Rule {
			SURFACE_QUADRATURE_POINT_RULES[element][rule] = make([]Point_Rule, element_facet_count(element))
			for facet_index in 0 ..< element_facet_count(element) {
				facet_type := element_facet_type(element, facet_index)
				facet_rule := INTERIOR_QUADRATURE_POINT_RULES[facet_type][rule]

				pr: Point_Rule
				pr.element_type = element
				pr.flags = {.Weighted, .Facet_Supported}
				pr.weights = facet_rule.weights

				pr.points = make([]Vec3, len(facet_rule.points))
				for point, i in facet_rule.points {
					pr.points[i] = facet_reference_to_volume_reference(element, point, facet_index)
				}

				tabulate_all_basis(&pr)

				SURFACE_QUADRATURE_POINT_RULES[element][rule][facet_index] = pr

			}
		}
	}

	tabulate_all_basis :: proc(pr: ^Point_Rule) {
		for family in Basis_Family {
			for order in Basis_Order {
				pr.basis_tables[family][order] = basis_tabulate(pr.element_type, {family, order}, pr.points)
			}
		}
	}
}
