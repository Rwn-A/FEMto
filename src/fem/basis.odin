/*
 Basis interfaces.

 Broken into two:
    - Basis interface: contains basic specific information
      used outside of hot loops for DOF maps and similar.

    - Local space interface: Just holds values and knows how to map to
      physical space.

 Gradient space is for basis families that behave like an H1 space locally.
 (gradient well defined, uses the jacobian inverse to transform gradients)
*/
package fem

VECTOR_CMPNTS :: 3

Basis_Family :: enum {
	Lagrange,
}

Basis_Descriptor :: struct {
	family: Basis_Family,
	order:  Basis_Order,
}

Grad_Space_Table :: struct {
	values:          []f64,
	gradients:       []Vec3,
	n_points, n_dof: int,
}

Grad_Space_Shape :: enum {
	Scalar,
	Vector,
}

Grad_Space :: struct(shape: Grad_Space_Shape) {
	using tbl: Grad_Space_Table,
}

basis_support :: proc(element_type: Element_Type, bd: Basis_Descriptor) -> []Basis_Support {
	switch bd.family {
	case .Lagrange:
		return REFERENCE_ELEMENTS[element_type].lagrange[bd.order].support
	case:
		panic("Unreachable")
	}
}

basis_count :: proc(element_type: Element_Type, bd: Basis_Descriptor) -> int {
	switch bd.family {
	case .Lagrange:
		return REFERENCE_ELEMENTS[element_type].lagrange[bd.order].arity
	case:
		panic("Unreachable")
	}
}

// used for internally for precomputing.
basis_tabulate :: proc(
	element_type: Element_Type,
	bd: Basis_Descriptor,
	points: []Vec3,
	allocator := context.allocator,
) -> Basis_Table {
	context.allocator = allocator

	switch bd.family {
	case .Lagrange:
		n_dof := basis_count(element_type, bd)
		n_points := len(points)
		tbl := Grad_Space_Table {
			values    = make([]f64, n_dof * n_points),
			gradients = make([]Vec3, n_dof * n_points),
			n_points  = n_points,
			n_dof     = n_dof,
		}

		for point, point_idx in points {
			for dof, dof_idx in 0 ..< n_dof {
				idx := grad_space_idx(tbl, point_idx, dof_idx)
				val, grad := REFERENCE_ELEMENTS[element_type].lagrange[bd.order].reference_values(dof_idx, point)
				tbl.values[idx] = val
				tbl.gradients[idx] = grad
			}
		}

		return tbl
	case:
		panic("unreachable")
	}

	return {}
}


// Returns the point rule your field must be evaluated over to project onto dofs.
basis_dof_functional_rule :: proc(
	element: Element_Type,
	bd: Basis_Descriptor,
	basis_dof: int, //TODO: rename basis_dof
) -> Point_Rule {
	switch bd.family {
	case .Lagrange:
		// this only works for lagrange and is a bit hacky but it avoids needing a bunch of 1-node rules.
		// basically take the subcell rule and strip down only the information for one of the points that corresponds to the
		// chosen dof.
		rule := SUBCELL_POINT_RULES[element][bd.order]

		rule.points = rule.points[basis_dof:basis_dof + 1]

		// only lets you eval the provided basis on the rule, seems reasonable but idk
		tbl := rule.basis_tables[bd.family][bd.order].(Grad_Space_Table)
		// dof in this case is really our point, so all basis values just for that point, whcih for lagrange is just 0s excpet for a single 1 but alas.
		vals := tbl.values[basis_dof * basis_count(element, bd):(basis_dof + 1) * basis_count(element, bd)]

		grads := tbl.gradients[basis_dof * basis_count(element, bd):(basis_dof + 1) * basis_count(element, bd)]

		rule.basis_tables[bd.family][bd.order] = Grad_Space_Table{vals, grads, 1, basis_count(element, bd)}

		return rule
	case:
		unreachable()
	}
}

scalar_basis_dof_functional :: proc(
	element: Element_Type,
	bd: Basis_Descriptor,
	dof: int,
	field_values: []f64,
) -> f64 {
	rule := basis_dof_functional_rule(element, bd, dof)

	assert(len(field_values) == len(rule.points))

	switch bd.family {
	case .Lagrange:
		return field_values[0]
	case:
		unreachable()
	}
}
vector_basis_dof_functional :: proc(
	element: Element_Type,
	bd: Basis_Descriptor,
	basis_dof: int,
	component: int, // component is for repeated vector basis, not native vector basis like RT.
	field_values: []Vec3,
) -> f64 {
	rule := basis_dof_functional_rule(element, bd, basis_dof)

	assert(len(field_values) == len(rule.points))

	switch bd.family {
	case .Lagrange:
		return field_values[0][component]
	case:
		unreachable()
	}
}


basis_dof_functional :: proc {
	scalar_basis_dof_functional,
	vector_basis_dof_functional,
}

basis_grad_space :: proc(phys: Mapped_Element, bd: Basis_Descriptor, $shape: Grad_Space_Shape) -> Grad_Space(shape) {
	tbl, exists := phys.rule.basis_tables[bd.family][bd.order].(Grad_Space_Table)
	assert(exists, "basis family at requested order was not tabulated as a gradient space on given point rule.")
	return {tbl}
}

grad_space_points :: proc(s: Grad_Space_Table) -> int {
	return s.n_points
}

grad_space_scalar_value :: proc(s: Grad_Space(.Scalar), qp, dof: int) -> f64 {
	return s.values[grad_space_idx(s, qp, dof)]
}

grad_space_scalar_gradient :: proc(s: Grad_Space(.Scalar), phys: Mapped_Element, qp, dof: int) -> Vec3 {
	return phys.im.pullback[qp] * s.gradients[grad_space_idx(s, qp, dof)]
}

grad_space_scalar_arity :: proc(s: Grad_Space(.Scalar)) -> int {
	return s.n_dof
}

grad_space_vector_arity :: proc(s: Grad_Space(.Vector)) -> int {
	return s.n_dof * VECTOR_CMPNTS
}

grad_space_vector_value :: proc(s: Grad_Space(.Vector), qp, dof: int) -> (v: Vec3) {
	sdof, cmpnt := grad_space_decompose(s, dof)
	v[cmpnt] = grad_space_scalar_value(cast(Grad_Space(.Scalar))s, qp, sdof)
	return v
}

grad_space_vector_gradient :: proc(s: Grad_Space(.Vector), phys: Mapped_Element, qp, dof: int) -> (m: Mat3) {
	sdof, cmpnt := grad_space_decompose(s, dof)

	scalar_gradient := grad_space_scalar_gradient(cast(Grad_Space(.Scalar))s, phys, qp, sdof)

	m[cmpnt, 0] = scalar_gradient[0]
	m[cmpnt, 1] = scalar_gradient[1]
	m[cmpnt, 2] = scalar_gradient[2]

	return m
}

grad_space_vector_symmetric_gradient :: proc(
	s: Grad_Space(.Vector),
	phys: Mapped_Element,
	qp, dof: int,
) -> (
	g: Voigt6,
) {
	m := space_gradient(s, phys, qp, dof)

	g[0] = m[0, 0]
	g[1] = m[1, 1]
	g[2] = m[2, 2]

	g[3] = (m[1, 2] + m[2, 1]) * 0.5
	g[4] = (m[0, 2] + m[2, 0]) * 0.5
	g[5] = (m[0, 1] + m[1, 0]) * 0.5

	return g
}

space_points :: proc {
	grad_space_points,
}

space_arity :: proc {
	grad_space_scalar_arity,
	grad_space_vector_arity,
}

space_value :: proc {
	grad_space_scalar_value,
	grad_space_vector_value,
}

space_gradient :: proc {
	grad_space_scalar_gradient,
	grad_space_vector_gradient,
}

space_symmetric_gradient :: proc {
	grad_space_vector_symmetric_gradient,
}


// does not respect space shape, used internally, would be confusing if exposed.
@(private)
grad_space_idx :: proc(tbl: Grad_Space_Table, point, dof: int) -> int {
	assert(point < tbl.n_points)
	return point * tbl.n_dof + dof
}

@(private)
grad_space_decompose :: proc(s: Grad_Space($E), dof: int) -> (scalar_dof, component: int) {
	#assert(E != .Scalar, "Tried to decompose scalar space into components, likely unintentional.")
	#partial switch E {
	case .Vector:
		return dof / VECTOR_CMPNTS, dof % VECTOR_CMPNTS
	case:
		panic("Unreachable")
	}
}
