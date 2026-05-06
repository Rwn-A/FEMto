// A mapped element is an element with geometry mapped to physical space.
// Geometry is mapped over a point rule.
// Multiple point rules (ex. volume and surface quadrature) require multiple mappings.
// Mapped element assumes 3D ambient dimension, lower dimensional meshes are considered to be
// embedded in 3D.
package fem

import "core:math/linalg"
import "core:slice"

Element :: struct {
	type:   Element_Type,
	order:  Element_Order,
	affine: bool,
	nodes:  []Vec3,
}

Mapped_Element :: struct {
	rule:    ^Point_Rule,
	element: Mesh_Element,
	im:      #soa[]Interior_Geometry,
	fm:      #soa[]Facet_Geometry,
}


Interior_Geometry :: struct {
	jacobian, pullback: Mat3,
	measure:            f64,
	point:              Vec3,
}

Facet_Geometry :: struct {
	jacobian, pullback: Mat3,
	measure:            f64,
	normal:             Vec3,
}

Quadrature_Region :: enum {
	Interior,
	Surface,
}

Quadrature_Descriptor :: struct {
	region: Quadrature_Region,
	rule:   Quadrature_Rule,
	facet:  int,
}

map_quadrature :: proc(element: Mesh_Element, qd: Quadrature_Descriptor, allocator := context.allocator) -> Mapped_Element {
	switch qd.region {
	case .Interior:
		return map_element(element, &INTERIOR_QUADRATURE_POINT_RULES[element.type][qd.rule], allocator)
	case .Surface:
		return map_element(element, &SURFACE_QUADRATURE_POINT_RULES[element.type][qd.rule][qd.facet], allocator)
	case:
		panic("Unreachable")
	}
}


// remap is slightly faster in the case of an affine element, if the element is non-affine it is
// of comparable speed to just creating a new mapping.
remap_quadrature :: proc(
	original: Mapped_Element,
	qd: Quadrature_Descriptor,
	allocator := context.allocator,
) -> Mapped_Element {
	switch qd.region {
	case .Interior:
		return remap_element(original, &INTERIOR_QUADRATURE_POINT_RULES[original.element.type][qd.rule], allocator)
	case .Surface:
		return remap_element(
			original,
			&SURFACE_QUADRATURE_POINT_RULES[original.element.type][qd.rule][qd.facet],
			allocator,
		)
	case:
		panic("Unreachable")
	}
}

dV :: proc(phys: Mapped_Element, qp: int) -> f64 {
	assert(
		.Weighted in phys.rule.flags,
		"Differential volume is only valid when called on a quadrature mapped element.",
	)
	return phys.im.measure[qp] * phys.rule.weights[qp]
}

dS :: proc(phys: Mapped_Element, qp: int) -> f64 {
	assert(
		.Weighted & .Facet_Supported in phys.rule.flags,
		"Differential surface area is only valid when called on a facet quadrature mapped element.",
	)
	return phys.fm.measure[qp] * phys.rule.weights[qp]
}

infer_quadrature :: proc(highest_basis: Basis_Order) -> Quadrature_Rule {
	needed := (int(highest_basis) + 1) * 2 + 1

	switch needed {
	case 0, 1, 2, 3:
		return .Quad_3
	case:
		return .Quad_5
	}
}

// Use for custom point rules, if doing quadrature or visualization prefer specific functions.
map_element :: proc(element: Mesh_Element, rule: ^Point_Rule, allocator := context.allocator) -> (phys: Mapped_Element) {
	context.allocator = allocator

	num_points := len(rule.points)

	phys.rule = rule
	phys.element = element
	phys.im = make(#soa[]Interior_Geometry, num_points)

	space_table := rule.basis_tables[.Lagrange][element.order].(Grad_Space_Table)

	for idx in 0 ..< num_points {
		geo := Interior_Geometry{}

		geo.jacobian = interior_jacobian(element, space_table, idx)
		geo.pullback = pullback(element_dim(element.type), geo.jacobian)
		geo.measure = measure(element_dim(element.type), geo.jacobian)

		if element.affine {
			fill_affine_interior(phys.im, geo)
			break
		}

		phys.im[idx] = geo
	}

	for idx in 0 ..< num_points {
		phys.im[idx].point = physical_point(element, space_table, idx)
	}


	if .Facet_Supported not_in rule.flags {return phys}

	map_facet(&phys)

	return phys
}

remap_element :: proc(
	original: Mapped_Element,
	rule: ^Point_Rule,
	allocator := context.allocator,
) -> (
	phys: Mapped_Element,
) {
	context.allocator = allocator

	if !original.element.affine {return map_element(original.element, rule)}

	num_points := len(rule.points)

	phys.rule = rule
	phys.element = original.element
	phys.im = make(#soa[]Interior_Geometry, num_points)

	fill_affine_interior(phys.im, original.im[0])

	space_table := rule.basis_tables[.Lagrange][phys.element.order].(Grad_Space_Table)
	for idx in 0 ..< num_points {
		phys.im[idx].point = physical_point(phys.element, space_table, idx)
	}

	if .Facet_Supported not_in rule.flags {return phys}

	map_facet(&phys)

	return phys
}

mapped_destroy :: proc(phys: ^Mapped_Element, allocator := context.allocator) {
	delete(phys.im, allocator)

	if .Facet_Supported in phys.rule.flags {
		delete(phys.fm, allocator)
	}
}

element_reduce_to_linear :: proc(element: Element) -> (reduced: Element) {
	reduced = element
	reduced.order = .Linear
	reduced.nodes = element.nodes[:len(element_sub_entity_nodes_all(element.type, .Linear, .D0))]
	return reduced
}


@(private = "file")
map_facet :: proc(phys: ^Mapped_Element) {
	n_points := len(phys.rule.points)
	facet := phys.rule.facet_index

	phys.fm = make(#soa[]Facet_Geometry, n_points)
	facet_type := element_facet_type(phys.element.type, facet)
	facet_dim := element_dim(facet_type)

	for idx in 0 ..< n_points {
		geo := Facet_Geometry{}

		geo.jacobian = facet_jacobian(phys.element, phys.im.jacobian[idx], facet)
		geo.pullback = pullback(facet_dim, geo.jacobian)
		geo.measure = measure(facet_dim, geo.jacobian)
		geo.normal = facet_normal(facet_dim, geo.jacobian, phys.im.jacobian[idx])

		if phys.element.affine {
			slice.fill(phys.fm.jacobian[:n_points], geo.jacobian)
			slice.fill(phys.fm.measure[:n_points], geo.measure)
			slice.fill(phys.fm.pullback[:n_points], geo.pullback)
			slice.fill(phys.fm.normal[:n_points], geo.normal)
			break
		} else {
			phys.fm[idx] = geo
		}
	}
}

@(private = "file")
fill_affine_interior :: proc(im: #soa[]Interior_Geometry, geo: Interior_Geometry) {
	n := len(im)
	slice.fill(im.jacobian[:n], geo.jacobian)
	slice.fill(im.measure[:n], geo.measure)
	slice.fill(im.pullback[:n], geo.pullback)
}

physical_point :: proc(element: Element, lagrange_ref: Grad_Space_Table, point: int) -> (p: Vec3) {
	for node, basis_idx in element.nodes {
		p += lagrange_ref.values[grad_space_idx(lagrange_ref, point, basis_idx)] * node
	}
	return p
}

interior_jacobian :: proc(element: Element, lagrange_ref: Grad_Space_Table, point: int) -> (j: Mat3) {
	for node, basis_idx in element.nodes {
		j += linalg.outer_product(node, lagrange_ref.gradients[grad_space_idx(lagrange_ref, point, basis_idx)])
	}
	return j
}

facet_jacobian :: proc(element: Element, interior_jacobian: Mat3, facet: int) -> (j: Mat3) {
	if element.type == .Line {
		j[0, 0] = 1
		return
	}
	for tangent, i in element_facet_tangents(element.type, facet) {j[i] = interior_jacobian * tangent}
	return j
}


// The measure is the determinant of the jacobian, or pseudo determinant in the 1D, 2D cases.
measure :: proc(dim: Dimension, jacobian: Mat3) -> f64 {
	switch dim {
	case .D0:
		return 1
	case .D1:
		return linalg.length(jacobian[0])
	case .D2:
		return linalg.length(linalg.cross(jacobian[0], jacobian[1]))
	case .D3:
		return linalg.determinant(jacobian)
	case:
		unreachable()
	}
}

// The pullback is the inverse transpose of the jacobian or the pseudo inverse.
// Because the simulation is always ambiently 3D, the returned matrix has some zero
// entries for 1D and 2D meshes.
pullback :: proc(dim: Dimension, jacobian: Mat3) -> (transform: Mat3) {
	switch dim {
	case .D0:
		return {}
	case .D1:
		len_sq := linalg.length2(jacobian[0])
		transform[0] = jacobian[0] / len_sq
		return linalg.transpose(transform)
	case .D2:
		J: matrix[3, 2]f64
		J[0] = jacobian[0]
		J[1] = jacobian[1]
		G_inv := linalg.inverse(linalg.transpose(J) * J)
		mat2x3 := linalg.transpose(J * G_inv)
		transform = matrix[3, 3]f64{
			mat2x3[0, 0], mat2x3[0, 1], mat2x3[0, 2],
			mat2x3[1, 0], mat2x3[1, 1], mat2x3[1, 2],
			0, 0, 0,
		}
		return linalg.transpose(transform)
	case .D3:
		return linalg.transpose(linalg.inverse(jacobian))
	case:
		unreachable()
	}
}

// Always points outward from the element jacobian pased in `parent_jacobian`.
facet_normal :: proc(facet_dim: Dimension, facet_jacobian: Mat3, parent_jacobian: Mat3) -> Vec3 {
	assert(facet_dim != .D3)
	switch facet_dim {
	case .D0:
		return linalg.normalize(parent_jacobian[0])
	case .D1:
		edge_tangent := linalg.normalize(facet_jacobian[0])
		surface_normal := linalg.normalize(linalg.cross(parent_jacobian[0], parent_jacobian[1]))
		return linalg.normalize(linalg.cross(edge_tangent, surface_normal))
	case .D2:
		normal := linalg.cross(facet_jacobian[0], facet_jacobian[1])
		return linalg.normalize(normal)
	case .D3:
		unreachable()
	case:
		unreachable()
	}
}
