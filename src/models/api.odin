// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
Interface between core fem library and physical systems that we choose to model.

Will change often as new PDEs break our abstractions, or we look to add broader dynamic behaviour.
*/
package models

import "../fem"
import "../la"

Lagrange_Scalar_Field :: struct {
	using layout:        fem.DOF_Layout,
	order:               fem.Order,
}

Lagrange_Vector_Field :: distinct Lagrange_Scalar_Field

Constraint :: struct {
	data:      rawptr,
	procedure: proc(data: rawptr) -> f64,
}

Variational_BC :: struct {
	data:      rawptr,
	procedure: proc(data: rawptr) -> f64,
}

Source :: struct {
  data: rawptr,
  procuedre: proc(data: rawptr) -> f64
}

lagrange_scalar_field :: proc(
	mesh: fem.Mesh,
	order: fem.Order,
	constrained_facets: map[fem.Boundary_ID]struct{},
	allocator := context.allocator,
) -> (
	field: Lagrange_Scalar_Field,
) {
	field = Lagrange_Scalar_Field {
		order  = .Linear,
		layout = fem.create_layout(mesh, .CG, fem.Lagrange_Scalar, order, constrained_facets, allocator),
	}
	return field
}

lagrange_vector_field :: proc(
	mesh: fem.Mesh,
	order: fem.Order,
	constrained_facets: map[fem.Boundary_ID]struct{},
	allocator := context.allocator,
) -> (
	field: Lagrange_Vector_Field,
) {
	field = Lagrange_Vector_Field {
		order  = .Linear,
		layout = fem.create_layout(mesh, .CG, fem.Lagrange_Vector, order, constrained_facets, allocator),
	}
	return field
}


lagrange_field_scalar_basis :: proc(field: Lagrange_Scalar_Field, element: fem.Element) -> fem.Lagrange_Scalar {
	return fem.basis_create(fem.Lagrange_Scalar, element, field.order)
}

lagrange_field_vector_basis :: proc(field: Lagrange_Vector_Field, element: fem.Element) -> fem.Lagrange_Vector {
	return fem.basis_create(fem.Lagrange_Vector, element, field.order)
}

field_basis :: proc {
	lagrange_field_scalar_basis,
	lagrange_field_vector_basis,
}

lagrange_scalar_field_visualizer :: proc(field: ^Lagrange_Scalar_Field, name: string) -> fem.Visualization_Field {
	procedure := proc(ctx: ^fem.Sample_Context, data: rawptr) -> (r: [fem.MAX_OUTPUT_FIELD_COMPONENTS]f64) {
		field := cast(^Lagrange_Scalar_Field)data
		basis := field_basis(field^, ctx.element)
		for b in 0 ..< basis.arity {
			coeff_index := field.mapping[ctx.element_id][b]
			r[0] += fem.basis_value(basis, ctx, b) * field.coeffs[coeff_index]
		}
		return r
	}

	return {friendly_name = name, components = 1, data = field, value_provider = procedure}
}

lagrange_vector_field_visualizer :: proc(field: ^Lagrange_Vector_Field, name: string) -> fem.Visualization_Field {
	procedure := proc(ctx: ^fem.Sample_Context, data: rawptr) -> (r: [fem.MAX_OUTPUT_FIELD_COMPONENTS]f64) {
		field := cast(^Lagrange_Vector_Field)data
		basis := field_basis(field^, ctx.element)
		for b in 0 ..< basis.arity {
			coeff_index := field.mapping[ctx.element_id][b]
			rv := fem.basis_value(basis, ctx, b) * field.coeffs[coeff_index]
			r[0] += rv[0]
			r[1] += rv[1]
			r[2] += rv[2]
		}
		return r
	}

	return {friendly_name = name, components = 1, data = field, value_provider = procedure}
}


field_visualizer :: proc {
	lagrange_scalar_field_visualizer,
	lagrange_vector_field_visualizer,
}

field_set_constraints :: proc(mesh: fem.Mesh, field: $T, constraints: map[fem.Boundary_ID]Constraint) {
	for element, id in mesh.elements {
		basis := field_basis(field, element)
		for adjacency, facet_index in element.adjacency {
			if bnd_id, ok := adjacency.(fem.Boundary_ID); ok {
				constraint := constraints[bnd_id] or_continue
				for basis_index in 0 ..< basis.arity {
					global_index := field.mapping[id][basis_index]
					if !field.constrained[global_index] {continue}
					field.coeffs[global_index] = constraint.procedure(constraint.data)
				}
			}
		}
	}
}
