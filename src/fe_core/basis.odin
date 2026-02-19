// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
/*
Basis interface.

Basis_Info is information used by code doesn't evaluate the basis, it can be constructed through runtime information if needed.

Basis itself is a compile time interface through static dispatch for performance.
*/
package fem

Basis_Support :: struct {
	entity_dim:         Dimension, //which entity supports the basis
	entity_index:       int, // which local entity index on the element is this (vertex 0, face 1, etc.)
	entity_basis_index: int, // which basis on the specific entity is this.
	entity_basis_arity: int, // how many total dofs exist on this entity.
}

Basis_Info :: struct {
	order:             Order,
	geometry_required: Geometry_Options,
	arity:             int,
	// should not be accessed directly, use basis_support
	_support:          []Basis_Support,
	// how many dofs per basis function, generally 1 for scalar basis or true vector valued basis (RT, Nedelec etc).
	// For vector basis built from multiple scalar basis components may be > 1.
	components:        int,
}

Basis_Family :: enum {
	LS,
	LV,
}

basis_create :: proc($T: typeid, element: Element, order: Order) -> T {
	when T == Basis_LS {
		return ls_create(element, order)
	} else when T == Basis_LV {
		return lv_create(element, order)} else {#panic("Basis create recieved a type that is not a basis.")
	}
}

basis_support :: proc(bi: Basis_Info, basis: int) -> Basis_Support {
	assert(basis >= 0 && basis < bi.arity)
	index, cmpnt := basis_decompose_component(bi, basis)
	scalar_support := bi._support[index]
	scalar_support.entity_basis_index = (scalar_support.entity_basis_index * bi.components) + cmpnt
	scalar_support.entity_basis_arity *= bi.components
	return scalar_support
}

@(private)
basis_decompose_component :: #force_inline proc(bi: Basis_Info, basis: int) -> (scalar_index, cmpnt: int) {
	return basis / bi.components, basis % bi.components
}


basis_value :: proc {
	ls_value,
	lv_value,
}

basis_gradient :: proc {
	ls_gradient,
	lv_gradient,
}
