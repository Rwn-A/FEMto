// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package models

import "core:log"
import "core:math/linalg"

import "../fem"
import "../la"

import "core:time"

conduction :: proc(
	mesh: fem.Mesh,
	element_id: fem.Entity_ID,
	temperature: Lagrange_Scalar_Field,
	conductivity: f64,
	weak_boundaries: map[fem.Boundary_ID]Variational_BC,
	allocator := context.allocator,
) -> (la.Dense_Matrix, la.Vec) {
	element := mesh.elements[element_id]

	basis := field_basis(temperature, element)
	si, ctx := fem.sample_context_create(mesh, fem.Entity_ID(element_id))

	k := la.dense_create(basis.arity, basis.arity, allocator)
	f := make(la.Vec, basis.arity, allocator)

	// weak bcs
	mask := fem.create_interface_facet_filter(element)
	for fem.iterate_facet_sample_group(&si, &ctx, .Facet_Integration_3, mask) {
		bound_id := element.adjacency[ctx.fs.facet_index].(fem.Boundary_ID)
		if bound_id not_in weak_boundaries {continue}
		for test in 0 ..< basis.arity {
			term := weak_boundaries[bound_id]
			res := term.procedure(term.data) * fem.dS(&ctx) * fem.basis_value(basis, &ctx, test)
			f[test] += res
		}
	}

	fem.sample_iterator_reset(&si)

	// volume terms
	for fem.iterate_sample_group(&si, &ctx, .Integration_3) {
		for test in 0 ..< basis.arity {
			test_grad := fem.basis_gradient(basis, &ctx, test)
			for trial in 0 ..< basis.arity {
				trial_grad := fem.basis_gradient(basis, &ctx, trial)
				integrand := linalg.dot(test_grad, trial_grad) * conductivity * fem.dV(&ctx)
				idx := la.dense_mat_idx(k, test, trial)
				k.values[idx] += integrand
			}
		}
	}

	return k, f
}


isothermal :: proc(T: f64) -> Constraint {
	return {transmute(rawptr)T, proc(data: rawptr) -> f64 {return transmute(f64)data}}
}

fixed_flux :: proc(q: f64) -> Variational_BC {
	return {transmute(rawptr)q, proc(data: rawptr) -> f64 {return transmute(f64)data}}
}

volumetric_flux :: proc(q: f64) -> Source {
	return {transmute(rawptr)q, proc(data: rawptr) -> f64 {return transmute(f64)data}}
}

