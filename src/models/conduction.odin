// SPDX-FileCopyrightText: 2026 Rowan Apps, Tor Rabien
// SPDX-License-Identifier: MIT
package models

import "core:log"
import "core:math/linalg"

import "../fem"
import "../la"

conduction :: proc(
	mesh: fem.Mesh,
	temperature: Lagrange_Scalar_Field,
	conductivity: f64,
	weak_boundaries: map[fem.Boundary_ID]Variational_BC,
	allocator := context.allocator,
) -> bool {
	context.allocator = allocator

	stiffness_sparsity := fem.create_sparsity_from_layout(mesh, temperature, temperature)
	stiffness := la.sparse_mat_from_sparsity(stiffness_sparsity)
	load := make(la.Vec, len(temperature.coeffs))

	defer {
		la.sparse_mat_destroy(stiffness)
		delete(load)
		delete(stiffness_sparsity.row_ptrs)
		delete(stiffness_sparsity.columns)
	}


	for element, id in mesh.elements {
		basis := field_basis(temperature, element)
		si, ctx := fem.sample_context_create(mesh, fem.Entity_ID(id))

		mask := fem.create_interface_facet_filter(element)
		for fem.iterate_facet_sample_group(&si, &ctx, .Facet_Integration_3, mask) {
			bound_id := element.adjacency[ctx.fs.facet_index].(fem.Boundary_ID)
			if bound_id not_in weak_boundaries {continue}
			for test in 0 ..< basis.arity {
				global_test := temperature.mapping[id][test]
				if temperature.current_constraints[global_test] {continue} 	// possible on boundary between dirichlet and neumman
				term := weak_boundaries[bound_id]
				res := term.procedure(term.data) * fem.dS(&ctx) * fem.basis_value(basis, &ctx, test)
				load[global_test] += res

			}
		}

		fem.sample_iterator_reset(&si)

		for fem.iterate_sample_group(&si, &ctx, .Integration_3) {
			for test in 0 ..< basis.arity {
				global_test := temperature.mapping[id][test]
				if temperature.current_constraints[global_test] {continue}
				test_grad := fem.basis_gradient(basis, &ctx, test)
				for trial in 0 ..< basis.arity {
					global_trial := temperature.mapping[id][trial]
					trial_grad := fem.basis_gradient(basis, &ctx, trial)
					integrand := linalg.dot(test_grad, trial_grad) * conductivity * fem.dV(&ctx)
					if temperature.current_constraints[global_trial] {
						load[global_test] -= integrand * temperature.coeffs[global_trial]
						continue
					}
					idx := la.sparse_mat_idx(stiffness, global_test, global_trial)
					stiffness.values[idx] += integrand
				}
			}
		}
	}

	iter, resid, converged := la.sparse_cg_solve(stiffness, temperature.coeffs, load)

	log.infof("Sparse solve: converged? %v, residual: %v, iterations: %v", converged, resid, iter)

	return converged
}


isothermal :: proc(T: f64) -> Constraint {
	return {transmute(rawptr)T, proc(data: rawptr) -> f64 {return transmute(f64)data}}
}

fixed_flux :: proc(q: f64) -> Variational_BC {
	return {transmute(rawptr)q, proc(data: rawptr) -> f64 {return transmute(f64)data}}
}
