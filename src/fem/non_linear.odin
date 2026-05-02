package fem

import "core:mem"
import "core:slice"

Nonlinear_Iterator :: struct {
	max_iters, current_iter: int,
	tolerance:               f64,
	initial_residual:        f64,
	trivial_soln_tolerance:  f64,
}

Nonlinear_Iter_State :: struct {
	residual:         Vector,
	tangent:          Sparse_Matrix,
	solution, update: Vector,
}

nli_create_state :: proc(sys: System, allocator := context.allocator) -> Nonlinear_Iter_State {
	return {
		system_vector(sys, allocator),
		system_matrix(sys, allocator),
		system_vector(sys, allocator),
		system_vector(sys, allocator),
	}
}

nli_create :: proc(tolerance: f64, max_iterations: int, trivial_soln_tol := 1e-14) -> Nonlinear_Iterator {
	return Nonlinear_Iterator {
		tolerance = tolerance,
		max_iters = max_iterations,
		trivial_soln_tolerance = trivial_soln_tol,
	}
}


nli_step :: proc(ni: ^Nonlinear_Iterator, state: Nonlinear_Iter_State) -> (int, bool, bool) {
	if ni.current_iter == 0 {

	} else if ni.current_iter == 1 {
		ni.initial_residual = nrm2(state.residual)
		if ni.initial_residual < ni.trivial_soln_tolerance {
			return {}, false, false
		}
	} else {
		if nrm2(state.residual) < ni.tolerance * ni.initial_residual {
			return {}, false, false
		}
		axpy(state.update, state.solution, 1.0)
		slice.zero(state.residual)
		slice.zero(state.tangent.values)
	}
	should_solve := ni.current_iter > 0
	if ni.current_iter >= ni.max_iters {return {}, false, false}
	defer ni.current_iter += 1
	return ni.current_iter, should_solve, true
}

// optional reasonable tolerance to use for a linear solve
nli_linear_tolerance :: proc(ni: Nonlinear_Iterator) -> f64 {
	// linear tolerance slightly tighter then non linear, can be changed
	// with fancier methods down the road.
	return ni.tolerance * 0.1
}
// Nonlinear_Iterator :: struct {
// 	u_cur:        Vector,
// 	R:            Vector,
// 	J:            Sparse_Matrix,
// 	du:           Vector,
// 	tolerance:    f64,
// 	max_iters:    int,
// 	current_iter: int,
// 	converged:    bool,
// 	allocator:    mem.Allocator,
// }


// // clones ics
// nli_create :: proc(
// 	sys: System,
// 	ics: Vector,
// 	tolerance: f64,
// 	max_iters: int,
// 	allocator := context.allocator,
// ) -> Nonlinear_Iterator {
// 	return {
// 		u_cur = slice.clone(ics, allocator),
// 		R = system_vector(sys, allocator),
// 		J = system_matrix(sys, allocator),
// 		du = system_vector(sys, allocator),
// 		tolerance = tolerance,
// 		max_iters = max_iters,
// 		allocator = allocator,
// 	}
// }

// nli_destroy :: proc(ni: ^Nonlinear_Iterator) {
// 	delete(ni.u_cur, ni.allocator)
// 	delete(ni.R, ni.allocator)
// 	delete(ni.du, ni.allocator)
// 	sparse_matrix_destroy(ni.J, ni.allocator)
// }

// nli_current :: proc(ni: Nonlinear_Iterator) -> Vector {
// 	return ni.u_cur
// }

// nli_residual :: proc(ni: Nonlinear_Iterator) -> Vector {
// 	return ni.R
// }

// nli_du :: proc(ni: Nonlinear_Iterator) -> Vector {
// 	return ni.du
// }

// nli_jacobian :: proc(ni: Nonlinear_Iterator) -> Sparse_Matrix {
// 	return ni.J
// }


// nli_step :: proc(ni: ^Nonlinear_Iterator) -> (int, bool) {
// 	if ni.current_iter >= ni.max_iters || ni.converged {return {}, false}
// 	defer ni.current_iter += 1
// 	return ni.current_iter, true
// }

// // convergence check, used to exit the non linear loop as nli_should_continue(ni) or_break
// nli_should_continue :: proc(ni: ^Nonlinear_Iterator) -> bool {
// 	ni.converged = nrm2(ni.R) < ni.tolerance
// 	return !ni.converged
// }

// // updates the current solution with the current du, to be called after a linear solve, before the next iteration.
// nli_update :: proc(ni: ^Nonlinear_Iterator) {
// 	if ni.converged {return} 	// dont update if we already convegred
// 	axpy(ni.du, ni.u_cur, 1.0)
// 	slice.zero(ni.R)
// 	slice.zero(ni.J.values)
// }
