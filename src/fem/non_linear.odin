package fem

import "core:log"

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

nli_create_state :: proc(sys: System, ics: Vector, allocator := context.allocator) -> Nonlinear_Iter_State {
	return {
		solution = slice.clone(ics),
		tangent = system_matrix(sys, allocator),
		residual = system_vector(sys, allocator),
		update = system_vector(sys, allocator),
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
	should_solve := ni.current_iter > 0

	if ni.current_iter == 0 {
		slice.zero(state.residual)
		slice.zero(state.tangent.values)
	}

	if ni.current_iter == 1 {
		ni.initial_residual = nrm2(state.residual)

		if ni.initial_residual < ni.trivial_soln_tolerance {
			return {}, false, false
		}
	}

	if ni.current_iter > 1 {
		if nrm2(state.residual) < ni.tolerance * ni.initial_residual {
			return {}, false, false
		}
	}

	if ni.current_iter >= ni.max_iters {return {}, false, false}
	defer ni.current_iter += 1

	return ni.current_iter, should_solve, true
}

nli_reached_max :: proc(ni: ^Nonlinear_Iterator) -> bool {
	return ni.current_iter >= ni.max_iters
}

nli_update :: proc(ni: ^Nonlinear_Iterator, state: Nonlinear_Iter_State) {
	axpy(state.update, state.solution, 1.0)
	slice.zero(state.residual)
	slice.zero(state.tangent.values)
}

// optional reasonable tolerance to use for a linear solve
nli_linear_tolerance :: proc(ni: Nonlinear_Iterator) -> f64 {
	// linear tolerance slightly tighter then non linear, can be changed
	// with fancier methods down the road.
	return ni.tolerance * 0.1
}
