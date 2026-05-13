package fem

import "core:log"

import "core:mem"
import "core:slice"

Nonlinear_Iterator :: struct {
    max_iters:              int,
    current_iter:           int,
    tolerance:              f64,
    initial_residual:       f64,
    trivial_soln_tolerance: f64,
}

nli_create :: proc(tolerance: f64, max_iterations: int, trivial_soln_tol := 1e-14) -> Nonlinear_Iterator {
    return Nonlinear_Iterator{
        tolerance              = tolerance,
        max_iters              = max_iterations,
        trivial_soln_tolerance = trivial_soln_tol,
    }
}

nli_step :: proc(ni: ^Nonlinear_Iterator, solution: Vector, tangent: Sparse_Matrix, residual, update: Vector) -> (should_solve: bool, ok: bool) {
    defer ni.current_iter += 1

    switch ni.current_iter {
    case 0:
        slice.zero(tangent.values)
        slice.zero(residual)
        return false, true

    case 1:
        ni.initial_residual = nrm2(residual)
        if ni.initial_residual < ni.trivial_soln_tolerance { return false, false }

    case:
        if nrm2(residual) < ni.tolerance * ni.initial_residual { return false, false }
    }

    if ni.current_iter >= ni.max_iters { return false, false }

    return true, true
}

nli_do_update :: proc(solution: Vector, tangent: Sparse_Matrix, residual, update: Vector) {
    axpy(update, solution, 1.0)
    slice.zero(residual)
    slice.zero(tangent.values)
}

nli_reached_max :: proc(ni: Nonlinear_Iterator) -> bool {
    return ni.current_iter >= ni.max_iters
}

nli_linear_tolerance :: proc(ni: Nonlinear_Iterator) -> f64 {
    return ni.tolerance * 0.1
}