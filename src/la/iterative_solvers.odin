// SPDX-FileCopyrightText: 2026 Rowan Apps
// SPDX-License-Identifier: MIT
package la

import "core:math"

// temporary
sparse_cg_solve :: proc(
	A: Block_Sparse_Matrix,
	x, b: Block_Vector,
	max_iter: int,
	rel_tol: f64,
	allocator := context.allocator,
) -> (converged: bool, iters: int, residual: f64) {
    n := len(x.values)

    r := block_vector_zero_clone(x)
    p := block_vector_zero_clone(x)
	Ap := block_vector_zero_clone(x)

	copy(r.values, b.values)
	gemv(A, x, r, a = -1.0, b = 1.0)

	copy(p.values, r.values)

    rsold := dot(r, r)
    b_norm := nrm2(b)

    if b_norm < 1e-15 {b_norm = 1.0}

    for iter := 0; iter < max_iter; iter += 1 {
		gemv(A, p, Ap, a = 1.0, b = 0.0)
		pAp := dot(p, Ap)

		if abs(pAp) < 1e-15 {
			residual = math.sqrt(rsold) / b_norm
			return false, iter, residual
		}

		alpha := rsold / pAp
		axpy(p, x, alpha)
		axpy(Ap, r, -alpha)

		rsnew := dot(r, r)
		residual = math.sqrt(rsnew) / b_norm

		if residual < rel_tol {
			return true, iter + 1, residual
		}

		beta := rsnew / rsold
		axpby(r, p, 1.0, beta)
		rsold = rsnew
	}

	residual = math.sqrt(rsold) / b_norm
	return false, max_iter, residual

}
