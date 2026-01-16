// SPDX-FileCopyrightText: 2025 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Iterative methods for solving large sparse matrices.
 ----------------------------------------------------
*/
package la

import "core:log"

import "core:math"

// Temporary.
sparse_cg_solve :: proc(
	A: Sparse_Matrix,
	x, b: Vec,
	max_iter: int = 1000,
	tol: f64 = 1e-10,
	allocator := context.allocator,
) -> (
	iterations: int,
	residual: f64,
	converged: bool,
) {
	assert(len(x) == len(b))
	assert(len(x) == sparse_mat_rows(A))

	n := len(x)

	r := make(Vec, n, allocator)
	defer delete(r, allocator)

	p := make(Vec, n, allocator)
	defer delete(p, allocator)

	Ap := make(Vec, n, allocator)
	defer delete(Ap, allocator)

	copy(r, b)
	sparse_mv(A, x, r, a = -1.0, b = 1.0)
	copy(p, r)

	rsold := dot(r, r)
	b_norm := nrm2(b)

	if b_norm < 1e-15 {b_norm = 1.0}

	for iter := 0; iter < max_iter; iter += 1 {
		sparse_mv(A, p, Ap, a = 1.0, b = 0.0)
		pAp := dot(p, Ap)

		if abs(pAp) < 1e-15 {
			residual = math.sqrt(rsold) / b_norm
			return iter, residual, false
		}

		alpha := rsold / pAp

		axpy(p, x, alpha)
		axpy(Ap, r, -alpha)

		rsnew := dot(r, r)
		residual = math.sqrt(rsnew) / b_norm
		if residual < tol {
			return iter + 1, residual, true
		}

		beta := rsnew / rsold
		axpby(r, p, 1.0, beta)
		rsold = rsnew
	}
	residual = math.sqrt(rsold) / b_norm
	return max_iter, residual, false
}
