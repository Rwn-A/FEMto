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
) -> (
	converged: bool,
	iters: int,
	residual: f64,
) {
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

// Full transparancy this was AI generated.
sparse_bicgstab_solve :: proc(
	A: Block_Sparse_Matrix,
	x, b: Block_Vector,
	max_iter: int,
	rel_tol: f64,
	allocator := context.allocator,
) -> (
	converged: bool,
	iters: int,
	residual: f64,
) {
	n := len(x.values)
	r     := block_vector_zero_clone(x)
	r_hat := block_vector_zero_clone(x)
	p     := block_vector_zero_clone(x)
	v     := block_vector_zero_clone(x)
	s     := block_vector_zero_clone(x)
	t     := block_vector_zero_clone(x)

	copy(r.values, b.values)
	gemv(A, x, r, a = -1.0, b = 1.0)

	copy(r_hat.values, r.values)

	b_norm := nrm2(b)
	if b_norm < 1e-15 { b_norm = 1.0 }

	rho_old : f64 = 1.0
	alpha   : f64 = 1.0
	omega   : f64 = 1.0


	for iter := 0; iter < max_iter; iter += 1 {
		rho_new := dot(r_hat, r)

		if abs(rho_new) < 1e-300 {
			residual = nrm2(r) / b_norm
			return false, iter, residual
		}

		beta := (rho_new / rho_old) * (alpha / omega)


		axpy(v, p, -omega)
		axpby(r, p, 1.0, beta)

		gemv(A, p, v, a = 1.0, b = 0.0)

		r_hat_v := dot(r_hat, v)
		if abs(r_hat_v) < 1e-300 {
			residual = nrm2(r) / b_norm
			return false, iter, residual
		}

		alpha = rho_new / r_hat_v

		copy(s.values, r.values)
		axpy(v, s, -alpha)

		s_norm := nrm2(s)
		if s_norm / b_norm < rel_tol {
			axpy(p, x, alpha)
			residual = s_norm / b_norm
			return true, iter + 1, residual
		}

		gemv(A, s, t, a = 1.0, b = 0.0)

		t_dot_t := dot(t, t)
		if abs(t_dot_t) < 1e-300 {
			axpy(p, x, alpha)
			residual = s_norm / b_norm
			return false, iter, residual
		}

		omega = dot(t, s) / t_dot_t


		axpy(p, x, alpha)
		axpy(s, x, omega)

		copy(r.values, s.values)
		axpy(t, r, -omega)

		residual = nrm2(r) / b_norm
		if residual < rel_tol {
			return true, iter + 1, residual
		}

		if abs(omega) < 1e-300 {
			return false, iter, residual
		}

		rho_old = rho_new
	}

	return false, max_iter, residual
}
