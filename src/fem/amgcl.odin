// bindings for the amgcl library for solving sparse linear systems.
package fem

import "core:c"
import "core:log"

when ODIN_OS == .Windows {
	foreign import amgcl_lib "../../.build/amgcl.lib"
} else {
	foreign import amgcl_lib "../../.build/amgcl.a"
}

Solver :: distinct rawptr

Solver_Kind :: enum c.int {
	CG_SA         = 0,
	BiCGStab_SA   = 1,
	GMRES_SA      = 2,
	FGMRES_SA     = 3,
	CG_RS         = 4,
	BiCGStab_RS   = 5,
	GMRES_RS      = 6,
	FGMRES_RS     = 7,
	CG_ILU0       = 8,
	BiCGStab_ILU0 = 9,
	GMRES_ILU0    = 10,
	FGMRES_ILU0   = 11,
}

Params :: struct {
	tolerance:     f64,
	max_iters:     c.int,
	verbose:       c.int,
	block_size:    c.int,
	gmres_m:       c.int,
	coarse_enough: c.int,
	kind:          Solver_Kind,
}

@(default_calling_convention = "c", link_prefix = "amgcl_")
foreign amgcl_lib {
	params_default :: proc(p: ^Params) ---
	create :: proc(n: c.int, row_ptr, col_ind: [^]c.int, values: [^]f64, p: ^Params) -> Solver ---
	solve :: proc(s: Solver, rhs, x: [^]f64, out_iters: ^c.int, out_residual: ^f64) -> c.int ---
	destroy :: proc(s: Solver) ---
}

create_solver :: proc(
	n: int,
	row_ptr: []c.int,
	col_ind: []c.int,
	values: []f64,
	kind: Solver_Kind,
	block_size: int = 1,
	tol: f64 = 1e-8,
	max_iters: int = 500,
	gmres_m: int = 30,
	coarse_enough: int = 50,
	verbose: bool = false,
) -> Solver {
	p: Params
	params_default(&p)
	p.kind = kind
	p.tolerance = tol
	p.max_iters = c.int(max_iters)
	p.block_size = c.int(block_size)
	p.gmres_m = c.int(gmres_m)
	p.coarse_enough = c.int(coarse_enough)
	p.verbose = 1 if verbose else 0
	return create(c.int(n), raw_data(row_ptr), raw_data(col_ind), raw_data(values), &p)
}

Solve_Result :: struct {
	iters:     int,
	residual:  f64,
	converged: bool,
}

solve_system :: proc(s: Solver, rhs: []f64, x: []f64) -> Solve_Result {
	iters: c.int
	residual: f64
	status := solve(s, raw_data(rhs), raw_data(x), &iters, &residual)
	return {int(iters), residual, status == 0}
}

// sparse_solve is the one-shot convenience path — creates, solves, destroys.
sparse_solve :: proc(
	A: Sparse_Matrix,
	x, b: Vector,
	kind: Solver_Kind = .GMRES_ILU0,
	tol: f64 = 1e-8,
	max_iters: int = 500,
	gmres_m: int = 30,
	coarse_enough: int = 50,
) -> (
	r: Solve_Result,
	ok: bool,
) {
	s := create_solver(
		sparse_matrix_rows(A),
		A.row_ptrs,
		A.columns,
		A.values,
		kind,
		1,
		tol,
		max_iters,
		gmres_m,
		coarse_enough,
	)
	if s == nil {
		log.error("amgcl: solver creation failed")
		return {}, false
	}
	defer destroy(s)
	r = solve_system(s, b, x)
	ok = r.converged
	return
}
