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
	FGMRES_SA     = 2,
	BiCGStab_ILU0 = 3,
	FGMRES_ILU0   = 4,
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

Solver_Options :: struct {
	kind:          Solver_Kind,
	tol:           f64,
	max_iters:     int,
	gmres_m:       int,
	coarse_enough: int,
	block_size:    int,
	verbose:       bool,
}

DEFAULT_SOLVER_OPTIONS :: Solver_Options {
	kind          = .FGMRES_SA,
	tol           = 1e-8,
	max_iters     = 500,
	gmres_m       = 30,
	coarse_enough = 50,
	block_size    = 1,
	verbose       = false,
}

create_solver :: proc(
	n: int,
	row_ptr: []c.int,
	col_ind: []c.int,
	values: []f64,
	opts: Solver_Options = DEFAULT_SOLVER_OPTIONS,
) -> Solver {
	p: Params
	params_default(&p)
	p.kind = opts.kind
	p.tolerance = opts.tol
	p.max_iters = c.int(opts.max_iters)
	p.gmres_m = c.int(opts.gmres_m)
	p.coarse_enough = c.int(opts.coarse_enough)
	p.block_size = c.int(opts.block_size)
	p.verbose = 1 if opts.verbose else 0
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


sparse_solve :: proc(
	A: Sparse_Matrix,
	x, b: Vector,
	opts: Solver_Options = DEFAULT_SOLVER_OPTIONS,
) -> (
	r: Solve_Result,
	ok: bool,
) {
	s := create_solver(sparse_matrix_rows(A), A.row_ptrs, A.columns, A.values, opts)
	if s == nil {
		log.error("amgcl: solver creation failed")
		return {}, false
	}
	defer destroy(s)
	r = solve_system(s, b, x)
	ok = r.converged
	return
}
