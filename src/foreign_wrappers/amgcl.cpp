#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/backend/builtin.hpp>
#include <memory>
#include <tuple>
#include <cstddef>
#include <cstdio>
#include <new>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct amgcl_solver amgcl_solver_t;

// Each value is a distinct (iterative solver, preconditioner) pair.
// Every combination here has its own C++ template instantiation.
typedef enum {
    // AMG preconditioned — smoothed aggregation hierarchy
    AMGCL_CG_SA        = 0,  // CG    + AMG(smoothed_aggregation, spai0)
    AMGCL_BICGSTAB_SA  = 1,  // BiCGStab + AMG(smoothed_aggregation, spai0)
    AMGCL_GMRES_SA     = 2,  // GMRES + AMG(smoothed_aggregation, spai0)
    AMGCL_FGMRES_SA    = 3,  // FGMRES + AMG(smoothed_aggregation, spai0)

    // AMG preconditioned — ruge-stuben hierarchy (no block_size)
    AMGCL_CG_RS        = 4,  // CG    + AMG(ruge_stuben, spai0)
    AMGCL_BICGSTAB_RS  = 5,  // BiCGStab + AMG(ruge_stuben, spai0)
    AMGCL_GMRES_RS     = 6,  // GMRES + AMG(ruge_stuben, spai0)
    AMGCL_FGMRES_RS    = 7,  // FGMRES + AMG(ruge_stuben, spai0)

    // ILU preconditioned — no hierarchy
    AMGCL_CG_ILU0      = 8,  // CG    + ILU(0) as preconditioner
    AMGCL_BICGSTAB_ILU0 = 9, // BiCGStab + ILU(0)
    AMGCL_GMRES_ILU0   = 10, // GMRES + ILU(0)
    AMGCL_FGMRES_ILU0  = 11, // FGMRES + ILU(0)
} amgcl_solver_kind_t;

typedef struct {
    double               tolerance;
    int                  max_iters;
    int                  verbose;
    int                  block_size;    // SA only: 1=scalar, 2=2D, 3=3D elasticity
    int                  gmres_m;       // GMRES/FGMRES restart parameter (default 30)
    int                  coarse_enough; // AMG: coarsest level size threshold (default 50)
    amgcl_solver_kind_t  kind;
} amgcl_params_t;

void            amgcl_params_default(amgcl_params_t *p);
amgcl_solver_t *amgcl_create(int n, const int *row_ptr, const int *col_ind,
                              const double *values, const amgcl_params_t *p);
int             amgcl_solve(amgcl_solver_t *s, const double *rhs, double *x,
                            int *out_iters, double *out_residual);
void            amgcl_destroy(amgcl_solver_t *s);

#ifdef __cplusplus
}
#endif

/* ======================================================================= */
/* Implementation                                                           */
/* ======================================================================= */

using Backend = amgcl::backend::builtin<double>;

// --- SA variants ---
using Solver_CG_SA = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
    amgcl::solver::cg<Backend>>;

using Solver_BiCGStab_SA = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
    amgcl::solver::bicgstab<Backend>>;

using Solver_GMRES_SA = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
    amgcl::solver::gmres<Backend>>;

using Solver_FGMRES_SA = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
    amgcl::solver::fgmres<Backend>>;

// --- RS variants ---
using Solver_CG_RS = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
    amgcl::solver::cg<Backend>>;

using Solver_BiCGStab_RS = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
    amgcl::solver::bicgstab<Backend>>;

using Solver_GMRES_RS = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
    amgcl::solver::gmres<Backend>>;

using Solver_FGMRES_RS = amgcl::make_solver
    amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::spai0>,
    amgcl::solver::fgmres<Backend>>;

// --- ILU0 variants ---
using Solver_CG_ILU0 = amgcl::make_solver
    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
    amgcl::solver::cg<Backend>>;

using Solver_BiCGStab_ILU0 = amgcl::make_solver
    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
    amgcl::solver::bicgstab<Backend>>;

using Solver_GMRES_ILU0 = amgcl::make_solver
    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
    amgcl::solver::gmres<Backend>>;

using Solver_FGMRES_ILU0 = amgcl::make_solver
    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
    amgcl::solver::fgmres<Backend>>;

struct amgcl_solver {
    double              tolerance;
    int                 max_iters;
    int                 verbose;
    int                 n;
    amgcl_solver_kind_t kind;

    // Only one will be non-null at runtime
    std::unique_ptr<Solver_CG_SA>        cg_sa;
    std::unique_ptr<Solver_BiCGStab_SA>  bicgstab_sa;
    std::unique_ptr<Solver_GMRES_SA>     gmres_sa;
    std::unique_ptr<Solver_FGMRES_SA>    fgmres_sa;
    std::unique_ptr<Solver_CG_RS>        cg_rs;
    std::unique_ptr<Solver_BiCGStab_RS>  bicgstab_rs;
    std::unique_ptr<Solver_GMRES_RS>     gmres_rs;
    std::unique_ptr<Solver_FGMRES_RS>    fgmres_rs;
    std::unique_ptr<Solver_CG_ILU0>      cg_ilu0;
    std::unique_ptr<Solver_BiCGStab_ILU0> bicgstab_ilu0;
    std::unique_ptr<Solver_GMRES_ILU0>   gmres_ilu0;
    std::unique_ptr<Solver_FGMRES_ILU0>  fgmres_ilu0;
};

void amgcl_params_default(amgcl_params_t *p) {
    if (!p) return;
    p->tolerance    = 1e-8;
    p->max_iters    = 500;
    p->verbose      = 0;
    p->block_size   = 1;
    p->gmres_m      = 30;
    p->coarse_enough = 50;
    p->kind         = AMGCL_CG_SA;
}

amgcl_solver_t *amgcl_create(int n, const int *row_ptr, const int *col_ind,
                              const double *values, const amgcl_params_t *params)
{
    amgcl_params_t p;
    amgcl_params_default(&p);
    if (params) p = *params;
    if (p.block_size < 1) p.block_size = 1;
    if (p.block_size > 4) p.block_size = 4;
    if (p.gmres_m    < 1) p.gmres_m    = 30;

    amgcl_solver *h = new(std::nothrow) amgcl_solver();
    if (!h) return nullptr;
    h->tolerance = p.tolerance;
    h->max_iters = p.max_iters;
    h->verbose   = p.verbose;
    h->n         = n;
    h->kind      = p.kind;

    // Helper lambda — avoids repeating the tuple construction 12 times
    auto make_A = [&]() {
        return std::make_tuple(
            n,
            amgcl::make_iterator_range(row_ptr, row_ptr + n + 1),
            amgcl::make_iterator_range(col_ind, col_ind + row_ptr[n]),
            amgcl::make_iterator_range(values,  values  + row_ptr[n]));
    };

    try {
        switch (p.kind) {
        // --- SA ---
        case AMGCL_CG_SA: {
            Solver_CG_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            h->cg_sa.reset(new Solver_CG_SA(make_A(), sp));
            break;
        }
        case AMGCL_BICGSTAB_SA: {
            Solver_BiCGStab_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            h->bicgstab_sa.reset(new Solver_BiCGStab_SA(make_A(), sp));
            break;
        }
        case AMGCL_GMRES_SA: {
            Solver_GMRES_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            h->gmres_sa.reset(new Solver_GMRES_SA(make_A(), sp));
            break;
        }
        case AMGCL_FGMRES_SA: {
            Solver_FGMRES_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            h->fgmres_sa.reset(new Solver_FGMRES_SA(make_A(), sp));
            break;
        }
        // --- RS ---
        case AMGCL_CG_RS: {
            Solver_CG_RS::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarse_enough = p.coarse_enough;
            h->cg_rs.reset(new Solver_CG_RS(make_A(), sp));
            break;
        }
        case AMGCL_BICGSTAB_RS: {
            Solver_BiCGStab_RS::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarse_enough = p.coarse_enough;
            h->bicgstab_rs.reset(new Solver_BiCGStab_RS(make_A(), sp));
            break;
        }
        case AMGCL_GMRES_RS: {
            Solver_GMRES_RS::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            sp.precond.coarse_enough = p.coarse_enough;
            h->gmres_rs.reset(new Solver_GMRES_RS(make_A(), sp));
            break;
        }
        case AMGCL_FGMRES_RS: {
            Solver_FGMRES_RS::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            sp.precond.coarse_enough = p.coarse_enough;
            h->fgmres_rs.reset(new Solver_FGMRES_RS(make_A(), sp));
            break;
        }
        // --- ILU0 ---
        case AMGCL_CG_ILU0: {
            Solver_CG_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            h->cg_ilu0.reset(new Solver_CG_ILU0(make_A(), sp));
            break;
        }
        case AMGCL_BICGSTAB_ILU0: {
            Solver_BiCGStab_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            h->bicgstab_ilu0.reset(new Solver_BiCGStab_ILU0(make_A(), sp));
            break;
        }
        case AMGCL_GMRES_ILU0: {
            Solver_GMRES_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            h->gmres_ilu0.reset(new Solver_GMRES_ILU0(make_A(), sp));
            break;
        }
        case AMGCL_FGMRES_ILU0: {
            Solver_FGMRES_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            h->fgmres_ilu0.reset(new Solver_FGMRES_ILU0(make_A(), sp));
            break;
        }
        default:
            delete h;
            return nullptr;
        }
    } catch (...) {
        delete h;
        return nullptr;
    }

    if (p.verbose) {
        // Print preconditioner info — each branch knows its active pointer
        switch (p.kind) {
        case AMGCL_CG_SA:         std::cout << h->cg_sa->precond()        << std::endl; break;
        case AMGCL_BICGSTAB_SA:   std::cout << h->bicgstab_sa->precond()  << std::endl; break;
        case AMGCL_GMRES_SA:      std::cout << h->gmres_sa->precond()      << std::endl; break;
        case AMGCL_FGMRES_SA:     std::cout << h->fgmres_sa->precond()     << std::endl; break;
        case AMGCL_CG_RS:         std::cout << h->cg_rs->precond()        << std::endl; break;
        case AMGCL_BICGSTAB_RS:   std::cout << h->bicgstab_rs->precond()  << std::endl; break;
        case AMGCL_GMRES_RS:      std::cout << h->gmres_rs->precond()      << std::endl; break;
        case AMGCL_FGMRES_RS:     std::cout << h->fgmres_rs->precond()     << std::endl; break;
        case AMGCL_CG_ILU0:       std::cout << h->cg_ilu0->precond()      << std::endl; break;
        case AMGCL_BICGSTAB_ILU0: std::cout << h->bicgstab_ilu0->precond()<< std::endl; break;
        case AMGCL_GMRES_ILU0:    std::cout << h->gmres_ilu0->precond()    << std::endl; break;
        case AMGCL_FGMRES_ILU0:   std::cout << h->fgmres_ilu0->precond()   << std::endl; break;
        default: break;
        }
    }
    return h;
}

// Macro to cut the repetition in amgcl_solve — each branch is identical
// except for which unique_ptr field it dereferences.
#define SOLVE_WITH(field)                                          \
    if (!h->field) return -1;                                      \
    std::tie(iters, resid) = (*h->field)(r, X);                    \
    break;

int amgcl_solve(amgcl_solver_t *h, const double *rhs, double *x,
                int *out_iters, double *out_residual)
{
    if (!h) return -1;
    std::size_t iters = 0;
    double      resid = 0.0;

    try {
        auto r = amgcl::make_iterator_range(rhs, rhs + h->n);
        auto X = amgcl::make_iterator_range(x,   x   + h->n);

        switch (h->kind) {
        case AMGCL_CG_SA:         SOLVE_WITH(cg_sa)
        case AMGCL_BICGSTAB_SA:   SOLVE_WITH(bicgstab_sa)
        case AMGCL_GMRES_SA:      SOLVE_WITH(gmres_sa)
        case AMGCL_FGMRES_SA:     SOLVE_WITH(fgmres_sa)
        case AMGCL_CG_RS:         SOLVE_WITH(cg_rs)
        case AMGCL_BICGSTAB_RS:   SOLVE_WITH(bicgstab_rs)
        case AMGCL_GMRES_RS:      SOLVE_WITH(gmres_rs)
        case AMGCL_FGMRES_RS:     SOLVE_WITH(fgmres_rs)
        case AMGCL_CG_ILU0:       SOLVE_WITH(cg_ilu0)
        case AMGCL_BICGSTAB_ILU0: SOLVE_WITH(bicgstab_ilu0)
        case AMGCL_GMRES_ILU0:    SOLVE_WITH(gmres_ilu0)
        case AMGCL_FGMRES_ILU0:   SOLVE_WITH(fgmres_ilu0)
        default: return -1;
        }
    } catch (...) {
        return -1;
    }

    if (out_iters)    *out_iters    = static_cast<int>(iters);
    if (out_residual) *out_residual = resid;
    if (h->verbose)
        printf("amgcl: %d iters  residual = %e  kind = %d\n",
               static_cast<int>(iters), resid, (int)h->kind);

    return resid <= h->tolerance ? 0 : 1;
}

#undef SOLVE_WITH

void amgcl_destroy(amgcl_solver_t *h) { delete h; }