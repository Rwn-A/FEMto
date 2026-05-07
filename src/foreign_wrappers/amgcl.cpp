#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/backend/builtin.hpp>
#include <functional>
#include <memory>
#include <tuple>
#include <cstddef>
#include <cstdio>
#include <new>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct amgcl_solver amgcl_solver_t;

typedef enum {
    AMGCL_CG_SA         = 0,
    AMGCL_BICGSTAB_SA   = 1,
    AMGCL_FGMRES_SA     = 2,
    AMGCL_BICGSTAB_ILU0 = 3,
    AMGCL_FGMRES_ILUT   = 4,
} amgcl_solver_kind_t;

typedef struct {
    double               tolerance;
    int                  max_iters;
    int                  verbose;
    int                  block_size;
    int                  gmres_m;
    int                  coarse_enough;
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
using Range   = amgcl::iterator_range<const double *>;
using RangeX  = amgcl::iterator_range<double *>;

using Solver_CG_SA         = amgcl::make_solver<amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::gauss_seidel>, amgcl::solver::cg<Backend>>;
using Solver_BiCGStab_SA   = amgcl::make_solver<amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::gauss_seidel>, amgcl::solver::bicgstab<Backend>>;
using Solver_FGMRES_SA     = amgcl::make_solver<amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::gauss_seidel>, amgcl::solver::fgmres<Backend>>;
using Solver_BiCGStab_ILU0 = amgcl::make_solver<amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>, amgcl::solver::bicgstab<Backend>>;
using Solver_FGMRES_ILU0   = amgcl::make_solver<amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>, amgcl::solver::fgmres<Backend>>;

struct amgcl_solver {
    double   tolerance;
    int      max_iters;
    int      verbose;
    int      n;

    std::function<std::tuple<std::size_t, double>(Range, RangeX)> solve_fn;
    std::function<void()>                                          print_fn;
};

void amgcl_params_default(amgcl_params_t *p) {
    if (!p) return;
    p->tolerance     = 1e-8;
    p->max_iters     = 500;
    p->verbose       = 0;
    p->block_size    = 1;
    p->gmres_m       = 30;
    p->coarse_enough = 50;
    p->kind          = AMGCL_CG_SA;
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

    auto make_A = [&]() {
        return std::make_tuple(
            n,
            amgcl::make_iterator_range(row_ptr, row_ptr + n + 1),
            amgcl::make_iterator_range(col_ind, col_ind + row_ptr[n]),
            amgcl::make_iterator_range(values,  values  + row_ptr[n]));
    };

    try {
        switch (p.kind) {
        case AMGCL_CG_SA: {
            Solver_CG_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            auto s = std::make_shared<Solver_CG_SA>(make_A(), sp);
            h->solve_fn = [s](Range r, RangeX x) { return (*s)(r, x); };
            h->print_fn = [s]{ std::cout << s->precond() << std::endl; };
            break;
        }
        case AMGCL_BICGSTAB_SA: {
            Solver_BiCGStab_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            auto s = std::make_shared<Solver_BiCGStab_SA>(make_A(), sp);
            h->solve_fn = [s](Range r, RangeX x) { return (*s)(r, x); };
            h->print_fn = [s]{ std::cout << s->precond() << std::endl; };
            break;
        }
        case AMGCL_FGMRES_SA: {
            Solver_FGMRES_SA::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            sp.precond.coarsening.aggr.block_size = p.block_size;
            sp.precond.coarse_enough = p.coarse_enough;
            auto s = std::make_shared<Solver_FGMRES_SA>(make_A(), sp);
            h->solve_fn = [s](Range r, RangeX x) { return (*s)(r, x); };
            h->print_fn = [s]{ std::cout << s->precond() << std::endl; };
            break;
        }
        case AMGCL_BICGSTAB_ILU0: {
            Solver_BiCGStab_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            auto s = std::make_shared<Solver_BiCGStab_ILU0>(make_A(), sp);
            h->solve_fn = [s](Range r, RangeX x) { return (*s)(r, x); };
            h->print_fn = [s]{ std::cout << s->precond() << std::endl; };
            break;
        }
        case AMGCL_FGMRES_ILUT: {
            Solver_FGMRES_ILU0::params sp;
            sp.solver.tol = p.tolerance; sp.solver.maxiter = p.max_iters;
            sp.solver.M   = p.gmres_m;
            auto s = std::make_shared<Solver_FGMRES_ILU0>(make_A(), sp);
            h->solve_fn = [s](Range r, RangeX x) { return (*s)(r, x); };
            h->print_fn = [s]{ std::cout << s->precond() << std::endl; };
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

    if (p.verbose && h->print_fn)
        h->print_fn();

    return h;
}

int amgcl_solve(amgcl_solver_t *h, const double *rhs, double *x,
                int *out_iters, double *out_residual)
{
    if (!h || !h->solve_fn) return -1;

    std::size_t iters = 0;
    double      resid = 0.0;

    try {
        auto r = amgcl::make_iterator_range(rhs, rhs + h->n);
        auto X = amgcl::make_iterator_range(x,   x   + h->n);
        std::tie(iters, resid) = h->solve_fn(r, X);
    } catch (...) {
        return -1;
    }

    if (out_iters)    *out_iters    = static_cast<int>(iters);
    if (out_residual) *out_residual = resid;
    if (h->verbose)
        printf("amgcl: %d iters  residual = %e\n", static_cast<int>(iters), resid);

    return resid <= h->tolerance ? 0 : 1;
}

void amgcl_destroy(amgcl_solver_t *h) { delete h; }