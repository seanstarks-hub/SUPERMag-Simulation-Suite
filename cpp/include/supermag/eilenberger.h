#ifndef SUPERMAG_EILENBERGER_H
#define SUPERMAG_EILENBERGER_H

#include "error.h"
#include "solver_options.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Eilenberger clean-limit solver
 * Riccati parameterization for numerical stability.
 *
 * T: temperature (K), must be > 0.
 * opts: solver options (NULL = defaults). Uses: matsubara_max, omega_cut_factor.
 *
 * Output units:
 *   f_out: anomalous Green's function |f(x)|, dimensionless in [0, 1].
 *   x_out: position (nm) spanning [−d_S, d_F]. */
int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    double T,
    const supermag_solver_options_t *opts,
    int n_grid, double* f_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
