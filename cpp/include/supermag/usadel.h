#ifndef SUPERMAG_USADEL_H
#define SUPERMAG_USADEL_H

#include "error.h"
#include "solver_options.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Usadel solver mode */
typedef enum {
    SUPERMAG_USADEL_LINEARIZED = 0,  /* Analytic linearized profile */
    SUPERMAG_USADEL_NONLINEAR  = 1   /* Self-consistent Newton iteration */
} supermag_usadel_mode_t;

/* Usadel diffusive-limit solver
 * Self-consistency loop for superconducting order parameter Delta(x)
 * with Kupriyanov-Lukichev boundary conditions at S/F interfaces.
 *
 * T: temperature (K), must be > 0.
 * mode: LINEARIZED or NONLINEAR.
 * opts: solver options (NULL = defaults). Uses: matsubara_max, omega_cut_factor,
 *       max_iter, conv_tol.
 *
 * Output units:
 *   Delta_out: superconducting gap Δ(x) (meV).
 *   x_out:     position (nm) spanning [−d_S, d_F]. */
int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    double T,
    supermag_usadel_mode_t mode,
    const supermag_solver_options_t *opts,
    int n_grid, double* Delta_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
