#ifndef SUPERMAG_USADEL_H
#define SUPERMAG_USADEL_H

#include "error.h"

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
 * mode: LINEARIZED or NONLINEAR. */
int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    double T,
    supermag_usadel_mode_t mode,
    int n_grid, double* Delta_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
