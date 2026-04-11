#ifndef SUPERMAG_EILENBERGER_H
#define SUPERMAG_EILENBERGER_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Eilenberger clean-limit solver
 * Riccati parameterization for numerical stability.
 *
 * T: temperature (K).  If T <= 0 the solver uses 0.5*Tc0. */
int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    double T,
    int n_grid, double* f_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
