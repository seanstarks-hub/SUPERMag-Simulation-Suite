#ifndef SUPERMAG_EILENBERGER_H
#define SUPERMAG_EILENBERGER_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: Eilenberger clean-limit solver
 * Riccati parameterization for numerical stability. */
int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    int n_grid, double* f_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
