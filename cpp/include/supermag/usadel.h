#ifndef SUPERMAG_USADEL_H
#define SUPERMAG_USADEL_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Usadel diffusive-limit solver
 * Self-consistency loop for superconducting order parameter Delta(x)
 * with Kupriyanov-Lukichev boundary conditions at S/F interfaces. */
int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    int n_grid, double* Delta_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
