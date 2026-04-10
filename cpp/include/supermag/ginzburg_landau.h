#ifndef SUPERMAG_GINZBURG_LANDAU_H
#define SUPERMAG_GINZBURG_LANDAU_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: Ginzburg-Landau free energy functional minimization.
 * Vortex states, mixed-state configurations, domain structures near Tc. */
int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    double* psi_real, double* psi_imag
);

#ifdef __cplusplus
}
#endif

#endif
