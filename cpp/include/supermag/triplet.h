#ifndef SUPERMAG_TRIPLET_H
#define SUPERMAG_TRIPLET_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: Spin-triplet superconductivity solver.
 * Long-range triplet correlations from inhomogeneous magnetization.
 * Spin-active interface boundary conditions. */
int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    int n_grid, double* f_triplet_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
