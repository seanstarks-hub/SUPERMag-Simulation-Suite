#ifndef SUPERMAG_TRIPLET_H
#define SUPERMAG_TRIPLET_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Spin-triplet superconductivity solver.
 * Long-range triplet correlations from inhomogeneous magnetization.
 * Spin-active interface boundary conditions.
 *
 * xi_F: short-range (ferromagnetic) coherence length (nm). If <= 0 uses 1.0.
 * xi_N: long-range (triplet) coherence length (nm).        If <= 0 uses 10.0. */
int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    double xi_F, double xi_N,
    int n_grid, double* f_triplet_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
