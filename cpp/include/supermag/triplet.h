#ifndef SUPERMAG_TRIPLET_H
#define SUPERMAG_TRIPLET_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Triplet solver mode */
typedef enum {
    SUPERMAG_TRIPLET_PHENOMENOLOGICAL = 0,  /* sin(α)·exp(-x/ξ_N) model */
    SUPERMAG_TRIPLET_DIFFUSIVE        = 1   /* Coupled diffusive f_0/f_1 equations */
} supermag_triplet_mode_t;

/* Spin-triplet superconductivity solver.
 * Long-range triplet correlations from inhomogeneous magnetization.
 * Solves linearized coupled Usadel equations for singlet (f_0) and
 * long-range triplet (f_1) components.
 *
 * E_ex_per_layer: per-layer exchange energy (meV), length n_layers. NULL → use global.
 * D_per_layer:    per-layer diffusion coeff (nm^2/ps), length n_layers. NULL → default.
 * xi_F: short-range (ferromagnetic) coherence length (nm). If <= 0 uses 1.0.
 * xi_N: long-range (triplet) coherence length (nm).        If <= 0 uses 10.0.
 * T:    temperature (K), must be > 0.
 * mode: PHENOMENOLOGICAL or DIFFUSIVE. */
int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    const double* E_ex_per_layer, const double* D_per_layer,
    double xi_F, double xi_N, double T,
    supermag_triplet_mode_t mode,
    int n_grid, double* f_triplet_out, double* x_out
);

#ifdef __cplusplus
}
#endif

#endif
