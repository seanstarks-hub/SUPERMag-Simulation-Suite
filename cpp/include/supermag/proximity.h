#ifndef SUPERMAG_PROXIMITY_H
#define SUPERMAG_PROXIMITY_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Equation model selection */
typedef enum {
    SUPERMAG_MODEL_THIN_S   = 0,
    SUPERMAG_MODEL_FOMINOV  = 1
} supermag_model_t;

/* Phase selection (kernel dispatch) */
typedef enum {
    SUPERMAG_PHASE_ZERO = 0,   /* coth kernel, 0-junction */
    SUPERMAG_PHASE_PI   = 1    /* tanh kernel, pi-junction */
} supermag_phase_t;

/* Depairing channels — all additive */
typedef struct {
    double ag;          /* Abrikosov-Gor'kov magnetic impurity pair-breaking */
    double zeeman;      /* Zeeman splitting */
    double orbital;     /* Orbital pair-breaking */
    double spin_orbit;  /* Spin-orbit scattering */
} supermag_depairing_t;

/* Full parameter set for a proximity calculation */
typedef struct {
    double Tc0;         /* Bulk superconductor critical temperature (K) */
    double d_S;         /* Superconductor layer thickness (nm) */
    double d_F;         /* Ferromagnet layer thickness (nm) */
    double xi_S;        /* Superconductor coherence length (nm) */
    double xi_F;        /* Ferromagnet coherence length (nm) */
    double gamma;       /* Interface transparency parameter (dimensionless) */
    double gamma_B;     /* Interface barrier parameter (dimensionless, Fominov model) */
    double E_ex;        /* Exchange energy (meV) */
    double D_F;         /* Diffusion coefficient in F layer (nm^2/ps) */
    supermag_model_t model;
    supermag_phase_t phase;
} supermag_proximity_params_t;

/* Total depairing parameter (sum of all channels) */
double supermag_depairing_total(const supermag_depairing_t *dp);

/**
 * Solve for Tc at a single d_F value.
 * Internally picks kernel, builds equation F(T), runs Brent root-finder.
 */
int supermag_proximity_solve_tc(
    const supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,  /* NULL if no depairing */
    double *tc_out
);

/**
 * Solve for Tc across an array of d_F values.
 * @param d_F_array   Array of ferromagnet thicknesses (nm)
 * @param n_dF        Number of elements
 * @param tc_out      Preallocated output array of length n_dF
 */
int supermag_proximity_solve_tc_batch(
    const supermag_proximity_params_t *params,
    const double *d_F_array,
    int n_dF,
    const supermag_depairing_t *depairing,
    double *tc_out
);

/**
 * Compute pair amplitude profile F(x) in the F layer.
 * F(x) = exp(-x/xi_F) * cos(x/xi_F)   for phase=0
 * F(x) = exp(-x/xi_F) * sin(x/xi_F)   for phase=pi
 *
 * @param d_F       Ferromagnet thickness (nm)
 * @param xi_F      Ferromagnet coherence length (nm)
 * @param phase     SUPERMAG_PHASE_ZERO or SUPERMAG_PHASE_PI
 * @param n_points  Number of grid points
 * @param x_out     Output: x positions (nm), preallocated
 * @param F_out     Output: F(x) values, preallocated
 */
int supermag_proximity_pair_amplitude(
    double d_F, double xi_F, supermag_phase_t phase,
    int n_points, double *x_out, double *F_out
);

#ifdef __cplusplus
}
#endif

#endif
