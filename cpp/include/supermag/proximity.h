#ifndef SUPERMAG_PROXIMITY_H
#define SUPERMAG_PROXIMITY_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Equation model selection (physics model for the Tc equation) */
typedef enum {
    SUPERMAG_MODEL_THIN_S       = 0,   /* Thin-S limit digamma equation */
    SUPERMAG_MODEL_FOMINOV      = 1,   /* Fominov PRB 66 single-mode */
    SUPERMAG_MODEL_FOMINOV_MULTI = 2   /* Fominov multimode (multiple sub-gap channels) */
} supermag_model_t;

/* Geometry selection (orthogonal to equation model) */
typedef enum {
    SUPERMAG_GEOM_BILAYER  = 0,   /* Standard S/F bilayer */
    SUPERMAG_GEOM_TRILAYER = 1,   /* S/N/F trilayer (geom_params → trilayer_params_t) */
    SUPERMAG_GEOM_GRADED   = 2,   /* Graded exchange (geom_params → graded_params_t) */
    SUPERMAG_GEOM_DOMAINS  = 3    /* Magnetic domains (geom_params → domain_params_t) */
} supermag_geometry_t;

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

/* Spin-active interface parameters */
typedef struct {
    double mixing_angle;  /* Spin-mixing angle θ (radians) */
    double polarization;  /* Interface spin polarization P ∈ [0,1] */
} supermag_spin_active_t;

/* Interface roughness */
typedef struct {
    double delta_roughness;  /* RMS roughness amplitude (nm) */
} supermag_roughness_t;

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
    double D_S;         /* Diffusion coefficient in S layer (nm^2/ps) */
    supermag_model_t model;
    supermag_phase_t phase;
    supermag_geometry_t geometry;
    const void *geom_params;         /* NULL for BILAYER; cast to trilayer/graded/domain_params_t */
    const supermag_spin_active_t *spin_active;  /* NULL if no spin-active interface */
} supermag_proximity_params_t;

/* Total depairing parameter (sum of all channels) */
double supermag_depairing_total(const supermag_depairing_t *dp);

/* Input struct for physics-based depairing computation.
 * Units: Gamma_s/Gamma_so in meV, H in Tesla, D in nm^2/ps,
 *        thickness in nm, Tc0 in Kelvin. */
typedef struct {
    double Gamma_s;     /* Spin-flip scattering rate (meV) */
    double H;           /* Applied magnetic field (T) */
    double D;           /* Diffusion coefficient (nm^2/ps) */
    double thickness;   /* Film thickness (nm) */
    double Gamma_so;    /* Spin-orbit scattering rate (meV) */
    double Tc0;         /* Bulk critical temperature (K) */
} supermag_depairing_input_t;

/**
 * Compute physics-based depairing channels from laboratory-measurable inputs.
 * Fills all four channels of the output depairing struct:
 *   AG      = Gamma_s / (2 * kB * Tc0)          [EQ-7A]
 *   Zeeman  = (muB * H)^2 / (2*pi*kB*Tc0)^2     [EQ-7B]
 *   Orbital = D*(e*H)^2*thickness^2 / (3*hbar^2*(2*pi*kB*Tc0))  [EQ-7C]
 *   SO      = Gamma_so / (2 * kB * Tc0)          [EQ-7D]
 */
int supermag_depairing_compute(
    const supermag_depairing_input_t *input,
    supermag_depairing_t *output
);

/* (Individual depairing and optimizer functions moved to depairing.h and optimizer.h) */

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

/* ── S/N/F Trilayer Kernel ──────────────────────────────────────── */

/* Parameters for the normal-metal interlayer in an S/N/F trilayer */
typedef struct {
    double d_N;         /* Normal metal thickness (nm) */
    double xi_N;        /* Normal metal coherence length (nm) */
    double R_B;         /* Interface barrier resistance (Ohm·nm^2), 0 = transparent */
} supermag_trilayer_params_t;

/**
 * Compute the effective proximity kernel for an S/N/F trilayer.  [EQ-13]
 * Composes the F-layer kernel with N-layer propagation via
 * continued-fraction formula.
 *
 * @param d_F       Ferromagnet thickness (nm)
 * @param xi_F      Ferromagnet coherence length (nm)
 * @param tri       Trilayer parameters (N-layer)
 * @param phase     SUPERMAG_PHASE_ZERO or SUPERMAG_PHASE_PI
 * @param K_real    Output: real part of effective kernel
 * @param K_imag    Output: imaginary part of effective kernel
 */
int supermag_proximity_kernel_snf(
    double d_F, double xi_F,
    const supermag_trilayer_params_t *tri,
    supermag_phase_t phase,
    double *K_real, double *K_imag
);

/* ── Graded Ferromagnet Kernel ──────────────────────────────────── */

/* Exchange energy profile across the ferromagnet */
typedef enum {
    SUPERMAG_GRADE_LINEAR      = 0,
    SUPERMAG_GRADE_EXPONENTIAL = 1,
    SUPERMAG_GRADE_STEP        = 2
} supermag_grade_profile_t;

/* Parameters for a graded (non-uniform) exchange energy profile */
typedef struct {
    double E_ex_surface;    /* Exchange energy at S/F interface (meV) */
    double E_ex_bulk;       /* Exchange energy at F vacuum surface (meV) */
    supermag_grade_profile_t profile;
    int n_slabs;            /* Number of slabs for discretization (>= 1) */
} supermag_graded_params_t;

/**
 * Compute the effective proximity kernel for a graded ferromagnet.  [EQ-14]
 * Slices the F layer into n_slices sub-layers with position-dependent
 * exchange energy, then cascades transfer matrices from vacuum inward.
 *
 * @param d_F         Total ferromagnet thickness (nm)
 * @param xi_F_ref    Reference coherence length at E_ex_left (nm)
 * @param grade       Graded parameters
 * @param phase       SUPERMAG_PHASE_ZERO or SUPERMAG_PHASE_PI
 * @param K_real      Output: real part of effective kernel
 * @param K_imag      Output: imaginary part of effective kernel
 */
int supermag_proximity_kernel_graded(
    double d_F, double xi_F_ref,
    const supermag_graded_params_t *grade,
    supermag_phase_t phase,
    double *K_real, double *K_imag
);

/* ── Magnetic Domain Kernel ─────────────────────────────────────── */

/* Parameters for alternating magnetic domains in the F layer */
typedef struct {
    double domain_width;    /* Width of each domain (nm) */
    double domain_wall;     /* Width of domain wall transition region (nm). 0 = sharp. */
} supermag_domain_params_t;

/**
 * Compute the effective proximity kernel for a domain-structured
 * ferromagnet with alternating ±E_ex magnetization.  [EQ-15]
 *
 * @param d_F       Total ferromagnet thickness (= d_domain * n_domains)
 * @param xi_F      Ferromagnet coherence length (nm)
 * @param E_ex      Exchange energy magnitude (meV)
 * @param dom       Domain parameters
 * @param phase     SUPERMAG_PHASE_ZERO or SUPERMAG_PHASE_PI
 * @param K_real    Output: real part of effective kernel
 * @param K_imag    Output: imaginary part of effective kernel
 */
int supermag_proximity_kernel_domains(
    double d_F, double xi_F, double E_ex,
    const supermag_domain_params_t *dom,
    supermag_phase_t phase,
    double *K_real, double *K_imag
);

/* Optimizer / inverse / fit functions moved to optimizer.h */

#ifdef __cplusplus
}
#endif

#endif
