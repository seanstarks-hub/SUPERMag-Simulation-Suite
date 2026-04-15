#ifndef SUPERMAG_GINZBURG_LANDAU_H
#define SUPERMAG_GINZBURG_LANDAU_H

#include "error.h"
#include "solver_options.h"

#ifdef __cplusplus
extern "C" {
#endif

/* GL solver mode */
typedef enum {
    SUPERMAG_GL_SCALAR = 0,  /* Standard scalar GL (no gauge field) */
    SUPERMAG_GL_GAUGE  = 1   /* Gauge-covariant with Peierls phases */
} supermag_gl_mode_t;

/* Ginzburg-Landau free energy functional minimization.
 * Vortex states, mixed-state configurations, domain structures near Tc.
 * Uses κ for the coherence length: ξ² = 1/(2κ²) in GL units.
 *
 * mode: SCALAR ignores H_applied; GAUGE includes Peierls phases.
 * H_applied: uniform applied magnetic field (GL units). Ignored in SCALAR mode.
 * opts: solver options (NULL = defaults). Uses: max_steps, conv_tol.
 * If psi_real[0] != 0.0, the initial condition is preserved (user-supplied).
 *
 * Output units:
 *   psi_real, psi_imag: GL order parameter (GL units).
 *   Equilibrium: |ψ|² = −α/β for uniform state (α < 0). */
int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    supermag_gl_mode_t mode, double H_applied,
    const supermag_solver_options_t *opts,
    double* psi_real, double* psi_imag
);

#ifdef __cplusplus
}
#endif

#endif
