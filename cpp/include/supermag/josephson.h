#ifndef SUPERMAG_JOSEPHSON_H
#define SUPERMAG_JOSEPHSON_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Josephson current-phase relation (CPR) for S/F/S junctions.
 * 0-pi transition detection as function of F-layer thickness and temperature.
 *
 * Tc0: bulk superconductor Tc (K).  Must be > 0.
 * gamma_B: interface barrier parameter (dimensionless). 0 = transparent.
 * phase_arr: if non-NULL, user-supplied phase grid (const, not overwritten);
 *            if NULL, auto-generated uniform grid φ ∈ [0, 2π).
 * current_out: output array of I(φ) (normalized to max = 1).
 * Ic_out: if non-NULL, filled with absolute critical current. */
int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
    double gamma_B,
    int n_phases, const double* phase_arr,
    double* current_out,
    double* Ic_out
);

#ifdef __cplusplus
}
#endif

#endif
