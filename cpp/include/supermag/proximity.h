#ifndef SUPERMAG_PROXIMITY_H
#define SUPERMAG_PROXIMITY_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute pair amplitude F(x) in the ferromagnet layer of an S/F bilayer.
 *
 * F(x) = F0 * exp(-x / xi_F) * cos(x / xi_F)
 *
 * @param F0       Pair amplitude at the S/F interface (dimensionless)
 * @param xi_F     Ferromagnet coherence length (nm)
 * @param d_F      Ferromagnet thickness (nm)
 * @param n_points Number of grid points
 * @param x_out    Output: x positions (nm), must be pre-allocated with n_points elements
 * @param F_out    Output: F(x) values, must be pre-allocated with n_points elements
 * @return         SUPERMAG_OK on success, error code otherwise
 */
int supermag_proximity_pair_amplitude(
    double F0, double xi_F, double d_F,
    int n_points, double* x_out, double* F_out
);

/**
 * Compute critical temperature Tc as a function of ferromagnet thickness d_F.
 *
 * Uses single-mode Buzdin approximation for S/F bilayer.
 * Tc(d_F) shows oscillatory behavior due to FFLO-like pair amplitude oscillation.
 *
 * @param Tc0      Bulk superconductor critical temperature (K)
 * @param d_S      Superconductor thickness (nm)
 * @param xi_S     Superconductor coherence length (nm)
 * @param xi_F     Ferromagnet coherence length (nm)
 * @param E_ex     Exchange energy (meV)
 * @param d_F_arr  Input: array of d_F values (nm)
 * @param n_dF     Number of d_F values
 * @param Tc_out   Output: Tc values (K), must be pre-allocated with n_dF elements
 * @return         SUPERMAG_OK on success, error code otherwise
 */
int supermag_proximity_critical_temp(
    double Tc0, double d_S, double xi_S, double xi_F, double E_ex,
    const double* d_F_arr, int n_dF, double* Tc_out
);

#ifdef __cplusplus
}
#endif

#endif
