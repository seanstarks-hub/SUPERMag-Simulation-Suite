#ifndef SUPERMAG_OPTIMIZER_H
#define SUPERMAG_OPTIMIZER_H

#include "error.h"
#include "proximity.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Optimize d_F to match target Tc using golden-section search.
 * Varies d_F within [d_F_lo, d_F_hi] to minimize |Tc(d_F) - Tc_target|.
 *
 * @param params      Base parameter set (d_F field is overwritten)
 * @param depairing   Depairing channels (NULL if none)
 * @param d_F_lo      Lower bound for d_F (nm)
 * @param d_F_hi      Upper bound for d_F (nm)
 * @param Tc_target   Target critical temperature (K)
 * @param d_F_out     Output: optimal d_F value
 */
int supermag_optimize_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double d_F_lo, double d_F_hi,
    double Tc_target,
    double *d_F_out
);

/**
 * Inverse solve: given Tc_target, find d_F that produces it via Brent's method.
 *
 * @param params     Base parameter set (d_F field is overwritten)
 * @param depairing  Depairing channels (NULL if none)
 * @param Tc_target  Target Tc (K)
 * @param d_F_lo     Lower bound for d_F (nm)
 * @param d_F_hi     Upper bound for d_F (nm)
 * @param d_F_out    Output: d_F that gives Tc_target
 */
int supermag_inverse_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double Tc_target,
    double d_F_lo, double d_F_hi,
    double *d_F_out
);

/**
 * Nelder-Mead least-squares fit of proximity parameters to experimental data.
 * Optimizes selected parameters (flagged by fit_flags) to minimize chi^2
 * between computed Tc(d_F) and experimental Tc data.
 *
 * @param params        Base parameter set (fitted fields overwritten)
 * @param depairing     Depairing channels (NULL if none)
 * @param d_F_data      Experimental d_F values (nm), length n_data
 * @param Tc_data       Experimental Tc values (K), length n_data
 * @param n_data        Number of data points (>= 2)
 * @param fit_gamma     If nonzero, include gamma in fit
 * @param fit_gamma_B   If nonzero, include gamma_B in fit
 * @param fit_E_ex      If nonzero, include E_ex in fit
 * @param fit_xi_F      If nonzero, include xi_F in fit
 * @param chi2_out      Output: final chi^2 value (NULL to skip)
 */
int supermag_fit_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    const double *d_F_data, const double *Tc_data, int n_data,
    int fit_gamma, int fit_gamma_B, int fit_E_ex, int fit_xi_F,
    double *chi2_out
);

#ifdef __cplusplus
}
#endif

#endif
