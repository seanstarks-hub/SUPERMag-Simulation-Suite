#include "supermag/proximity.h"
#include <cmath>

extern "C" {

int supermag_proximity_critical_temp(
    double Tc0, double d_S, double xi_S, double xi_F, double E_ex,
    const double* d_F_arr, int n_dF, double* Tc_out)
{
    if (!d_F_arr || !Tc_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_dF < 1 || d_S <= 0.0 || xi_S <= 0.0 || xi_F <= 0.0 || E_ex <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Single-mode Buzdin approximation for Tc suppression in S/F bilayer.
    // The pair-breaking parameter includes oscillatory contribution from
    // FFLO-like pair amplitude in the F layer.
    //
    // Tc(d_F) = Tc0 * (1 - eta(d_F))
    // where eta captures the pair-breaking effect of the ferromagnet.
    //
    // eta ~ (xi_S / d_S) * Re[ (1 - exp(-2*kappa_F*d_F)) / (kappa_F * xi_S) ]
    // kappa_F = (1 + i) / xi_F  (complex inverse decay length in ferromagnet)
    //
    // This gives oscillatory Tc(d_F) with period ~ pi * xi_F and
    // exponential envelope ~ exp(-d_F / xi_F).

    double inv_xi_F = 1.0 / xi_F;

    for (int i = 0; i < n_dF; ++i) {
        double d_F = d_F_arr[i];
        if (d_F < 0.0) {
            Tc_out[i] = Tc0;
            continue;
        }

        // Complex kappa_F = (1+i)/xi_F, so kappa_F * d_F = (1+i)*d_F/xi_F
        double arg = d_F * inv_xi_F;

        // exp(-2*kappa_F*d_F) = exp(-2*arg) * (cos(2*arg) - i*sin(2*arg))
        double exp_decay = std::exp(-2.0 * arg);
        double real_exp = exp_decay * std::cos(2.0 * arg);
        double imag_exp = -exp_decay * std::sin(2.0 * arg);

        // (1 - exp(-2*kappa_F*d_F))
        double real_num = 1.0 - real_exp;
        double imag_num = -imag_exp;

        // kappa_F * xi_S = (1+i) * xi_S / xi_F
        double kappa_xi_S_real = xi_S / xi_F;
        double kappa_xi_S_imag = xi_S / xi_F;

        // Divide (real_num + i*imag_num) by (kappa_xi_S_real + i*kappa_xi_S_imag)
        double denom = kappa_xi_S_real * kappa_xi_S_real + kappa_xi_S_imag * kappa_xi_S_imag;
        double ratio_real = (real_num * kappa_xi_S_real + imag_num * kappa_xi_S_imag) / denom;

        // Pair-breaking parameter
        double eta = (xi_S / d_S) * ratio_real;

        // Clamp: Tc cannot go below 0 or above Tc0
        double Tc = Tc0 * (1.0 - eta);
        if (Tc < 0.0) Tc = 0.0;
        if (Tc > Tc0) Tc = Tc0;

        Tc_out[i] = Tc;
    }

    return SUPERMAG_OK;
}

}
