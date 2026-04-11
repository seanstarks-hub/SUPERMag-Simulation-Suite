// Josephson current-phase relation solver
// Computes CPR for S/F/S junctions using Buzdin formula.
// Detects 0-pi transitions via sign change of critical current.

#include "supermag/josephson.h"
#include <cmath>

extern "C" {

int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T,
    int n_phases, double* phase_arr, double* current_out)
{
    if (!phase_arr || !current_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_phases <= 0 || xi_F <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2; // meV/K

    // Phase array
    for (int i = 0; i < n_phases; ++i)
        phase_arr[i] = 2.0 * pi * i / n_phases;

    // Complex wave vector: q = (1+i)/xi_F
    // exp(-q*d_F) = exp(-d_F/xi_F) * [cos(d_F/xi_F) - i*sin(d_F/xi_F)]
    double ratio = d_F / xi_F;
    double decay = std::exp(-ratio);
    // I_c ∝ exp(-d_F/xi_F) * cos(d_F/xi_F - π/4) [Buzdin]
    double I_c = decay * std::cos(ratio - pi / 4.0);

    // Temperature factor from BCS gap  [KNOWN-LIMIT-2: hardcoded Tc_ref]
    double Tc_ref = 9.2;
    double Delta_0 = 1.764 * kB_meV * Tc_ref;
    double t_ratio = T / Tc_ref;
    if (t_ratio > 0.9999) t_ratio = 0.9999;
    double Delta_T = Delta_0 * std::sqrt(1.0 - t_ratio);
    double T_factor = std::tanh(Delta_T / (2.0 * kB_meV * std::fmax(T, 0.01)));

    double amplitude = I_c * T_factor;

    // Compute CPR and find max for normalization
    double max_I = 0.0;
    for (int i = 0; i < n_phases; ++i) {
        current_out[i] = amplitude * std::sin(phase_arr[i]);
        if (std::fabs(current_out[i]) > max_I)
            max_I = std::fabs(current_out[i]);
    }

    // Normalize
    if (max_I > 0) {
        for (int i = 0; i < n_phases; ++i)
            current_out[i] /= max_I;
    }

    return SUPERMAG_OK;
}

}
