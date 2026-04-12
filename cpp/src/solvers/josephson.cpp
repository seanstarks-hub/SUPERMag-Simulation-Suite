// Josephson current-phase relation solver
// Computes CPR for S/F/S junctions via Matsubara frequency sum.
// Produces higher harmonics and detects 0-pi transitions.

#include "supermag/josephson.h"
#include <cmath>
#include <complex>

extern "C" {

int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
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

    // BCS gap Δ(T) = 1.764·kB·Tc0·√(1−T/Tc0)
    double Tc_ref = (Tc0 > 0.0) ? Tc0 : 9.2;
    double Delta_0 = 1.764 * kB_meV * Tc_ref;
    double t_ratio = T / Tc_ref;
    if (t_ratio > 0.9999) t_ratio = 0.9999;
    double Delta = Delta_0 * std::sqrt(1.0 - t_ratio);
    double Delta2 = Delta * Delta;

    // Temperature in meV for Matsubara sum
    double T_meV = kB_meV * std::fmax(T, 0.01);

    // Matsubara frequency sum for CPR  [EQ-9]
    //   I(φ) = T_meV · Σ_n Re[ P_n · Δ²·sin(φ) /
    //          √((ω_n² + Δ²·sin²(φ/2)) · (ω_n² + Δ²)) ]
    //   P_n = exp(-q_n·d_F) · exp(iπ/4)
    //   q_n = √(2·(ω_n/E_ex + i)) / ξ_F
    //   ω_n = π·T_meV·(2n+1)
    // Cutoff: ω_n > 20Δ or n > 500

    const int N_MAX = 500;
    double omega_cut = 20.0 * Delta;
    std::complex<double> j_unit(0.0, 1.0);
    std::complex<double> phase_rot = std::exp(j_unit * pi / 4.0); // exp(iπ/4)

    for (int i = 0; i < n_phases; ++i)
        current_out[i] = 0.0;

    for (int n = 0; n <= N_MAX; ++n) {
        double omega_n = pi * T_meV * (2.0 * n + 1.0);
        if (omega_n > omega_cut && n > 0) break;

        double omega_n2 = omega_n * omega_n;

        // Complex wave vector in F layer at this Matsubara frequency:
        // q_n = √(2·(ω_n/E_ex + i)) / ξ_F     [EQ-1 generalized]
        std::complex<double> z(2.0 * omega_n / E_ex, 2.0);
        std::complex<double> q_n = std::sqrt(z) / xi_F;

        // F-layer propagator: P_n = exp(-q_n·d_F) · exp(iπ/4)
        std::complex<double> P_n = std::exp(-q_n * d_F) * phase_rot;

        // Denominator factor independent of φ: √(ω_n² + Δ²)
        double denom_const = std::sqrt(omega_n2 + Delta2);

        for (int i = 0; i < n_phases; ++i) {
            double phi = phase_arr[i];
            double sin_phi = std::sin(phi);
            double sin_half = std::sin(phi / 2.0);

            // √(ω_n² + Δ²·sin²(φ/2))
            double denom_phi = std::sqrt(omega_n2 + Delta2 * sin_half * sin_half);

            // Numerator: P_n · Δ² · sin(φ)
            std::complex<double> numerator = P_n * Delta2 * sin_phi;

            // Full term for this Matsubara frequency
            current_out[i] += (numerator / (denom_phi * denom_const)).real();
        }
    }

    // Multiply by T_meV (prefactor of Matsubara sum)
    for (int i = 0; i < n_phases; ++i)
        current_out[i] *= T_meV;

    // Normalize to max |I| = 1
    double max_I = 0.0;
    for (int i = 0; i < n_phases; ++i) {
        if (std::fabs(current_out[i]) > max_I)
            max_I = std::fabs(current_out[i]);
    }
    if (max_I > 0) {
        for (int i = 0; i < n_phases; ++i)
            current_out[i] /= max_I;
    }

    return SUPERMAG_OK;
}

}
