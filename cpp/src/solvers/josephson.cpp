// Josephson current-phase relation solver
// Computes CPR for S/F/S junctions via Matsubara frequency sum.
// Produces higher harmonics and detects 0-pi transitions.
// Extended version supports interface barrier (gamma_B) and absolute Ic.

#include "supermag/josephson.h"
#include "supermag/solver_options.h"
#include <cmath>
#include <complex>
#include <vector>

// Internal implementation
static int josephson_cpr_impl(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
    double gamma_B,
    int n_phases, const double* phase_arr_in,
    const supermag_solver_options_t *opts,
    double* current_out,
    double* Ic_out)
{
    if (!current_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_phases <= 0 || xi_F <= 0)
        return SUPERMAG_ERR_INVALID_DIM;
    if (Tc0 <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2; // meV/K

    // Phase array: use user-supplied or generate uniform grid
    std::vector<double> phase_buf;
    const double* phase_arr = phase_arr_in;
    if (!phase_arr) {
        phase_buf.resize(n_phases);
        for (int i = 0; i < n_phases; ++i)
            phase_buf[i] = 2.0 * pi * i / n_phases;
        phase_arr = phase_buf.data();
    }

    // BCS gap Δ(T) = 1.764·kB·Tc0·√(1−T/Tc0)
    double Delta_0 = 1.764 * kB_meV * Tc0;
    double t_ratio = T / Tc0;
    if (t_ratio > 0.9999) t_ratio = 0.9999;
    double Delta = Delta_0 * std::sqrt(1.0 - t_ratio);
    double Delta2 = Delta * Delta;

    // Temperature in meV for Matsubara sum
    double T_meV = kB_meV * std::fmax(T, 0.01);

    // Barrier damping factor  [2F-1]
    double barrier_damp = 1.0 / (1.0 + gamma_B);

    // Matsubara frequency sum for CPR  [EQ-9]
    supermag_solver_options_t defaults = supermag_default_solver_options();
    const supermag_solver_options_t *o = opts ? opts : &defaults;
    const int N_MAX = o->matsubara_max;
    double omega_cut = o->omega_cut_factor * Delta;
    std::complex<double> j_unit(0.0, 1.0);
    std::complex<double> phase_rot = std::exp(j_unit * pi / 4.0);

    for (int i = 0; i < n_phases; ++i)
        current_out[i] = 0.0;

    for (int n = 0; n <= N_MAX; ++n) {
        double omega_n = pi * T_meV * (2.0 * n + 1.0);
        if (omega_n > omega_cut && n > 0) break;

        double omega_n2 = omega_n * omega_n;

        // Complex wave vector in F layer
        std::complex<double> z(2.0 * omega_n / E_ex, 2.0);
        std::complex<double> q_n = std::sqrt(z) / xi_F;

        // F-layer propagator with barrier  [2F-1]
        std::complex<double> P_n = std::exp(-q_n * d_F) * phase_rot * barrier_damp;

        double denom_const = std::sqrt(omega_n2 + Delta2);

        for (int i = 0; i < n_phases; ++i) {
            double phi = phase_arr[i];
            double sin_phi = std::sin(phi);
            double sin_half = std::sin(phi / 2.0);

            double denom_phi = std::sqrt(omega_n2 + Delta2 * sin_half * sin_half);
            std::complex<double> numerator = P_n * Delta2 * sin_phi;

            current_out[i] += (numerator / (denom_phi * denom_const)).real();
        }
    }

    // Multiply by T_meV (prefactor of Matsubara sum)
    for (int i = 0; i < n_phases; ++i)
        current_out[i] *= T_meV;

    // Compute absolute critical current  [2F-2]
    double max_I = 0.0;
    for (int i = 0; i < n_phases; ++i) {
        if (std::fabs(current_out[i]) > max_I)
            max_I = std::fabs(current_out[i]);
    }

    if (Ic_out)
        *Ic_out = max_I;

    // Normalize to max |I| = 1
    if (max_I > 0) {
        for (int i = 0; i < n_phases; ++i)
            current_out[i] /= max_I;
    }

    return SUPERMAG_OK;
}

extern "C" {

int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
    double gamma_B,
    int n_phases, const double* phase_arr,
    const supermag_solver_options_t *opts,
    double* current_out,
    double* Ic_out)
{
    return josephson_cpr_impl(d_F, xi_F, E_ex, T, Tc0, gamma_B,
                               n_phases, phase_arr, opts,
                               current_out, Ic_out);
}

}
