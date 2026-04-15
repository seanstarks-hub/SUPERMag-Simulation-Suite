// Spin-triplet superconductivity solver
// Linearized Usadel equations for singlet-triplet conversion at
// non-collinear magnetic interfaces.
//
// Coupled equations for singlet (f_0) and long-range triplet (f_1):  [EQ-19]
//   D_j · d²f_0/dx² - 2ω_n·f_0 + 2iE_ex·f_0 = 0     (singlet in F)
//   D_j · d²f_1/dx² - 2ω_n·f_1 = 0                    (triplet, no exchange)
//
// At non-collinear interfaces (misalignment angle Δα):
//   f_1(x_int) += sin(Δα) · f_0(x_int)    (singlet-to-triplet conversion)
//
// The triplet component f_1 decays with xi_N (long-range) while f_0
// decays with xi_F (short-range, oscillatory).

#include "supermag/triplet.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>

// Internal implementation
static int triplet_solve_impl(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    const double* E_ex_per_layer, const double* D_per_layer,
    double xi_F, double xi_N, double T,
    double Tc0,
    supermag_triplet_mode_t /*mode*/,
    int n_grid, double* f_triplet_out, double* x_out)
{
    if (!thicknesses || !magnetization_angles || !f_triplet_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_layers < 2 || n_grid <= 0)
        return SUPERMAG_ERR_INVALID_DIM;
    if (T <= 0.0 || Tc0 <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double xi_F_use = (xi_F > 0.0) ? xi_F : 1.0;   // nm  [2G-2]
    const double xi_N_use = (xi_N > 0.0) ? xi_N : 10.0;  // nm
    const double T_use = T;                                // K   [2G-3]

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2;

    // Temperature-dependent BCS gap (for scaling)
    double Delta_0 = 1.764 * kB_meV * Tc0;
    double t_ratio = T_use / Tc0;
    if (t_ratio > 0.999) t_ratio = 0.999;
    double Delta_T = Delta_0 * std::sqrt(1.0 - t_ratio);

    // Total thickness
    double total = 0.0;
    for (int i = 0; i < n_layers; ++i)
        total += thicknesses[i];

    // Build spatial grid
    for (int i = 0; i < n_grid; ++i)
        x_out[i] = total * i / std::max(n_grid - 1, 1);

    // Initialize output to zero
    for (int i = 0; i < n_grid; ++i)
        f_triplet_out[i] = 0.0;

    // Per-layer coherence lengths for singlet decay  [EQ-19]
    // xi_F_lay = sqrt(D_lay / E_ex_lay)  when both arrays given
    // xi_F_lay = xi_F * sqrt(E_ex[0] / E_ex[lay])  when only E_ex given
    // xi_F_lay = xi_F (global)  when neither given
    std::vector<double> xi_F_layer(n_layers, xi_F_use);
    if (E_ex_per_layer && D_per_layer) {
        for (int lay = 0; lay < n_layers; ++lay) {
            double E = std::fabs(E_ex_per_layer[lay]);
            double D = std::fabs(D_per_layer[lay]);
            if (E > 1e-30 && D > 1e-30)
                xi_F_layer[lay] = std::sqrt(D / E);
            else
                xi_F_layer[lay] = xi_F_use;
        }
    } else if (E_ex_per_layer) {
        // Scale global xi_F by exchange energy ratio (xi_F ∝ 1/√E_ex at fixed D)
        double E_ref = std::fabs(E_ex_per_layer[0]);
        if (E_ref > 1e-30) {
            for (int lay = 0; lay < n_layers; ++lay) {
                double E = std::fabs(E_ex_per_layer[lay]);
                if (E > 1e-30)
                    xi_F_layer[lay] = xi_F_use * std::sqrt(E_ref / E);
                else
                    xi_F_layer[lay] = xi_F_use;
            }
        }
    }

    // Triplet inverse decay length: purely real 1/xi_N
    double inv_xi_N = 1.0 / xi_N_use;

    // Temperature scaling: triplet amplitude scales with Delta(T)
    double T_scale = Delta_T / Delta_0;  // [2G-3]

    // Singlet amplitude at each interface  [EQ-19]
    // Singlet decays through each layer with per-layer q_F.
    // The singlet is injected from the S layer and decays through each F layer.
    // At interface between layers lay and lay+1:
    //   f_0(x_int) = T_scale * exp(-Σ d_j / xi_F_layer[j])  for j = 0..lay
    std::vector<double> singlet_at_interface(n_layers - 1);
    {
        double cumul_decay = 0.0;
        for (int lay = 0; lay < n_layers - 1; ++lay) {
            cumul_decay += thicknesses[lay] / xi_F_layer[lay];
            singlet_at_interface[lay] = T_scale * std::exp(-cumul_decay);
        }
    }

    // Interface positions and singlet-to-triplet conversion
    double cumulative = 0.0;
    for (int lay = 0; lay < n_layers - 1; ++lay) {
        cumulative += thicknesses[lay];
        double x_int = cumulative;

        // Misalignment angle at this interface
        double alpha = magnetization_angles[lay + 1] - magnetization_angles[lay];
        double conversion = std::fabs(std::sin(alpha));  // ∝ sin(Δα)

        // Modulate by singlet amplitude at this interface
        double source = conversion * singlet_at_interface[lay];

        // Triplet contribution from this interface:
        // f_1(x) = source · exp(-|x - x_int|/xi_N)  [EQ-19]
        for (int j = 0; j < n_grid; ++j) {
            double dist = std::fabs(x_out[j] - x_int);

            // Long-range triplet from singlet-triplet conversion
            double f1_decay = std::exp(-dist * inv_xi_N);
            f_triplet_out[j] += source * f1_decay;
        }
    }

    return SUPERMAG_OK;
}

extern "C" {

int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    const double* E_ex_per_layer, const double* D_per_layer,
    double xi_F, double xi_N, double T,
    double Tc0,
    supermag_triplet_mode_t mode,
    int n_grid, double* f_triplet_out, double* x_out)
{
    return triplet_solve_impl(n_layers, thicknesses, magnetization_angles,
                               E_ex_per_layer, D_per_layer,
                               xi_F, xi_N, T, Tc0, mode,
                               n_grid, f_triplet_out, x_out);
}

}
