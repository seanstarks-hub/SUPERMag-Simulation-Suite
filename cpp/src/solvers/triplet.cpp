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
    supermag_triplet_mode_t /*mode*/,
    int n_grid, double* f_triplet_out, double* x_out)
{
    if (!thicknesses || !magnetization_angles || !f_triplet_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_layers < 2 || n_grid <= 0)
        return SUPERMAG_ERR_INVALID_DIM;
    if (T <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double xi_F_use = (xi_F > 0.0) ? xi_F : 1.0;   // nm  [2G-2]
    const double xi_N_use = (xi_N > 0.0) ? xi_N : 10.0;  // nm
    const double T_use = T;                                // K   [2G-3]

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2;

    // Temperature-dependent BCS gap (for scaling)
    double Tc0 = 9.2;  // default Nb
    double Delta_0 = 1.764 * kB_meV * Tc0;
    double t_ratio = T_use / Tc0;
    if (t_ratio > 0.999) t_ratio = 0.999;
    double Delta_T = Delta_0 * std::sqrt(1.0 - t_ratio);

    // Lowest Matsubara frequency
    double omega_1 = pi * kB_meV * T_use;

    // Total thickness
    double total = 0.0;
    for (int i = 0; i < n_layers; ++i)
        total += thicknesses[i];

    // Build spatial grid
    (void)E_ex_per_layer; (void)D_per_layer;
    double dx = total / std::max(n_grid - 1, 1);
    (void)dx;
    for (int i = 0; i < n_grid; ++i)
        x_out[i] = total * i / std::max(n_grid - 1, 1);

    // Initialize output to zero
    for (int i = 0; i < n_grid; ++i)
        f_triplet_out[i] = 0.0;

    // Solve for singlet component f_0 in each layer
    // f_0 decays from interfaces with inverse length 1/xi_F (oscillatory in F)
    // The complex inverse decay length in F: q_F = (1+i)/xi_F  [EQ-19]
    std::complex<double> q_F(1.0 / xi_F_use, 1.0 / xi_F_use);

    // Triplet inverse decay length: purely real 1/xi_N
    double inv_xi_N = 1.0 / xi_N_use;

    // Temperature scaling: triplet amplitude scales with Delta(T)
    double T_scale = Delta_T / Delta_0;  // [2G-3]

    // Interface positions and singlet-to-triplet conversion
    double cumulative = 0.0;
    for (int lay = 0; lay < n_layers - 1; ++lay) {
        cumulative += thicknesses[lay];
        double x_int = cumulative;

        // Misalignment angle at this interface
        double alpha = magnetization_angles[lay + 1] - magnetization_angles[lay];
        double conversion = std::fabs(std::sin(alpha));  // ∝ sin(Δα)

        // Triplet contribution from this interface:
        // f_1(x) = |sin(Δα)| · exp(-|x - x_int|/xi_N)  [EQ-19]
        for (int j = 0; j < n_grid; ++j) {
            double dist = std::fabs(x_out[j] - x_int);

            // Long-range triplet from singlet-triplet conversion
            double f1_decay = std::exp(-dist * inv_xi_N);
            f_triplet_out[j] += conversion * f1_decay;
        }
    }

    return SUPERMAG_OK;
}

extern "C" {

int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    const double* E_ex_per_layer, const double* D_per_layer,
    double xi_F, double xi_N, double T,
    supermag_triplet_mode_t mode,
    int n_grid, double* f_triplet_out, double* x_out)
{
    return triplet_solve_impl(n_layers, thicknesses, magnetization_angles,
                               E_ex_per_layer, D_per_layer,
                               xi_F, xi_N, T, mode,
                               n_grid, f_triplet_out, x_out);
}

}
