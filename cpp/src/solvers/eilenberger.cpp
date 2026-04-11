// Eilenberger clean-limit solver
// Riccati parameterization for quasiclassical Green's function.
// Forward Euler integration with angular (Fermi surface) averaging.

#include "supermag/eilenberger.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>

extern "C" {

int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    int n_grid, double* f_out, double* x_out)
{
    if (!f_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_grid <= 4 || n_grid > 100000 || xi_S <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2;
    const double Delta_0 = 1.764 * kB_meV * Tc0;
    // Use T = 0.5 * Tc0  [KNOWN-LIMIT-1]
    double T_use = 0.5 * Tc0;
    double Delta_T = Delta_0 * std::sqrt(1.0 - T_use / Tc0);

    // hbar*v_F = 2*pi*kB*Tc * xi_S
    double hbar_vF = 2.0 * pi * kB_meV * Tc0 * xi_S;

    // Lowest Matsubara frequency
    double omega_1 = pi * kB_meV * T_use;

    // BCS Riccati parameter
    double a_BCS = Delta_T / (omega_1 + std::sqrt(omega_1 * omega_1 + Delta_T * Delta_T));

    int n_S = n_grid / 2;
    int n_F = n_grid - n_S;

    // Spatial grid
    for (int i = 0; i < n_S; ++i)
        x_out[i] = -d_S + d_S * i / std::max(n_S - 1, 1);
    for (int i = 0; i < n_F; ++i)
        x_out[n_S + i] = d_F * i / std::max(n_F - 1, 1);

    double dx = (d_S + d_F) / std::max(n_grid - 1, 1);

    // Initialize f_out to zero
    for (int i = 0; i < n_grid; ++i)
        f_out[i] = 0.0;

    // Angular quadrature: Gauss-Legendre 8-point
    const int n_angles = 8;
    static const double gl_nodes[8] = {
        -0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425,
         0.1834346425,  0.5255324099,  0.7966664774,  0.9602898565
    };
    static const double gl_weights[8] = {
        0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
        0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363
    };

    double weight_sum = 0.0;
    for (int k = 0; k < n_angles; ++k)
        weight_sum += gl_weights[k];

    for (int k = 0; k < n_angles; ++k) {
        double cos_theta = gl_nodes[k];
        if (std::fabs(cos_theta) < 1e-10) continue;

        double vx = hbar_vF * cos_theta;

        // Riccati integration
        std::vector<std::complex<double>> a(n_grid);

        if (cos_theta > 0) {
            // Right-mover: left → right
            a[0] = std::complex<double>(a_BCS, 0.0);
            for (int j = 0; j < n_grid - 1; ++j) {
                double xj = x_out[j];
                double D = (xj < 0) ? Delta_T : 0.0;
                double E = (xj < 0) ? 0.0 : E_ex;
                std::complex<double> da =
                    (-2.0 * std::complex<double>(omega_1, E) * a[j]
                     + D * (1.0 - a[j] * a[j])) / vx;
                a[j + 1] = a[j] + da * dx;
                if (std::abs(a[j + 1]) > 2.0)
                    a[j + 1] *= 2.0 / std::abs(a[j + 1]);
            }
        } else {
            // Left-mover: right → left
            a[n_grid - 1] = std::complex<double>(0.0, 0.0);
            for (int j = n_grid - 1; j > 0; --j) {
                double xj = x_out[j];
                double D = (xj < 0) ? Delta_T : 0.0;
                double E = (xj < 0) ? 0.0 : E_ex;
                std::complex<double> da =
                    (-2.0 * std::complex<double>(omega_1, E) * a[j]
                     + D * (1.0 - a[j] * a[j])) / vx;
                a[j - 1] = a[j] - da * dx;
                if (std::abs(a[j - 1]) > 2.0)
                    a[j - 1] *= 2.0 / std::abs(a[j - 1]);
            }
        }

        // f = 2a / (1 + |a|²)
        for (int j = 0; j < n_grid; ++j) {
            double abs_a2 = std::norm(a[j]);
            std::complex<double> f_j = 2.0 * a[j] / (1.0 + abs_a2);
            f_out[j] += gl_weights[k] * std::abs(f_j);
        }
    }

    // Normalize by total weight
    for (int j = 0; j < n_grid; ++j)
        f_out[j] /= weight_sum;

    return SUPERMAG_OK;
}

}
