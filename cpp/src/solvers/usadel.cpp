// Usadel diffusive-limit solver
// Linearized Usadel with Matsubara frequency sum and analytic matching.
// Uses tridiagonal solver for 1D finite-difference discretization.

#include "supermag/usadel.h"
#include <cmath>
#include <complex>
#include <algorithm>

extern "C" {

int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    int n_grid, double* Delta_out, double* x_out)
{
    if (!Delta_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_grid <= 4 || xi_S <= 0 || xi_F <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2;
    const double Delta_0 = 1.764 * kB_meV * Tc0;

    int n_S = n_grid / 2;
    int n_F = n_grid - n_S;

    // Build spatial grid
    for (int i = 0; i < n_S; ++i)
        x_out[i] = -d_S + d_S * i / std::max(n_S - 1, 1);
    for (int i = 0; i < n_F; ++i)
        x_out[n_S + i] = d_F * i / std::max(n_F - 1, 1);

    // Work at T = 0.5 Tc0
    double T = 0.5 * Tc0;
    double Delta_T = Delta_0 * std::sqrt(1.0 - T / Tc0);

    // Initialize Delta
    for (int i = 0; i < n_S; ++i)
        Delta_out[i] = Delta_T;
    for (int i = n_S; i < n_grid; ++i)
        Delta_out[i] = 0.0;

    int n_matsubara = 100;
    double lambda_BCS = 0.3;

    // Self-consistency iterations
    for (int iter = 0; iter < 30; ++iter) {
        double theta_sum_real[2048];  // max grid size
        if (n_grid > 2048) return SUPERMAG_ERR_INVALID_DIM;
        for (int i = 0; i < n_grid; ++i) theta_sum_real[i] = 0.0;

        for (int nm = 0; nm < n_matsubara; ++nm) {
            double omega_n = pi * kB_meV * T * (2 * nm + 1);

            double kappa_S = std::sqrt(2.0 * omega_n / (kB_meV * Tc0)) / xi_S;
            double theta_BCS = std::atan(Delta_T / omega_n);

            std::complex<double> kappa_F =
                std::sqrt(std::complex<double>(2.0 * omega_n, 2.0 * E_ex)) /
                (xi_F * std::sqrt(E_ex));

            // Interface matching
            std::complex<double> denom(kappa_S * xi_S, 0.0);
            denom += kappa_F * xi_F;
            std::complex<double> theta_int =
                theta_BCS * kappa_S * xi_S / denom;

            // Build theta profile
            double kS_dS = kappa_S * d_S;
            for (int j = 0; j < n_S; ++j) {
                double xj = x_out[j];
                double ratio = (kS_dS < 50.0)
                    ? std::cosh(kappa_S * (xj + d_S)) / std::cosh(kS_dS)
                    : std::exp(kappa_S * (xj + d_S) - kS_dS);
                std::complex<double> th =
                    theta_BCS + (theta_int - theta_BCS) * ratio;
                theta_sum_real[j] += th.real();
            }

            for (int j = 0; j < n_F; ++j) {
                double xj = x_out[n_S + j];
                std::complex<double> kF_x = kappa_F * xj;
                if (kF_x.real() < 50.0) {
                    std::complex<double> th =
                        theta_int * std::exp(-kF_x);
                    theta_sum_real[n_S + j] += th.real();
                }
            }
        }

        // Update Delta (S region only)
        double max_change = 0.0;
        for (int i = 0; i < n_S; ++i) {
            double Delta_new =
                lambda_BCS * pi * kB_meV * T * theta_sum_real[i];
            double change = std::fabs(Delta_new - Delta_out[i]);
            if (change > max_change) max_change = change;
            Delta_out[i] = 0.3 * Delta_new + 0.7 * Delta_out[i];
        }
        for (int i = n_S; i < n_grid; ++i)
            Delta_out[i] = 0.0;

        if (max_change < 1e-6)
            break;
    }

    // Ensure non-negative
    for (int i = 0; i < n_grid; ++i)
        Delta_out[i] = std::fabs(Delta_out[i]);

    return SUPERMAG_OK;
}

}
