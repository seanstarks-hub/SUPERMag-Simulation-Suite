// Eilenberger clean-limit solver
// Riccati parameterization for quasiclassical Green's function.
// RK4 integration with angular (Fermi surface) averaging and
// Matsubara frequency summation.
//
// Riccati ODE:  [EQ-17]
//   da/dx = [ -2(ω_n + iE_ex)·a + Δ(1 - a²) ] / v_x
// where a is the Riccati amplitude, v_x = ℏv_F·cos(θ).
//
// Anomalous function:  f_n = 2a/(1 + |a|²)
// Total pair amplitude: f(x) = 2πT Σ_n Σ_θ w_θ |f_n|

#include "supermag/eilenberger.h"
#include "supermag/solver_options.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>

// RK4 right-hand side for the Riccati ODE
static inline std::complex<double> riccati_rhs(
    std::complex<double> a, double Delta, double omega_n, double E_ex, double vx)
{
    std::complex<double> omega_eff(omega_n, E_ex);
    return (-2.0 * omega_eff * a + Delta * (1.0 - a * a)) / vx;
}

extern "C" {

int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    double T,
    const supermag_solver_options_t *opts,
    int n_grid, double* f_out, double* x_out)
{
    if (!f_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_grid <= 4 || n_grid > 100000 || xi_S <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    const double pi = 3.14159265358979323846;
    const double kB_meV = 8.617333262e-2;
    const double Delta_0 = 1.764 * kB_meV * Tc0;
    // T must be positive (no default)
    if (T <= 0.0) return SUPERMAG_ERR_INVALID_DIM;
    double T_use = T;
    if (T_use >= Tc0) T_use = 0.999 * Tc0;
    double Delta_T = Delta_0 * std::sqrt(1.0 - T_use / Tc0);

    // hbar*v_F = 2*pi*kB*Tc * xi_S
    double hbar_vF = 2.0 * pi * kB_meV * Tc0 * xi_S;

    // Proportional grid splitting  [2C-4]
    int n_S = std::max(2, static_cast<int>(n_grid * d_S / (d_S + d_F)));
    int n_F = n_grid - n_S;
    if (n_F < 2) { n_F = 2; n_S = n_grid - 2; }

    // Spatial grid
    for (int i = 0; i < n_S; ++i)
        x_out[i] = -d_S + d_S * i / std::max(n_S - 1, 1);
    for (int i = 0; i < n_F; ++i)
        x_out[n_S + i] = d_F * i / std::max(n_F - 1, 1);

    double dx = (d_S + d_F) / std::max(n_grid - 1, 1);

    // Initialize f_out to zero
    for (int i = 0; i < n_grid; ++i)
        f_out[i] = 0.0;

    // Angular quadrature: Gauss-Legendre 16-point
    const int n_angles = 16;
    static const double gl_nodes[16] = {
        -0.9894009349, -0.9445750230, -0.8656312024, -0.7554044084,
        -0.6178762444, -0.4580167776, -0.2816035508, -0.0950125098,
         0.0950125098,  0.2816035508,  0.4580167776,  0.6178762444,
         0.7554044084,  0.8656312024,  0.9445750230,  0.9894009349
    };
    static const double gl_weights[16] = {
        0.0271524594, 0.0622535239, 0.0951585117, 0.1246289713,
        0.1495959888, 0.1691565194, 0.1826034150, 0.1894506105,
        0.1894506105, 0.1826034150, 0.1691565194, 0.1495959888,
        0.1246289713, 0.0951585117, 0.0622535239, 0.0271524594
    };

    double weight_sum = 0.0;
    for (int k = 0; k < n_angles; ++k)
        weight_sum += gl_weights[k];

    // Matsubara frequency sum  [2C-2]
    double T_meV = kB_meV * T_use;
    supermag_solver_options_t defaults = supermag_default_solver_options();
    const supermag_solver_options_t *o = opts ? opts : &defaults;
    const int N_FREQ_MAX = o->matsubara_max;
    double omega_cut = o->omega_cut_factor * std::max(Delta_T, 1e-6);

    // Accumulate contributions from all frequencies
    std::vector<double> f_accum(n_grid, 0.0);
    int n_freq_used = 0;

    for (int nf = 0; nf < N_FREQ_MAX; ++nf) {
        double omega_n = pi * T_meV * (2.0 * nf + 1.0);
        if (omega_n > omega_cut && nf > 0) break;
        ++n_freq_used;

        // BCS Riccati parameter at this frequency
        double a_BCS = Delta_T / (omega_n + std::sqrt(omega_n * omega_n + Delta_T * Delta_T));

        for (int k = 0; k < n_angles; ++k) {
            double cos_theta = gl_nodes[k];
            if (std::fabs(cos_theta) < 1e-10) continue;

            double vx = hbar_vF * cos_theta;

            // Riccati integration with RK4  [2C-1]
            std::vector<std::complex<double>> a(n_grid);

            if (cos_theta > 0) {
                // Right-mover: left → right
                a[0] = std::complex<double>(a_BCS, 0.0);
                for (int j = 0; j < n_grid - 1; ++j) {
                    double xj = x_out[j];
                    double D_val = (xj < 0) ? Delta_T : 0.0;
                    double E_val = (xj < 0) ? 0.0 : E_ex;

                    // RK4 step
                    auto k1 = riccati_rhs(a[j], D_val, omega_n, E_val, vx);
                    auto k2 = riccati_rhs(a[j] + 0.5*dx*k1, D_val, omega_n, E_val, vx);
                    auto k3 = riccati_rhs(a[j] + 0.5*dx*k2, D_val, omega_n, E_val, vx);
                    auto k4 = riccati_rhs(a[j] + dx*k3, D_val, omega_n, E_val, vx);
                    a[j + 1] = a[j] + (dx / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

                    // Physical normalization constraint  [2C-3]
                    if (std::abs(a[j + 1]) > 1.0)
                        a[j + 1] /= std::abs(a[j + 1]);
                }
            } else {
                // Left-mover: right → left
                a[n_grid - 1] = std::complex<double>(0.0, 0.0);
                for (int j = n_grid - 1; j > 0; --j) {
                    double xj = x_out[j];
                    double D_val = (xj < 0) ? Delta_T : 0.0;
                    double E_val = (xj < 0) ? 0.0 : E_ex;

                    // RK4 step (backward: dx → -dx)
                    auto k1 = riccati_rhs(a[j], D_val, omega_n, E_val, vx);
                    auto k2 = riccati_rhs(a[j] - 0.5*dx*k1, D_val, omega_n, E_val, vx);
                    auto k3 = riccati_rhs(a[j] - 0.5*dx*k2, D_val, omega_n, E_val, vx);
                    auto k4 = riccati_rhs(a[j] - dx*k3, D_val, omega_n, E_val, vx);
                    a[j - 1] = a[j] - (dx / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

                    // Physical normalization constraint  [2C-3]
                    if (std::abs(a[j - 1]) > 1.0)
                        a[j - 1] /= std::abs(a[j - 1]);
                }
            }

            // f_n = 2a / (1 + |a|²)
            for (int j = 0; j < n_grid; ++j) {
                double abs_a2 = std::norm(a[j]);
                std::complex<double> f_j = 2.0 * a[j] / (1.0 + abs_a2);
                f_accum[j] += gl_weights[k] * std::abs(f_j);
            }
        }
    }

    // Normalize by total angular weight and number of frequencies
    double norm = weight_sum * std::max(n_freq_used, 1);
    for (int j = 0; j < n_grid; ++j)
        f_out[j] = f_accum[j] / norm;

    return SUPERMAG_OK;
}

}
