// Ginzburg-Landau free energy minimization
// TDGL relaxation on 2D grid with periodic boundary conditions.
// Supports gauge-covariant kinetic term with vector potential for applied field.
//
// GL equation (TDGL relaxation):  [EQ-11]
//   ∂ψ/∂t = -αψ - β|ψ|²ψ + ξ²(∇ - iA)²ψ
// where ξ² = 1/(2κ²) in GL units.
//
// Vector potential in Landau gauge: A = (0, H·x, 0)
// Gauge-covariant Laplacian:
//   (∇ - iA)²ψ = ∂²ψ/∂x² + (∂/∂y - iA_y)²ψ
//
// Self-consistent field equation (Maxwell):  [EQ-18]
//   J = (κ²/2) Im[ψ*(∇ - iA)ψ]
//   ∇²A = -J

#include "supermag/ginzburg_landau.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>

// Internal implementation with H_applied parameter
static int gl_minimize_impl(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    double H_applied,
    double* psi_real, double* psi_imag)
{
    if (!psi_real || !psi_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (nx <= 0 || ny <= 0 || dx <= 0 || beta <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    int N = nx * ny;

    // Equilibrium |ψ|² = -α/β for uniform state (α < 0)
    double psi_eq = (alpha < 0) ? std::sqrt(-alpha / beta) : 0.0;

    // Preserve user initial condition if psi_real[0] != 0  [2E-4]
    bool user_ic = (std::fabs(psi_real[0]) > 1e-30);
    if (!user_ic) {
        // Initialize with bulk value + small perturbation
        for (int i = 0; i < N; ++i) {
            double noise_r = 0.05 * std::sin(i * 1.2345 + 0.1);
            double noise_i = 0.05 * std::cos(i * 2.3456 + 0.2);
            psi_real[i] = psi_eq * (1.0 + noise_r);
            psi_imag[i] = psi_eq * noise_i;
        }
    }

    // Coherence length from κ: ξ² = 1/(2κ²)  [2E-1]
    double xi2 = 1.0 / (2.0 * kappa * kappa);
    double dt = 0.1 * dx * dx / (4.0 * xi2 + 1e-15);
    double inv_dx2 = 1.0 / (dx * dx);

    // Vector potential in Landau gauge: A_y(x) = H · x  [2E-2]
    // A_x = 0, A_y = H * x_position
    std::vector<double> Ay(N, 0.0);
    if (std::fabs(H_applied) > 1e-30) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                double x_pos = ix * dx;
                Ay[iy * nx + ix] = H_applied * x_pos;
            }
        }
    }

    // Self-consistent vector potential update  [2E-3]
    std::vector<double> Ay_update(N, 0.0);

    // Adaptive convergence  [2E-5]
    const int max_steps = 5000;
    const double conv_tol = 1e-8;
    int converge_count = 0;

    // Double-buffer
    std::vector<double> buf_r(N), buf_i(N);

    for (int step = 0; step < max_steps; ++step) {
        double max_residual = 0.0;

        // Copy current state into read-only snapshot
        for (int k = 0; k < N; ++k) {
            buf_r[k] = psi_real[k];
            buf_i[k] = psi_imag[k];
        }

        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                int idx = iy * nx + ix;

                // Periodic neighbors
                int left  = iy * nx + ((ix - 1 + nx) % nx);
                int right = iy * nx + ((ix + 1) % nx);
                int down  = ((iy - 1 + ny) % ny) * nx + ix;
                int up    = ((iy + 1) % ny) * nx + ix;

                double ay = Ay[idx];
                (void)ay;

                if (std::fabs(H_applied) < 1e-30) {
                    // No field: standard Laplacian
                    double lap_r = (buf_r[left] + buf_r[right]
                                  + buf_r[down] + buf_r[up]
                                  - 4.0 * buf_r[idx]) * inv_dx2;
                    double lap_i = (buf_i[left] + buf_i[right]
                                  + buf_i[down] + buf_i[up]
                                  - 4.0 * buf_i[idx]) * inv_dx2;

                    double abs2 = buf_r[idx] * buf_r[idx] + buf_i[idx] * buf_i[idx];

                    double dpsi_r = -alpha * buf_r[idx] - beta * abs2 * buf_r[idx] + xi2 * lap_r;
                    double dpsi_i = -alpha * buf_i[idx] - beta * abs2 * buf_i[idx] + xi2 * lap_i;

                    psi_real[idx] = buf_r[idx] + dt * dpsi_r;
                    psi_imag[idx] = buf_i[idx] + dt * dpsi_i;

                    double res = std::sqrt(dpsi_r * dpsi_r + dpsi_i * dpsi_i);
                    if (res > max_residual) max_residual = res;
                } else {
                    // Gauge-covariant Laplacian: (∇ - iA)²ψ  [2E-2]
                    // x-direction (A_x = 0): standard second derivative
                    double d2x_r = (buf_r[left] + buf_r[right] - 2.0 * buf_r[idx]) * inv_dx2;
                    double d2x_i = (buf_i[left] + buf_i[right] - 2.0 * buf_i[idx]) * inv_dx2;

                    // y-direction with A_y: (∂_y - iA_y)²ψ
                    // Use Peierls phases: ψ_up * exp(-iA_y·dx), ψ_down * exp(+iA_y·dx)
                    double phase_up = -ay * dx;
                    double phase_down = ay * dx;

                    // exp(iφ) = cos(φ) + i·sin(φ)
                    double c_up = std::cos(phase_up), s_up = std::sin(phase_up);
                    double c_dn = std::cos(phase_down), s_dn = std::sin(phase_down);

                    // Rotated neighbor values
                    double up_r = c_up * buf_r[up] - s_up * buf_i[up];
                    double up_i = s_up * buf_r[up] + c_up * buf_i[up];
                    double dn_r = c_dn * buf_r[down] - s_dn * buf_i[down];
                    double dn_i = s_dn * buf_r[down] + c_dn * buf_i[down];

                    double d2y_r = (up_r + dn_r - 2.0 * buf_r[idx]) * inv_dx2;
                    double d2y_i = (up_i + dn_i - 2.0 * buf_i[idx]) * inv_dx2;

                    double cov_lap_r = d2x_r + d2y_r;
                    double cov_lap_i = d2x_i + d2y_i;

                    double abs2 = buf_r[idx] * buf_r[idx] + buf_i[idx] * buf_i[idx];

                    double dpsi_r = -alpha * buf_r[idx] - beta * abs2 * buf_r[idx] + xi2 * cov_lap_r;
                    double dpsi_i = -alpha * buf_i[idx] - beta * abs2 * buf_i[idx] + xi2 * cov_lap_i;

                    psi_real[idx] = buf_r[idx] + dt * dpsi_r;
                    psi_imag[idx] = buf_i[idx] + dt * dpsi_i;

                    double res = std::sqrt(dpsi_r * dpsi_r + dpsi_i * dpsi_i);
                    if (res > max_residual) max_residual = res;
                }
            }
        }

        // Self-consistent A update every 50 steps when field is applied  [2E-3]
        if (std::fabs(H_applied) > 1e-30 && step % 50 == 49) {
            double kappa2_half = kappa * kappa * 0.5;
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    int idx = iy * nx + ix;
                    int up = ((iy + 1) % ny) * nx + ix;
                    int down = ((iy - 1 + ny) % ny) * nx + ix;

                    // Supercurrent J_y = Im[ψ* (∂_y - iA_y)ψ]
                    double dpsi_dy_r = (psi_real[up] - psi_real[down]) / (2.0 * dx);
                    double dpsi_dy_i = (psi_imag[up] - psi_imag[down]) / (2.0 * dx);

                    // (∂_y - iA_y)ψ = (dpsi_dy_r + Ay*ψ_i) + i(dpsi_dy_i - Ay*ψ_r)
                    double cov_r = dpsi_dy_r + Ay[idx] * psi_imag[idx];
                    double cov_i = dpsi_dy_i - Ay[idx] * psi_real[idx];

                    // J = Im[ψ* · cov_deriv] = ψ_r·cov_i - ψ_i·cov_r
                    double Jy = kappa2_half *
                        (psi_real[idx] * cov_i - psi_imag[idx] * cov_r);

                    // Small correction to A_y (relaxation toward ∇²A = -J)
                    Ay[idx] += 0.01 * dt * Jy;
                }
            }
        }

        // Convergence check  [2E-5]
        if (max_residual < conv_tol) {
            converge_count++;
            if (converge_count >= 3) break;  // converged for 3 consecutive checks
        } else {
            converge_count = 0;
        }
    }

    return SUPERMAG_OK;
}

extern "C" {

int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    supermag_gl_mode_t mode, double H_applied,
    double* psi_real, double* psi_imag)
{
    double H = (mode == SUPERMAG_GL_GAUGE) ? H_applied : 0.0;
    return gl_minimize_impl(alpha, beta, kappa, nx, ny, dx, H, psi_real, psi_imag);
}

}
