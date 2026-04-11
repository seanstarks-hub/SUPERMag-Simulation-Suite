// Ginzburg-Landau free energy minimization
// TDGL relaxation on 2D grid with periodic boundary conditions.

#include "supermag/ginzburg_landau.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>

extern "C" {

int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    double* psi_real, double* psi_imag)
{
    if (!psi_real || !psi_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (nx <= 0 || ny <= 0 || dx <= 0 || beta <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    int N = nx * ny;

    // Equilibrium |ψ|² = -α/β for uniform state (α < 0)
    double psi_eq = (alpha < 0) ? std::sqrt(-alpha / beta) : 0.0;

    // Initialize with bulk value + small perturbation
    // Use deterministic "random" seed for reproducibility
    for (int i = 0; i < N; ++i) {
        // Simple deterministic perturbation
        double noise_r = 0.05 * std::sin(i * 1.2345 + 0.1);
        double noise_i = 0.05 * std::cos(i * 2.3456 + 0.2);
        psi_real[i] = psi_eq * (1.0 + noise_r);
        psi_imag[i] = psi_eq * noise_i;
    }

    double xi2 = dx * dx;
    double dt = 0.1 * dx * dx / (4.0 * xi2 + 1e-15);
    int n_steps = 2000;
    double inv_dx2 = 1.0 / (dx * dx);

    // Double-buffer: compute Laplacian from snapshot, then update  [EQ-11]
    std::vector<double> buf_r(N), buf_i(N);

    for (int step = 0; step < n_steps; ++step) {
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

                // Laplacian from snapshot (consistent Euler step)
                double lap_r = (buf_r[left] + buf_r[right]
                              + buf_r[down] + buf_r[up]
                              - 4.0 * buf_r[idx]) * inv_dx2;
                double lap_i = (buf_i[left] + buf_i[right]
                              + buf_i[down] + buf_i[up]
                              - 4.0 * buf_i[idx]) * inv_dx2;

                double abs2 = buf_r[idx] * buf_r[idx]
                            + buf_i[idx] * buf_i[idx];

                // dψ/dt = -αψ - β|ψ|²ψ + ξ²∇²ψ
                double dpsi_r = -alpha * buf_r[idx]
                              - beta * abs2 * buf_r[idx]
                              + xi2 * lap_r;
                double dpsi_i = -alpha * buf_i[idx]
                              - beta * abs2 * buf_i[idx]
                              + xi2 * lap_i;

                psi_real[idx] = buf_r[idx] + dt * dpsi_r;
                psi_imag[idx] = buf_i[idx] + dt * dpsi_i;

                double res = std::sqrt(dpsi_r * dpsi_r + dpsi_i * dpsi_i);
                if (res > max_residual) max_residual = res;
            }
        }

        if (step % 100 == 99 && max_residual < 1e-8)
            break;
    }

    return SUPERMAG_OK;
}

}
