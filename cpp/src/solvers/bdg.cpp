// BdG tight-binding Hamiltonian solver
// Constructs 2N×2N Nambu-doubled BdG matrix and diagonalizes via cyclic Jacobi.
// Optional eigenvector output for LDOS analysis.

#include "supermag/bdg.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

// Cyclic Jacobi eigenvalue algorithm for real symmetric matrices.  [2D-1]
// Systematically sweeps all off-diagonal elements in row-major order.
// Returns 0 if converged, -1 if max_iter reached.
// If V is non-null, accumulates the rotation matrix (eigenvectors).
static int cyclic_jacobi(double* A, int n, double* eigvals,
                          double* V, int max_iter = 200) {
    std::vector<double> M(n * n);
    std::memcpy(M.data(), A, n * n * sizeof(double));

    // Initialize eigenvector matrix to identity if requested  [2D-3]
    if (V) {
        std::memset(V, 0, n * n * sizeof(double));
        for (int i = 0; i < n; ++i)
            V[i * n + i] = 1.0;
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        // Check convergence: sum of squares of off-diagonal elements
        double off_norm = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                off_norm += M[i * n + j] * M[i * n + j];
        if (off_norm < 1e-24) {
            for (int i = 0; i < n; ++i)
                eigvals[i] = M[i * n + i];
            std::sort(eigvals, eigvals + n);
            return 0;
        }

        // Sweep all off-diagonal pairs systematically
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double apq = M[p * n + q];
                if (std::fabs(apq) < 1e-15) continue;

                double app = M[p * n + p], aqq = M[q * n + q];
                double theta = 0.5 * std::atan2(2.0 * apq, app - aqq);
                double c = std::cos(theta), s = std::sin(theta);

                // Apply Jacobi rotation to M
                for (int i = 0; i < n; ++i) {
                    double mip = M[i * n + p], miq = M[i * n + q];
                    M[i * n + p] = c * mip + s * miq;
                    M[i * n + q] = -s * mip + c * miq;
                }
                for (int j = 0; j < n; ++j) {
                    double mpj = M[p * n + j], mqj = M[q * n + j];
                    M[p * n + j] = c * mpj + s * mqj;
                    M[q * n + j] = -s * mpj + c * mqj;
                }

                // Accumulate eigenvector rotation
                if (V) {
                    for (int i = 0; i < n; ++i) {
                        double vip = V[i * n + p], viq = V[i * n + q];
                        V[i * n + p] = c * vip + s * viq;
                        V[i * n + q] = -s * vip + c * viq;
                    }
                }
            }
        }
    }

    // Did not converge — still return best estimate
    for (int i = 0; i < n; ++i)
        eigvals[i] = M[i * n + i];
    std::sort(eigvals, eigvals + n);
    return -1;
}

// Internal implementation shared by both API functions
static int bdg_solve_impl(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues, double* eigenvectors_out)
{
    if (!eigenvalues_out || !n_eigenvalues)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_sites <= 0 || n_sites > 2000)   // [2D-4] raised from 500
        return SUPERMAG_ERR_INVALID_DIM;

    int N = n_sites;
    int dim = 2 * N;
    *n_eigenvalues = dim;

    // Unify units: convert t_hop from eV to meV internally  [2D-2]
    double t_hop_meV = t_hop * 1e3;

    // Build BdG matrix (real symmetric for s-wave real Δ)  [EQ-10]
    // All quantities now in meV
    std::vector<double> H(dim * dim, 0.0);

    for (int i = 0; i < N; ++i) {
        // Electron block: on-site = -mu + E_ex
        H[i * dim + i] = -mu + E_ex;
        // Hole block: on-site = +mu + E_ex
        H[(N + i) * dim + (N + i)] = mu + E_ex;
        // Pairing
        H[i * dim + (N + i)] = Delta;
        H[(N + i) * dim + i] = Delta;
    }

    for (int i = 0; i < N - 1; ++i) {
        // Electron hopping
        H[i * dim + (i + 1)] = -t_hop_meV;
        H[(i + 1) * dim + i] = -t_hop_meV;
        // Hole hopping (sign flip for Nambu)
        H[(N + i) * dim + (N + i + 1)] = t_hop_meV;
        H[(N + i + 1) * dim + (N + i)] = t_hop_meV;
    }

    // Diagonalize with cyclic Jacobi
    int jrc = cyclic_jacobi(H.data(), dim, eigenvalues_out, eigenvectors_out);

    if (jrc != 0)
        return SUPERMAG_ERR_NO_CONVERGE;

    return SUPERMAG_OK;
}

extern "C" {

int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues,
    double* eigenvectors_out)
{
    return bdg_solve_impl(n_sites, t_hop, Delta, E_ex, mu,
                          eigenvalues_out, n_eigenvalues, eigenvectors_out);
}

}
