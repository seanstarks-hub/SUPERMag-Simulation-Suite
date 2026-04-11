// BdG tight-binding Hamiltonian solver
// Constructs 2N×2N Nambu-doubled BdG matrix and diagonalizes.

#include "supermag/bdg.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

// Simple Jacobi eigenvalue algorithm for real symmetric matrices
// (production code would use LAPACK dsyev)
// Returns 0 if converged, -1 if max_iter reached without convergence.
static int jacobi_eigvals(double* A, int n, double* eigvals, int max_iter = 200) {
    // Work on a copy
    std::vector<double> M(n * n);
    std::memcpy(M.data(), A, n * n * sizeof(double));

    for (int iter = 0; iter < max_iter; ++iter) {
        // Find largest off-diagonal element
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (std::fabs(M[i * n + j]) > max_off) {
                    max_off = std::fabs(M[i * n + j]);
                    p = i; q = j;
                }
        if (max_off < 1e-12) {
            for (int i = 0; i < n; ++i)
                eigvals[i] = M[i * n + i];
            std::sort(eigvals, eigvals + n);
            return 0;  // converged
        }

        // Compute rotation
        double app = M[p * n + p], aqq = M[q * n + q], apq = M[p * n + q];
        double theta = 0.5 * std::atan2(2.0 * apq, app - aqq);
        double c = std::cos(theta), s = std::sin(theta);

        // Apply rotation
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
    }

    // Did not converge — still return best estimate
    for (int i = 0; i < n; ++i)
        eigvals[i] = M[i * n + i];
    std::sort(eigvals, eigvals + n);
    return -1;  // not converged
}

extern "C" {

int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues)
{
    if (!eigenvalues_out || !n_eigenvalues)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_sites <= 0 || n_sites > 500)
        return SUPERMAG_ERR_INVALID_DIM;

    int N = n_sites;
    int dim = 2 * N;
    *n_eigenvalues = dim;

    // Convert meV → eV
    double Delta_eV = Delta * 1e-3;
    double E_ex_eV = E_ex * 1e-3;
    double mu_eV = mu * 1e-3;

    // Build BdG matrix (real symmetric for s-wave real Δ)  [EQ-10]
    std::vector<double> H(dim * dim, 0.0);

    for (int i = 0; i < N; ++i) {
        // Electron block: on-site = -mu + E_ex  [EQ-10]
        H[i * dim + i] = -mu_eV + E_ex_eV;
        // Hole block: on-site = +mu + E_ex
        H[(N + i) * dim + (N + i)] = mu_eV + E_ex_eV;
        // Pairing
        H[i * dim + (N + i)] = Delta_eV;
        H[(N + i) * dim + i] = Delta_eV;
    }

    for (int i = 0; i < N - 1; ++i) {
        // Electron hopping
        H[i * dim + (i + 1)] = -t_hop;
        H[(i + 1) * dim + i] = -t_hop;
        // Hole hopping (sign flip)
        H[(N + i) * dim + (N + i + 1)] = t_hop;
        H[(N + i + 1) * dim + (N + i)] = t_hop;
    }

    // Diagonalize
    int jrc = jacobi_eigvals(H.data(), dim, eigenvalues_out);

    // Convert eV → meV
    for (int i = 0; i < dim; ++i)
        eigenvalues_out[i] *= 1e3;

    if (jrc != 0)
        return SUPERMAG_ERR_NO_CONVERGE;

    return SUPERMAG_OK;
}

}
