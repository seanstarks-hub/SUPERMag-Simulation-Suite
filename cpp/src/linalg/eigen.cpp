/**
 * Eigenvalue solver utilities.
 *
 * Provides a Jacobi eigenvalue algorithm for real symmetric matrices.
 * Used by the BdG solver for small-to-moderate Hamiltonians.
 */

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace supermag {

/**
 * Jacobi eigenvalue algorithm for a real symmetric matrix.
 *
 * Overwrites the lower triangle of A in place. Eigenvalues are written
 * to eigenvalues[] in ascending order. Matrix is stored in row-major.
 *
 * @param A          Row-major n×n symmetric matrix (overwritten).
 * @param n          Matrix dimension.
 * @param eigenvalues Output array of length n.
 * @param max_iter   Maximum Jacobi sweeps (default 200).
 * @return           Number of sweeps performed, or -1 on failure.
 */
int jacobi_eigenvalues(double* A, int n, double* eigenvalues, int max_iter) {
    if (!A || !eigenvalues || n <= 0) return -1;
    if (n == 1) {
        eigenvalues[0] = A[0];
        return 0;
    }

    const double tol = 1e-12;

    for (int sweep = 0; sweep < max_iter; ++sweep) {
        // Compute off-diagonal norm
        double off_norm = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                off_norm += A[i * n + j] * A[i * n + j];

        if (off_norm < tol) {
            // Extract diagonal eigenvalues and sort
            for (int i = 0; i < n; ++i)
                eigenvalues[i] = A[i * n + i];
            std::sort(eigenvalues, eigenvalues + n);
            return sweep;
        }

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double Apq = A[p * n + q];
                if (std::abs(Apq) < tol * 0.01) continue;

                double tau = (A[q * n + q] - A[p * n + p]) / (2.0 * Apq);
                double t = (tau >= 0 ? 1.0 : -1.0) /
                           (std::abs(tau) + std::sqrt(1.0 + tau * tau));
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s = t * c;

                // Apply rotation
                A[p * n + p] -= t * Apq;
                A[q * n + q] += t * Apq;
                A[p * n + q] = 0.0;
                A[q * n + p] = 0.0;

                for (int r = 0; r < n; ++r) {
                    if (r == p || r == q) continue;
                    double Arp = A[r * n + p];
                    double Arq = A[r * n + q];
                    A[r * n + p] = c * Arp - s * Arq;
                    A[p * n + r] = A[r * n + p];
                    A[r * n + q] = s * Arp + c * Arq;
                    A[q * n + r] = A[r * n + q];
                }
            }
        }
    }

    // Did not converge — still extract diagonal
    for (int i = 0; i < n; ++i)
        eigenvalues[i] = A[i * n + i];
    std::sort(eigenvalues, eigenvalues + n);
    return -1;
}

} // namespace supermag
