#include "supermag/error.h"
#include <cmath>
#include <vector>

extern "C" {

/**
 * Solve a tridiagonal system Ax = d using the Thomas algorithm.
 * a[0..n-1]: sub-diagonal (a[0] unused)
 * b[0..n-1]: main diagonal
 * c[0..n-1]: super-diagonal (c[n-1] unused)
 * d[0..n-1]: right-hand side (overwritten with solution)
 * n: system size
 */
int supermag_tridiag_solve(
    const double* a, const double* b, const double* c,
    double* d, int n)
{
    if (!a || !b || !c || !d)
        return SUPERMAG_ERR_NULL_PTR;
    if (n < 1)
        return SUPERMAG_ERR_INVALID_DIM;
    if (n == 1) {
        if (std::abs(b[0]) < 1e-300) return SUPERMAG_ERR_NO_CONVERGE;
        d[0] /= b[0];
        return SUPERMAG_OK;
    }

    // Working copies of b (modified during forward sweep)
    std::vector<double> b_prime(n);
    b_prime[0] = b[0];
    if (std::abs(b_prime[0]) < 1e-300) return SUPERMAG_ERR_NO_CONVERGE;

    // Forward elimination
    for (int i = 1; i < n; ++i) {
        double m = a[i] / b_prime[i - 1];
        b_prime[i] = b[i] - m * c[i - 1];
        d[i] = d[i] - m * d[i - 1];
        if (std::abs(b_prime[i]) < 1e-300) return SUPERMAG_ERR_NO_CONVERGE;
    }

    // Back substitution
    d[n - 1] /= b_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        d[i] = (d[i] - c[i] * d[i + 1]) / b_prime[i];
    }

    return SUPERMAG_OK;
}

}
