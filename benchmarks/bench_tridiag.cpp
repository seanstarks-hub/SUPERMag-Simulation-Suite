// Benchmark for tridiagonal solver
#include "supermag/error.h"
#include <chrono>
#include <cstdio>
#include <vector>

extern "C" int supermag_tridiag_solve(
    const double* a, const double* b, const double* c, double* d, int n);

int main() {
    const int n = 10000;
    const int n_iterations = 1000;

    std::vector<double> a(n, -1.0), b(n, 2.0), c(n, -1.0), d(n, 1.0);
    a[0] = 0; c[n-1] = 0;

    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < n_iterations; ++iter) {
        std::vector<double> d_copy = d;
        supermag_tridiag_solve(a.data(), b.data(), c.data(), d_copy.data(), n);
    }
    auto end = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::printf("tridiag_solve: %d iterations x %d unknowns\n", n_iterations, n);
    std::printf("  Total: %.2f ms\n", ms);
    std::printf("  Per call: %.3f us\n", ms * 1000.0 / n_iterations);
    return 0;
}
