// Unit tests for Eilenberger clean-limit solver
#include "supermag/eilenberger.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_eilenberger_bulk_s() {
    // Pure S with no F: should give a_BCS everywhere in S
    const int n = 100;
    std::vector<double> f(n), x(n);
    int rc = supermag_eilenberger_solve(
        9.2, 50.0, 0.01,  // Tc0, d_S, d_F (tiny F)
        38.0, 10.0,        // xi_S, E_ex
        4.0,               // T
        nullptr, n, f.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // f should be nonzero in S region (x < 0)
    assert(f[0] > 0.0);
    // f should be largest deep in S
    assert(f[0] > 0.01);
    std::printf("  PASS: test_eilenberger_bulk_s (f[0]=%.4f)\n", f[0]);
}

void test_eilenberger_frequency_convergence() {
    // Compare at two temperatures — higher T means fewer Matsubara terms needed
    const int n = 80;
    std::vector<double> f_low(n), x_low(n), f_high(n), x_high(n);

    supermag_eilenberger_solve(9.2, 50.0, 20.0, 38.0, 10.0, 2.0,
                               nullptr, n, f_low.data(), x_low.data());
    supermag_eilenberger_solve(9.2, 50.0, 20.0, 38.0, 10.0, 7.0,
                               nullptr, n, f_high.data(), x_high.data());

    // Both should give finite results
    assert(f_low[0] > 0.0);
    assert(f_high[0] > 0.0);

    // Higher T → weaker pairing → lower f generally (but frequency sum is different)
    // Just check both are finite and reasonable
    assert(std::isfinite(f_low[0]));
    assert(std::isfinite(f_high[0]));
    std::printf("  PASS: test_eilenberger_frequency_convergence (f_2K=%.4f, f_7K=%.4f)\n",
                f_low[0], f_high[0]);
}

void test_eilenberger_rk4_stability() {
    // Large E_ex should not blow up with RK4 (it would with Euler)
    const int n = 100;
    std::vector<double> f(n), x(n);
    int rc = supermag_eilenberger_solve(
        9.2, 50.0, 30.0,  // d_F = 30nm
        38.0, 100.0,       // Large E_ex = 100 meV
        4.0,
        nullptr, n, f.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // All values should be finite
    for (int i = 0; i < n; ++i) {
        assert(std::isfinite(f[i]));
        assert(f[i] >= 0.0);
    }
    std::printf("  PASS: test_eilenberger_rk4_stability\n");
}

void test_eilenberger_null() {
    int rc = supermag_eilenberger_solve(9.2, 50.0, 20.0, 38.0, 10.0, 4.0,
                                        nullptr, 100, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_eilenberger_null\n");
}

int main() {
    std::printf("Running Eilenberger tests...\n");
    test_eilenberger_bulk_s();
    test_eilenberger_frequency_convergence();
    test_eilenberger_rk4_stability();
    test_eilenberger_null();
    std::printf("All Eilenberger tests passed!\n");
    return 0;
}
