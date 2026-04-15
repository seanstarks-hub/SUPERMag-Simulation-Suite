// Unit tests for triplet superconductivity solver
#include "supermag/triplet.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_triplet_collinear() {
    // Parallel magnetizations → f_triplet = 0
    const int n_layers = 3;
    double thick[] = {5.0, 10.0, 5.0};
    double angles[] = {0.0, 0.0, 0.0};  // all collinear

    const int n_grid = 100;
    std::vector<double> f(n_grid), x(n_grid);
    int rc = supermag_triplet_solve(n_layers, thick, angles,
                                    nullptr, nullptr,
                                    1.0, 10.0, 4.0, 9.2,
                                    SUPERMAG_TRIPLET_DIFFUSIVE,
                                    n_grid, f.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // All triplet amplitudes should be zero (sin(0) = 0)
    for (int i = 0; i < n_grid; ++i)
        assert(std::fabs(f[i]) < 1e-12);

    std::printf("  PASS: test_triplet_collinear\n");
}

void test_triplet_perpendicular() {
    // α = π/2 → maximum triplet conversion
    const double pi = 3.14159265358979323846;
    const int n_layers = 3;
    double thick[] = {5.0, 10.0, 5.0};
    double angles[] = {0.0, pi / 2.0, 0.0};  // 90° misalignment

    const int n_grid = 100;
    std::vector<double> f_perp(n_grid), x_perp(n_grid);
    int rc = supermag_triplet_solve(n_layers, thick, angles,
                                    nullptr, nullptr,
                                    1.0, 10.0, 4.0, 9.2,
                                    SUPERMAG_TRIPLET_DIFFUSIVE,
                                    n_grid, f_perp.data(), x_perp.data());
    assert(rc == SUPERMAG_OK);

    // Should have nonzero triplet at interfaces
    double max_f = 0.0;
    for (int i = 0; i < n_grid; ++i)
        if (f_perp[i] > max_f) max_f = f_perp[i];
    assert(max_f > 0.0);

    // Compare with smaller angle — 90° should give larger triplet
    double angles_small[] = {0.0, pi / 6.0, 0.0};  // 30°
    std::vector<double> f_small(n_grid), x_small(n_grid);
    supermag_triplet_solve(n_layers, thick, angles_small,
                           nullptr, nullptr,
                           1.0, 10.0, 4.0, 9.2,
                           SUPERMAG_TRIPLET_DIFFUSIVE,
                           n_grid, f_small.data(), x_small.data());
    double max_small = 0.0;
    for (int i = 0; i < n_grid; ++i)
        if (f_small[i] > max_small) max_small = f_small[i];

    assert(max_f > max_small);
    std::printf("  PASS: test_triplet_perpendicular (max_90=%.4f > max_30=%.4f)\n",
                max_f, max_small);
}

void test_triplet_decay_length() {
    // f_triplet should decay with xi_N away from interfaces
    const double pi = 3.14159265358979323846;
    const int n_layers = 2;
    double thick[] = {10.0, 40.0};  // interface at x=10
    double angles[] = {0.0, pi / 2.0};

    const int n_grid = 200;
    std::vector<double> f(n_grid), x(n_grid);
    double xi_N = 10.0;
    int rc = supermag_triplet_solve(n_layers, thick, angles,
                                    nullptr, nullptr,
                                    1.0, xi_N, 4.0, 9.2,
                                    SUPERMAG_TRIPLET_DIFFUSIVE,
                                    n_grid, f.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // Find amplitude at interface and at distance 2*xi_N
    int idx_int = 0, idx_far = 0;
    for (int i = 0; i < n_grid; ++i) {
        if (std::fabs(x[i] - 10.0) < std::fabs(x[idx_int] - 10.0)) idx_int = i;
        if (std::fabs(x[i] - 30.0) < std::fabs(x[idx_far] - 30.0)) idx_far = i;
    }

    // At 2·xi_N from interface, amplitude should be ~e^{-2} ≈ 0.14 of interface value
    if (f[idx_int] > 1e-10) {
        double ratio = f[idx_far] / f[idx_int];
        assert(ratio < 0.5);  // significant decay
        assert(ratio > 0.001);  // but not zero (long-range)
    }
    std::printf("  PASS: test_triplet_decay_length\n");
}

void test_triplet_xi_F_effect() {
    // xi_F should control short-range decay (now actually used)
    const double pi = 3.14159265358979323846;
    const int n_layers = 2;
    double thick[] = {10.0, 30.0};
    double angles[] = {0.0, pi / 2.0};

    const int n_grid = 100;
    std::vector<double> f1(n_grid), x1(n_grid), f2(n_grid), x2(n_grid);

    // Small xi_F → faster short-range decay
    supermag_triplet_solve(n_layers, thick, angles,
                           nullptr, nullptr,
                           0.5, 10.0, 4.0, 9.2,
                           SUPERMAG_TRIPLET_DIFFUSIVE,
                           n_grid, f1.data(), x1.data());
    // Large xi_F → slower short-range decay
    supermag_triplet_solve(n_layers, thick, angles,
                           nullptr, nullptr,
                           5.0, 10.0, 4.0, 9.2,
                           SUPERMAG_TRIPLET_DIFFUSIVE,
                           n_grid, f2.data(), x2.data());

    // With larger xi_F, the total triplet near the interface should differ
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < n_grid; ++i) { sum1 += f1[i]; sum2 += f2[i]; }
    assert(sum1 > 0.0);
    assert(sum2 > 0.0);
    assert(std::fabs(sum1 - sum2) > 1e-6);
    std::printf("  PASS: test_triplet_xi_F_effect (sum_0.5=%.4f, sum_5=%.4f)\n", sum1, sum2);
}

void test_triplet_temperature() {
    // Higher T → weaker triplet (through Delta(T) scaling)
    const double pi = 3.14159265358979323846;
    const int n_layers = 2;
    double thick[] = {10.0, 20.0};
    double angles[] = {0.0, pi / 2.0};

    const int n_grid = 100;
    std::vector<double> f_low(n_grid), x_low(n_grid), f_high(n_grid), x_high(n_grid);

    supermag_triplet_solve(n_layers, thick, angles,
                               nullptr, nullptr,
                               1.0, 10.0, 2.0, 9.2,
                               SUPERMAG_TRIPLET_DIFFUSIVE,
                               n_grid, f_low.data(), x_low.data());
    supermag_triplet_solve(n_layers, thick, angles,
                               nullptr, nullptr,
                               1.0, 10.0, 8.0, 9.2,
                               SUPERMAG_TRIPLET_DIFFUSIVE,
                               n_grid, f_high.data(), x_high.data());

    double max_low = 0, max_high = 0;
    for (int i = 0; i < n_grid; ++i) {
        if (f_low[i] > max_low) max_low = f_low[i];
        if (f_high[i] > max_high) max_high = f_high[i];
    }

    // Lower T → stronger triplet
    assert(max_low > max_high);
    std::printf("  PASS: test_triplet_temperature (max_2K=%.4f > max_8K=%.4f)\n",
                max_low, max_high);
}

void test_triplet_null() {
    int rc = supermag_triplet_solve(2, nullptr, nullptr,
                                    nullptr, nullptr,
                                    1.0, 10.0, 4.0, 9.2,
                                    SUPERMAG_TRIPLET_DIFFUSIVE,
                                    100, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_triplet_null\n");
}

int main() {
    std::printf("Running Triplet tests...\n");
    test_triplet_collinear();
    test_triplet_perpendicular();
    test_triplet_decay_length();
    test_triplet_xi_F_effect();
    test_triplet_temperature();
    test_triplet_null();
    std::printf("All Triplet tests passed!\n");
    return 0;
}
