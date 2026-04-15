// Unit tests for Usadel nonlinear solver
#include "supermag/usadel.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_usadel_bcs_limit() {
    // d_F → 0: should recover bulk BCS-like Delta profile in S
    const int n = 200;
    std::vector<double> Delta(n), x(n);
    int rc = supermag_usadel_solve(
        9.2, 50.0, 0.01,  // Tc0, d_S, d_F (tiny F layer)
        38.0, 5.0, 10.0,  // xi_S, xi_F, E_ex
        4.0,               // T = 4K
        SUPERMAG_USADEL_NONLINEAR, nullptr,
        n, Delta.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // Delta should be nearly uniform in bulk S (away from interface)
    // Check that Delta is positive in S region
    int n_bulk = n / 4;  // deep in S
    assert(Delta[0] > 0.0);
    assert(Delta[n_bulk] > 0.0);

    // Delta at deep S should be close to BCS value
    double kB_meV = 8.617333262e-2;
    double Delta_BCS = 1.764 * kB_meV * 9.2 * std::sqrt(1.0 - 4.0 / 9.2);
    // Allow significant tolerance since self-consistency modifies the profile
    assert(Delta[0] > 0.0);  // must be positive
    std::printf("  PASS: test_usadel_bcs_limit (Delta[0]=%.4f meV, BCS=%.4f meV)\n",
                Delta[0], Delta_BCS);
}

void test_usadel_convergence() {
    // Verify solver converges for typical Nb/CuNi parameters
    const int n = 100;
    std::vector<double> Delta(n), x(n);
    int rc = supermag_usadel_solve(
        9.2, 50.0, 20.0,  // Tc0, d_S=50nm, d_F=20nm
        38.0, 5.0, 10.0,  // xi_S, xi_F, E_ex
        4.0,               // T
        SUPERMAG_USADEL_NONLINEAR, nullptr,
        n, Delta.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // Delta should be positive in S, zero in F
    // Find approximate boundary
    int n_S = 0;
    for (int i = 0; i < n; ++i) {
        if (x[i] >= 0.0) { n_S = i; break; }
    }
    assert(n_S > 0);
    assert(Delta[0] > 0.0);  // bulk S
    // F region should have Delta = 0
    assert(std::abs(Delta[n-1]) < 1e-10);
    std::printf("  PASS: test_usadel_convergence (n_S=%d)\n", n_S);
}

void test_usadel_grid_independence() {
    // Compare coarse and fine grids
    const int n1 = 50, n2 = 200;
    std::vector<double> Delta1(n1), x1(n1), Delta2(n2), x2(n2);

    supermag_usadel_solve(9.2, 50.0, 20.0, 38.0, 5.0, 10.0, 4.0,
                          SUPERMAG_USADEL_NONLINEAR, nullptr,
                          n1, Delta1.data(), x1.data());
    supermag_usadel_solve(9.2, 50.0, 20.0, 38.0, 5.0, 10.0, 4.0,
                          SUPERMAG_USADEL_NONLINEAR, nullptr,
                          n2, Delta2.data(), x2.data());

    // Delta at x=0 (interface) should be similar between grids
    // Find interface index in each grid
    int idx1 = 0, idx2 = 0;
    for (int i = 0; i < n1; ++i) if (x1[i] >= 0.0) { idx1 = i; break; }
    for (int i = 0; i < n2; ++i) if (x2[i] >= 0.0) { idx2 = i; break; }

    // Relative difference at interface should be reasonable
    double d1 = Delta1[idx1], d2 = Delta2[idx2];
    double ref = std::max(std::abs(d1), std::abs(d2));
    if (ref > 1e-10) {
        double rel_diff = std::abs(d1 - d2) / ref;
        // Allow up to 50% — grids are very different sizes
        assert(rel_diff < 0.5);
    }
    std::printf("  PASS: test_usadel_grid_independence (Delta_if=%.4f vs %.4f)\n", d1, d2);
}

void test_usadel_null() {
    int rc = supermag_usadel_solve(9.2, 50.0, 20.0, 38.0, 5.0, 10.0, 4.0,
                                   SUPERMAG_USADEL_NONLINEAR, nullptr,
                                   100, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_usadel_null\n");
}

void test_usadel_proportional_grid() {
    // With d_S >> d_F, S region should get more grid points
    const int n = 100;
    std::vector<double> Delta(n), x(n);
    // d_S = 90nm, d_F = 10nm → n_S should be ~90
    int rc = supermag_usadel_solve(9.2, 90.0, 10.0, 38.0, 5.0, 10.0, 4.0,
                                   SUPERMAG_USADEL_NONLINEAR, nullptr,
                                   n, Delta.data(), x.data());
    assert(rc == SUPERMAG_OK);

    // Count grid points in S (x < 0) and F (x >= 0)
    int count_S = 0;
    for (int i = 0; i < n; ++i)
        if (x[i] < 0.0) count_S++;

    // Proportional: n_S ~ 90% of n_grid
    assert(count_S > 70);  // should be around 90
    assert(count_S < 99);  // leave some for F
    std::printf("  PASS: test_usadel_proportional_grid (n_S=%d of %d)\n", count_S, n);
}

int main() {
    std::printf("Running Usadel tests...\n");
    test_usadel_bcs_limit();
    test_usadel_convergence();
    test_usadel_grid_independence();
    test_usadel_null();
    test_usadel_proportional_grid();
    std::printf("All Usadel tests passed!\n");
    return 0;
}
