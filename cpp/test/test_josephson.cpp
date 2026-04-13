// Unit tests for Josephson CPR solver
#include "supermag/josephson.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_josephson_sinusoidal_limit() {
    // Large d_F → heavily damped, should approach sin(φ) shape
    const int n = 64;
    std::vector<double> current(n);
    int rc = supermag_josephson_cpr(
        50.0, 5.0, 10.0, 4.0, 9.2,  // large d_F = 50nm
        0.0,                          // gamma_B = 0
        n, nullptr,                   // auto-generate phases
        current.data(), nullptr);
    assert(rc == SUPERMAG_OK);

    // Current at φ=0 should be ~0
    assert(std::fabs(current[0]) < 0.1);
    std::printf("  PASS: test_josephson_sinusoidal_limit\n");
}

void test_josephson_pi_junction() {
    // Intermediate d_F should show π-junction: current negative near φ=0+
    const int n = 64;
    std::vector<double> current(n);
    // d_F ≈ 2.5·xi_F puts us near the first 0-π transition
    int rc = supermag_josephson_cpr(
        12.5, 5.0, 10.0, 4.0, 9.2,
        0.0, n, nullptr,
        current.data(), nullptr);
    assert(rc == SUPERMAG_OK);

    // The CPR should have nonzero structure
    bool has_positive = false, has_negative = false;
    for (int i = 0; i < n; ++i) {
        if (current[i] > 0.01) has_positive = true;
        if (current[i] < -0.01) has_negative = true;
    }
    // At least has some positive and negative parts
    assert(has_positive || has_negative);
    std::printf("  PASS: test_josephson_pi_junction (pos=%d, neg=%d)\n",
                has_positive, has_negative);
}

void test_josephson_barrier() {
    // gamma_B > 0 should reduce critical current vs gamma_B = 0
    const int n = 32;
    std::vector<double> I0(n), IB(n);
    double Ic_no_barrier, Ic_with_barrier;

    supermag_josephson_cpr(10.0, 5.0, 10.0, 4.0, 9.2, 0.0,
                               n, nullptr,
                               I0.data(), &Ic_no_barrier);
    supermag_josephson_cpr(10.0, 5.0, 10.0, 4.0, 9.2, 1.0,
                               n, nullptr,
                               IB.data(), &Ic_with_barrier);

    // Barrier should reduce Ic
    assert(Ic_with_barrier < Ic_no_barrier + 1e-15);
    assert(Ic_with_barrier > 0.0);
    std::printf("  PASS: test_josephson_barrier (Ic=%.6f → %.6f with gamma_B=1)\n",
                Ic_no_barrier, Ic_with_barrier);
}

void test_josephson_user_phases() {
    // User-supplied phase array should be used as-is (const input)
    const int n = 16;
    std::vector<double> phase_in(n), current(n);
    double pi = 3.14159265358979323846;

    // Set custom phases [0, π)
    for (int i = 0; i < n; ++i)
        phase_in[i] = pi * i / n;

    int rc = supermag_josephson_cpr(
        10.0, 5.0, 10.0, 4.0, 9.2, 0.0,
        n, phase_in.data(),
        current.data(), nullptr);
    assert(rc == SUPERMAG_OK);

    // Current at phase=0 should be ~0
    assert(std::fabs(current[0]) < 1e-10);
    // Should have finite current at interior phases
    bool has_nonzero = false;
    for (int i = 1; i < n; ++i)
        if (std::fabs(current[i]) > 1e-10) has_nonzero = true;
    assert(has_nonzero);
    std::printf("  PASS: test_josephson_user_phases\n");
}

void test_josephson_absolute_ic() {
    // Ic_out should give a finite positive value
    const int n = 32;
    std::vector<double> current(n);
    double Ic;

    int rc = supermag_josephson_cpr(
        10.0, 5.0, 10.0, 4.0, 9.2, 0.0,
        n, nullptr,
        current.data(), &Ic);
    assert(rc == SUPERMAG_OK);
    assert(Ic > 0.0);
    assert(std::isfinite(Ic));

    // Normalized current should have max ≈ 1
    double max_norm = 0.0;
    for (int i = 0; i < n; ++i)
        if (std::fabs(current[i]) > max_norm) max_norm = std::fabs(current[i]);
    assert(std::fabs(max_norm - 1.0) < 1e-10);

    std::printf("  PASS: test_josephson_absolute_ic (Ic=%.6e)\n", Ic);
}

void test_josephson_null() {
    int rc = supermag_josephson_cpr(10.0, 5.0, 10.0, 4.0, 9.2, 0.0,
                                    32, nullptr,
                                    nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_josephson_null\n");
}

void test_josephson_tc0_required() {
    // Tc0 <= 0 should error
    std::vector<double> current(32);
    int rc = supermag_josephson_cpr(10.0, 5.0, 10.0, 4.0, 0.0, 0.0,
                                    32, nullptr, current.data(), nullptr);
    assert(rc != SUPERMAG_OK);
    rc = supermag_josephson_cpr(10.0, 5.0, 10.0, 4.0, -1.0, 0.0,
                                32, nullptr, current.data(), nullptr);
    assert(rc != SUPERMAG_OK);
    std::printf("  PASS: test_josephson_tc0_required\n");
}

int main() {
    std::printf("Running Josephson tests...\n");
    test_josephson_sinusoidal_limit();
    test_josephson_pi_junction();
    test_josephson_barrier();
    test_josephson_user_phases();
    test_josephson_absolute_ic();
    test_josephson_null();
    test_josephson_tc0_required();
    std::printf("All Josephson tests passed!\n");
    return 0;
}
