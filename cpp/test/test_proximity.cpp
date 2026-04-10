// Unit tests for proximity solver (pair_amplitude and critical_temp)
#include "supermag/proximity.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_pair_amplitude_basic() {
    const int n = 100;
    std::vector<double> x(n), F(n);
    int rc = supermag_proximity_pair_amplitude(1.0, 5.0, 20.0, n, x.data(), F.data());
    assert(rc == SUPERMAG_OK);
    // F(0) should equal F0
    assert(std::abs(F[0] - 1.0) < 1e-12);
    // F should decay
    assert(std::abs(F[n-1]) < std::abs(F[0]));
    // x should span [0, d_F]
    assert(std::abs(x[0]) < 1e-12);
    assert(std::abs(x[n-1] - 20.0) < 1e-12);
    std::printf("  PASS: test_pair_amplitude_basic\n");
}

void test_pair_amplitude_null() {
    int rc = supermag_proximity_pair_amplitude(1.0, 5.0, 20.0, 10, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_pair_amplitude_null\n");
}

void test_pair_amplitude_invalid() {
    double x[2], F[2];
    int rc = supermag_proximity_pair_amplitude(1.0, 5.0, -1.0, 2, x, F);
    assert(rc == SUPERMAG_ERR_INVALID_DIM);
    std::printf("  PASS: test_pair_amplitude_invalid\n");
}

void test_critical_temp_basic() {
    const int n = 50;
    std::vector<double> d_F(n), Tc(n);
    for (int i = 0; i < n; ++i) d_F[i] = 0.5 + i * 0.5;

    int rc = supermag_proximity_critical_temp(9.2, 50.0, 38.0, 0.7, 256.0,
                                              d_F.data(), n, Tc.data());
    assert(rc == SUPERMAG_OK);
    // All Tc values should be non-negative and <= Tc0
    for (int i = 0; i < n; ++i) {
        assert(Tc[i] >= 0.0);
        assert(Tc[i] <= 9.2 + 1e-10);
    }
    std::printf("  PASS: test_critical_temp_basic\n");
}

void test_critical_temp_null() {
    double Tc;
    int rc = supermag_proximity_critical_temp(9.2, 50.0, 38.0, 0.7, 256.0, nullptr, 1, &Tc);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_critical_temp_null\n");
}

int main() {
    std::printf("Running proximity tests...\n");
    test_pair_amplitude_basic();
    test_pair_amplitude_null();
    test_pair_amplitude_invalid();
    test_critical_temp_basic();
    test_critical_temp_null();
    std::printf("All proximity tests passed!\n");
    return 0;
}
