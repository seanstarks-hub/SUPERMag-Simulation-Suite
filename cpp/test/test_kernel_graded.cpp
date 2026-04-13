// Unit tests for graded ferromagnet kernel.
#include "supermag/proximity.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <complex>

// Forward-declare internal kernel for comparison
namespace supermag {
std::complex<double> kernel_coth(double d_F, double xi_F);
std::complex<double> kernel_tanh(double d_F, double xi_F);
}

void test_graded_uniform() {
    // Uniform exchange (E_left == E_right) should match bilayer kernel
    double d_F = 10.0, xi_F = 5.0, E_ex = 10.0;

    supermag_graded_params_t grade;
    grade.E_ex_surface = E_ex;
    grade.E_ex_bulk = E_ex;
    grade.profile = SUPERMAG_GRADE_LINEAR;
    grade.n_slabs = 50;

    double Kr, Ki;
    int rc = supermag_proximity_kernel_graded(d_F, xi_F, &grade,
                                              SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);

    auto K_ref = supermag::kernel_coth(d_F, xi_F);
    assert(std::abs(Kr - K_ref.real()) < 1e-8);
    assert(std::abs(Ki - K_ref.imag()) < 1e-8);

    // Pi phase too
    rc = supermag_proximity_kernel_graded(d_F, xi_F, &grade,
                                          SUPERMAG_PHASE_PI, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    auto K_pi = supermag::kernel_tanh(d_F, xi_F);
    assert(std::abs(Kr - K_pi.real()) < 1e-8);
    assert(std::abs(Ki - K_pi.imag()) < 1e-8);

    std::printf("  PASS: test_graded_uniform\n");
}

void test_graded_convergence() {
    // Increasing n_slices should converge
    double d_F = 10.0, xi_F = 5.0;

    supermag_graded_params_t grade;
    grade.E_ex_surface = 10.0;
    grade.E_ex_bulk = 50.0;
    grade.profile = SUPERMAG_GRADE_LINEAR;

    double prev_Kr = 0.0, prev_Ki = 0.0;
    double Kr, Ki;

    for (int n : {10, 50, 100, 500}) {
        grade.n_slabs = n;
        int rc = supermag_proximity_kernel_graded(d_F, xi_F, &grade,
                                                  SUPERMAG_PHASE_ZERO, &Kr, &Ki);
        assert(rc == SUPERMAG_OK);
        assert(std::isfinite(Kr) && std::isfinite(Ki));

        if (n > 10) {
            // Change should decrease with more slices
            double delta_r = std::abs(Kr - prev_Kr);
            double delta_i = std::abs(Ki - prev_Ki);
            assert(delta_r < 0.1 && delta_i < 0.1);  // reasonable convergence
        }
        prev_Kr = Kr;
        prev_Ki = Ki;
    }
    std::printf("  PASS: test_graded_convergence\n");
}

void test_graded_profiles() {
    // All three profiles should produce finite results
    double d_F = 10.0, xi_F = 5.0;

    supermag_graded_params_t grade;
    grade.E_ex_surface = 10.0;
    grade.E_ex_bulk = 30.0;
    grade.n_slabs = 50;

    for (auto prof : {SUPERMAG_GRADE_LINEAR, SUPERMAG_GRADE_EXPONENTIAL, SUPERMAG_GRADE_STEP}) {
        grade.profile = prof;
        double Kr, Ki;
        int rc = supermag_proximity_kernel_graded(d_F, xi_F, &grade,
                                                  SUPERMAG_PHASE_ZERO, &Kr, &Ki);
        assert(rc == SUPERMAG_OK);
        assert(std::isfinite(Kr) && std::isfinite(Ki));
    }
    std::printf("  PASS: test_graded_profiles\n");
}

void test_graded_null() {
    double Kr, Ki;
    int rc = supermag_proximity_kernel_graded(10.0, 5.0, nullptr,
                                              SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);

    supermag_graded_params_t grade = {10.0, 30.0, SUPERMAG_GRADE_LINEAR, 50};
    rc = supermag_proximity_kernel_graded(10.0, 5.0, &grade,
                                          SUPERMAG_PHASE_ZERO, nullptr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_graded_null\n");
}

void test_graded_invalid() {
    supermag_graded_params_t grade = {10.0, 30.0, SUPERMAG_GRADE_LINEAR, 0};
    double Kr, Ki;
    // n_slabs = 0 → INVALID_DIM
    int rc = supermag_proximity_kernel_graded(10.0, 5.0, &grade,
                                              SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_ERR_INVALID_DIM);
    std::printf("  PASS: test_graded_invalid\n");
}

int main() {
    std::printf("Running graded ferromagnet kernel tests...\n");
    test_graded_uniform();
    test_graded_convergence();
    test_graded_profiles();
    test_graded_null();
    test_graded_invalid();
    std::printf("All graded ferromagnet kernel tests passed!\n");
    return 0;
}
