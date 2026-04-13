// Unit tests for S/N/F trilayer kernel.
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

void test_snf_reduces_to_sf() {
    // d_N → 0 should recover bilayer kernel
    double d_F = 10.0, xi_F = 5.0;

    supermag_trilayer_params_t tri;
    tri.d_N = 0.0;
    tri.xi_N = 10.0;
    tri.R_B = 0.0;

    double Kr, Ki;
    int rc = supermag_proximity_kernel_snf(d_F, xi_F, &tri,
                                           SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);

    auto K_bilayer = supermag::kernel_coth(d_F, xi_F);
    assert(std::abs(Kr - K_bilayer.real()) < 1e-12);
    assert(std::abs(Ki - K_bilayer.imag()) < 1e-12);

    // Pi phase too
    rc = supermag_proximity_kernel_snf(d_F, xi_F, &tri,
                                       SUPERMAG_PHASE_PI, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    auto K_pi = supermag::kernel_tanh(d_F, xi_F);
    assert(std::abs(Kr - K_pi.real()) < 1e-12);
    assert(std::abs(Ki - K_pi.imag()) < 1e-12);

    std::printf("  PASS: test_snf_reduces_to_sf\n");
}

void test_snf_thick_N() {
    // d_N >> xi_N: S layer sees only N metal, kernel → q_N = 1/xi_N
    double d_F = 10.0, xi_F = 5.0;

    supermag_trilayer_params_t tri;
    tri.d_N = 1000.0;   // very thick N
    tri.xi_N = 10.0;
    tri.R_B = 0.0;

    double Kr, Ki;
    int rc = supermag_proximity_kernel_snf(d_F, xi_F, &tri,
                                           SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);

    // For thick N: tanh(q_N*d_N) → 1, so
    //   K_SNF = q_N * (K_F + q_N) / (q_N + K_F) = q_N
    double q_N = 1.0 / tri.xi_N;
    assert(std::abs(Kr - q_N) < 1e-6);
    assert(std::abs(Ki) < 1e-6);
    std::printf("  PASS: test_snf_thick_N\n");
}

void test_snf_finite_N() {
    // Non-trivial case: finite N layer should modify kernel
    double d_F = 10.0, xi_F = 5.0;

    supermag_trilayer_params_t tri;
    tri.d_N = 5.0;
    tri.xi_N = 10.0;
    tri.R_B = 0.0;

    double Kr, Ki;
    int rc = supermag_proximity_kernel_snf(d_F, xi_F, &tri,
                                           SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    assert(std::isfinite(Kr));
    assert(std::isfinite(Ki));

    // Should differ from bilayer
    auto K_bilayer = supermag::kernel_coth(d_F, xi_F);
    double diff = std::abs(std::complex<double>(Kr, Ki) - K_bilayer);
    assert(diff > 1e-6);
    std::printf("  PASS: test_snf_finite_N\n");
}

void test_snf_null() {
    double Kr, Ki;
    int rc = supermag_proximity_kernel_snf(10.0, 5.0, nullptr,
                                           SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);

    supermag_trilayer_params_t tri = {5.0, 10.0, 0.0};
    rc = supermag_proximity_kernel_snf(10.0, 5.0, &tri,
                                       SUPERMAG_PHASE_ZERO, nullptr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_snf_null\n");
}

int main() {
    std::printf("Running S/N/F trilayer kernel tests...\n");
    test_snf_reduces_to_sf();
    test_snf_thick_N();
    test_snf_finite_N();
    test_snf_null();
    std::printf("All S/N/F trilayer kernel tests passed!\n");
    return 0;
}
