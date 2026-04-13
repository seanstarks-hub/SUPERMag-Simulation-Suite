// Unit tests for magnetic domain kernel.
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

void test_domains_single() {
    // Single domain (d_F == domain_width) should match bilayer kernel
    double d_F = 10.0, xi_F = 5.0, E_ex = 10.0;

    supermag_domain_params_t dom;
    dom.domain_width = d_F;  // N = d_F/domain_width = 1

    double Kr, Ki;
    int rc = supermag_proximity_kernel_domains(d_F, xi_F, E_ex, &dom,
                                               SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);

    auto K_ref = supermag::kernel_coth(d_F, xi_F);
    assert(std::abs(Kr - K_ref.real()) < 1e-10);
    assert(std::abs(Ki - K_ref.imag()) < 1e-10);

    // Pi phase
    rc = supermag_proximity_kernel_domains(d_F, xi_F, E_ex, &dom,
                                           SUPERMAG_PHASE_PI, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    auto K_pi = supermag::kernel_tanh(d_F, xi_F);
    assert(std::abs(Kr - K_pi.real()) < 1e-10);
    assert(std::abs(Ki - K_pi.imag()) < 1e-10);

    std::printf("  PASS: test_domains_single\n");
}

void test_domains_two_sharp() {
    // Two domains (d_F = 2*domain_width) — result should differ from bilayer
    double xi_F = 5.0, E_ex = 10.0;
    double d_dom = 5.0;
    double d_F = 2 * d_dom;

    supermag_domain_params_t dom;
    dom.domain_width = d_dom;  // N = 10/5 = 2

    double Kr, Ki;
    int rc = supermag_proximity_kernel_domains(d_F, xi_F, E_ex, &dom,
                                               SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    assert(std::isfinite(Kr) && std::isfinite(Ki));

    // Should differ from single-domain bilayer of same total thickness
    auto K_bilayer = supermag::kernel_coth(d_F, xi_F);
    double diff = std::abs(std::complex<double>(Kr, Ki) - K_bilayer);
    assert(diff > 1e-6);
    std::printf("  PASS: test_domains_two_sharp\n");
}

void test_domains_many() {
    // Many domains — should produce valid results
    double xi_F = 5.0, E_ex = 10.0;
    double d_dom = 10.0;
    double d_F = 3 * d_dom;

    supermag_domain_params_t dom;
    dom.domain_width = d_dom;  // N = 30/10 = 3

    double Kr, Ki;
    int rc = supermag_proximity_kernel_domains(d_F, xi_F, E_ex, &dom,
                                               SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_OK);
    assert(std::isfinite(Kr) && std::isfinite(Ki));

    // Compare with 2 domains — should differ
    supermag_domain_params_t dom2;
    dom2.domain_width = d_F / 2.0;  // N = 2
    double Kr2, Ki2;
    supermag_proximity_kernel_domains(d_F, xi_F, E_ex, &dom2,
                                      SUPERMAG_PHASE_ZERO, &Kr2, &Ki2);
    double diff = std::abs(std::complex<double>(Kr, Ki) - std::complex<double>(Kr2, Ki2));
    assert(diff > 1e-8);
    std::printf("  PASS: test_domains_many\n");
}

void test_domains_null() {
    double Kr, Ki;
    int rc = supermag_proximity_kernel_domains(10.0, 5.0, 10.0, nullptr,
                                               SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);

    supermag_domain_params_t dom = {10.0};
    rc = supermag_proximity_kernel_domains(10.0, 5.0, 10.0, &dom,
                                           SUPERMAG_PHASE_ZERO, nullptr, &Ki);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_domains_null\n");
}

void test_domains_invalid() {
    supermag_domain_params_t dom = {0.0};  // domain_width = 0 → invalid
    double Kr, Ki;
    int rc = supermag_proximity_kernel_domains(10.0, 5.0, 10.0, &dom,
                                               SUPERMAG_PHASE_ZERO, &Kr, &Ki);
    assert(rc == SUPERMAG_ERR_INVALID_DIM);
    std::printf("  PASS: test_domains_invalid\n");
}

int main() {
    std::printf("Running magnetic domain kernel tests...\n");
    test_domains_single();
    test_domains_two_sharp();
    test_domains_many();
    test_domains_null();
    test_domains_invalid();
    std::printf("All magnetic domain kernel tests passed!\n");
    return 0;
}
