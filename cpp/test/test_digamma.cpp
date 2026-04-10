// Unit tests for the complex digamma function
#include "../src/common/digamma.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <complex>

static const double EULER_GAMMA = 0.5772156649015329;

void test_digamma_at_one() {
    // psi(1) = -gamma_Euler
    auto result = supermag::digamma({1.0, 0.0});
    assert(std::abs(result.real() - (-EULER_GAMMA)) < 1e-10);
    assert(std::abs(result.imag()) < 1e-10);
    std::printf("  PASS: test_digamma_at_one\n");
}

void test_digamma_at_half() {
    // psi(1/2) = -gamma_Euler - 2*ln(2)
    double expected = -EULER_GAMMA - 2.0 * std::log(2.0);
    auto result = supermag::digamma({0.5, 0.0});
    assert(std::abs(result.real() - expected) < 1e-10);
    assert(std::abs(result.imag()) < 1e-10);
    std::printf("  PASS: test_digamma_at_half\n");
}

void test_digamma_recurrence() {
    // psi(z+1) = psi(z) + 1/z for arbitrary complex z
    std::complex<double> z(3.7, 1.2);
    auto psi_z = supermag::digamma(z);
    auto psi_z1 = supermag::digamma(z + 1.0);
    auto expected = psi_z + 1.0 / z;
    assert(std::abs(psi_z1 - expected) < 1e-10);
    std::printf("  PASS: test_digamma_recurrence\n");
}

void test_digamma_integer_values() {
    // psi(n) = -gamma + sum_{k=1}^{n-1} 1/k  (harmonic numbers)
    // psi(2) = -gamma + 1 = 1 - gamma
    auto result = supermag::digamma({2.0, 0.0});
    double expected = 1.0 - EULER_GAMMA;
    assert(std::abs(result.real() - expected) < 1e-10);
    // psi(3) = -gamma + 1 + 1/2
    result = supermag::digamma({3.0, 0.0});
    expected = 1.5 - EULER_GAMMA;
    assert(std::abs(result.real() - expected) < 1e-10);
    std::printf("  PASS: test_digamma_integer_values\n");
}

void test_digamma_complex() {
    // Test complex argument: verify recurrence holds for purely imaginary shift
    std::complex<double> z(0.5, 2.0);
    auto psi_z = supermag::digamma(z);
    auto psi_z1 = supermag::digamma(z + 1.0);
    auto diff = psi_z1 - psi_z - 1.0 / z;
    assert(std::abs(diff) < 1e-10);
    std::printf("  PASS: test_digamma_complex\n");
}

void test_digamma_large_argument() {
    // For large |z|, psi(z) ~ ln(z) - 1/(2z)
    std::complex<double> z(100.0, 0.0);
    auto result = supermag::digamma(z);
    double approx = std::log(100.0) - 0.005;
    assert(std::abs(result.real() - approx) < 0.001);
    std::printf("  PASS: test_digamma_large_argument\n");
}

int main() {
    std::printf("Running digamma tests...\n");
    test_digamma_at_one();
    test_digamma_at_half();
    test_digamma_recurrence();
    test_digamma_integer_values();
    test_digamma_complex();
    test_digamma_large_argument();
    std::printf("All digamma tests passed!\n");
    return 0;
}
