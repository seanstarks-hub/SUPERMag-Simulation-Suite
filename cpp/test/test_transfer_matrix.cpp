// Unit tests for the transfer matrix utility.
#include "../src/proximity/transfer_matrix.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <complex>

// Forward-declare internal kernel functions for cross-validation
namespace supermag {
std::complex<double> kernel_coth(double d_F, double xi_F);
std::complex<double> kernel_tanh(double d_F, double xi_F);
}

void test_identity() {
    auto I = supermag::mat_identity();
    assert(std::abs(I.m[0][0] - 1.0) < 1e-15);
    assert(std::abs(I.m[0][1]) < 1e-15);
    assert(std::abs(I.m[1][0]) < 1e-15);
    assert(std::abs(I.m[1][1] - 1.0) < 1e-15);

    // Zero-thickness layer should give identity
    auto q = std::complex<double>(1.0, 1.0);
    auto M = supermag::layer_transfer_matrix(q, 0.0);
    assert(std::abs(M.m[0][0] - 1.0) < 1e-15);
    assert(std::abs(M.m[0][1]) < 1e-15);
    assert(std::abs(M.m[1][0]) < 1e-15);
    assert(std::abs(M.m[1][1] - 1.0) < 1e-15);
    std::printf("  PASS: test_identity\n");
}

void test_determinant() {
    // det(M) = cosh²(qd) - sinh²(qd) = 1 for any q, d
    double xi_F = 5.0;
    auto q = std::complex<double>(1.0, 1.0) / xi_F;

    for (double d : {1.0, 5.0, 10.0, 50.0, 100.0}) {
        auto M = supermag::layer_transfer_matrix(q, d);
        auto det = supermag::mat_det(M);
        double err = std::abs(det - 1.0);
        // Tolerance scales with matrix-entry magnitudes due to catastrophic
        // cancellation in cosh²-sinh² for large |Re(qd)|.
        double scale = std::abs(M.m[0][0]) * std::abs(M.m[1][1])
                     + std::abs(M.m[0][1]) * std::abs(M.m[1][0]);
        double tol = std::max(1e-12, scale * 1e-12);
        assert(err < tol);
    }
    std::printf("  PASS: test_determinant\n");
}

void test_single_layer_coth() {
    // For a single layer with vacuum boundary condition (K_boundary = ∞),
    // we can recover kernel_coth by constructing the transfer matrix and
    // extracting K = M[1][0] / M[0][0] = q * sinh(qd) / cosh(qd) ... wait,
    // that gives tanh. The impedance convention:
    //
    // For coth kernel (pi-junction):
    //   K_coth = q * coth(qd) = K extracted from M with free-surface BC
    //
    // M[1][0] / M[0][0] = q*sinh(qd) / cosh(qd) = q*tanh(qd) = K_tanh
    // M[0][0] / M[0][1] = cosh(qd) / (sinh(qd)/q) = q*coth(qd) = K_coth
    //
    // So coth kernel = M[0][0] * q / M[1][0]... no.  Let's verify:
    // The correct relation depends on boundary conditions.
    // For a free surface: K = M[0][0]/M[0][1] = q*coth(qd)
    //
    // Actually, with our convention:
    //   M = [[cosh(qd), sinh(qd)/q], [q*sinh(qd), cosh(qd)]]
    //   extract_kernel = M[1][0]/M[0][0] = q*sinh(qd)/cosh(qd) = q*tanh(qd)
    //
    // So extract_kernel for a single layer gives the tanh kernel.
    // The coth kernel needs M[0][0]/M[0][1] = cosh(qd)/(sinh(qd)/q) = q*coth(qd)

    double d_F = 10.0, xi_F = 5.0;
    auto q = std::complex<double>(1.0, 1.0) / xi_F;
    auto M = supermag::layer_transfer_matrix(q, d_F);

    // Verify: M[0][0]/M[0][1] gives coth kernel
    auto K_coth_from_M = M.m[0][0] * q / (M.m[1][0]); // no, use the ratio directly
    // M[0][0]/M[0][1] = cosh/(sinh/q) = q*cosh/sinh = q*coth
    auto K_coth_ratio = M.m[0][0] / M.m[0][1];   // = cosh(qd) / (sinh(qd)/q) = q*coth(qd)

    auto K_coth_ref = supermag::kernel_coth(d_F, xi_F);
    assert(std::abs(K_coth_ratio - K_coth_ref) < 1e-12);

    std::printf("  PASS: test_single_layer_coth\n");
}

void test_single_layer_tanh() {
    // extract_kernel(M) = M[1][0]/M[0][0] = q*tanh(qd) = kernel_tanh
    double d_F = 10.0, xi_F = 5.0;
    auto M = supermag::layer_transfer_matrix(
        std::complex<double>(1.0, 1.0) / xi_F, d_F);

    auto K_tanh_from_M = supermag::extract_kernel(M);
    auto K_tanh_ref = supermag::kernel_tanh(d_F, xi_F);
    assert(std::abs(K_tanh_from_M - K_tanh_ref) < 1e-12);
    std::printf("  PASS: test_single_layer_tanh\n");
}

void test_composition() {
    // Two identical layers of thickness d should give same result
    // as one layer of thickness 2d (for transfer matrices).
    double d = 5.0, xi_F = 3.0;
    auto q = std::complex<double>(1.0, 1.0) / xi_F;

    auto M_single = supermag::layer_transfer_matrix(q, d);
    auto M_double = supermag::layer_transfer_matrix(q, 2.0 * d);
    auto M_composed = supermag::mat_multiply(M_single, M_single);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            assert(std::abs(M_composed.m[i][j] - M_double.m[i][j]) < 1e-10);

    std::printf("  PASS: test_composition\n");
}

void test_overflow_transfer_matrix() {
    // Very thick layer should not produce NaN/Inf
    auto q = std::complex<double>(1.0, 1.0) / 1.0;  // xi_F = 1
    auto M = supermag::layer_transfer_matrix(q, 2000.0);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            assert(std::isfinite(M.m[i][j].real()) && std::isfinite(M.m[i][j].imag()));

    auto K = supermag::extract_kernel(M);
    assert(std::isfinite(K.real()) && std::isfinite(K.imag()));
    std::printf("  PASS: test_overflow_transfer_matrix\n");
}

int main() {
    std::printf("Running transfer matrix tests...\n");
    test_identity();
    test_determinant();
    test_single_layer_coth();
    test_single_layer_tanh();
    test_composition();
    test_overflow_transfer_matrix();
    std::printf("All transfer matrix tests passed!\n");
    return 0;
}
