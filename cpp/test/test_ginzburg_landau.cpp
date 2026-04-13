// Unit tests for Ginzburg-Landau solver
#include "supermag/ginzburg_landau.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_gl_uniform_equilibrium() {
    // α < 0: equilibrium |ψ|² = -α/β
    const int nx = 20, ny = 20;
    int N = nx * ny;
    std::vector<double> psi_r(N, 0.0), psi_i(N, 0.0);

    double alpha = -1.0, beta = 1.0, kappa = 1.0, dx = 0.5;
    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                  SUPERMAG_GL_SCALAR, 0.0,
                                  psi_r.data(), psi_i.data());
    assert(rc == SUPERMAG_OK);

    // Check |ψ|² ≈ -α/β = 1.0 at center
    double expected = -alpha / beta;  // 1.0
    int center = (ny/2) * nx + nx/2;
    double abs2 = psi_r[center] * psi_r[center] + psi_i[center] * psi_i[center];
    assert(std::fabs(abs2 - expected) < 0.2);
    std::printf("  PASS: test_gl_uniform_equilibrium (|psi|^2=%.4f, expected=%.4f)\n",
                abs2, expected);
}

void test_gl_normal_state() {
    // α > 0: ψ → 0 (normal state)
    const int nx = 16, ny = 16;
    int N = nx * ny;
    std::vector<double> psi_r(N, 0.0), psi_i(N, 0.0);

    double alpha = 1.0, beta = 1.0, kappa = 1.0, dx = 0.5;
    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                  SUPERMAG_GL_SCALAR, 0.0,
                                  psi_r.data(), psi_i.data());
    assert(rc == SUPERMAG_OK);

    // ψ should be near zero everywhere
    double max_abs2 = 0.0;
    for (int i = 0; i < N; ++i) {
        double abs2 = psi_r[i] * psi_r[i] + psi_i[i] * psi_i[i];
        if (abs2 > max_abs2) max_abs2 = abs2;
    }
    assert(max_abs2 < 0.1);
    std::printf("  PASS: test_gl_normal_state (max|psi|^2=%.6f)\n", max_abs2);
}

void test_gl_kappa_effect() {
    // Varying κ should change the coherence length and thus the spatial profile
    const int nx = 20, ny = 20;
    int N = nx * ny;
    double alpha = -1.0, beta = 1.0, dx = 0.5;

    std::vector<double> psi_r1(N, 0.0), psi_i1(N, 0.0);
    std::vector<double> psi_r2(N, 0.0), psi_i2(N, 0.0);

    supermag_gl_minimize(alpha, beta, 0.5, nx, ny, dx,
                         SUPERMAG_GL_SCALAR, 0.0, psi_r1.data(), psi_i1.data());
    supermag_gl_minimize(alpha, beta, 2.0, nx, ny, dx,
                         SUPERMAG_GL_SCALAR, 0.0, psi_r2.data(), psi_i2.data());

    // Both should converge to something non-trivial
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < N; ++i) {
        sum1 += psi_r1[i] * psi_r1[i] + psi_i1[i] * psi_i1[i];
        sum2 += psi_r2[i] * psi_r2[i] + psi_i2[i] * psi_i2[i];
    }
    assert(sum1 > 0.0);
    assert(sum2 > 0.0);
    std::printf("  PASS: test_gl_kappa_effect (sum1=%.2f, sum2=%.2f)\n", sum1/N, sum2/N);
}

void test_gl_preserve_ic() {
    // If psi_real[0] != 0, should preserve user initial condition
    const int nx = 10, ny = 10;
    int N = nx * ny;
    std::vector<double> psi_r(N, 0.0), psi_i(N, 0.0);

    // Set a custom IC
    for (int i = 0; i < N; ++i) {
        psi_r[i] = 0.8;
        psi_i[i] = 0.0;
    }

    double alpha = -1.0, beta = 1.0, kappa = 1.0, dx = 0.5;
    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                  SUPERMAG_GL_SCALAR, 0.0,
                                  psi_r.data(), psi_i.data());
    assert(rc == SUPERMAG_OK);

    // Should have converged to equilibrium from user IC
    double abs2 = psi_r[N/2] * psi_r[N/2] + psi_i[N/2] * psi_i[N/2];
    assert(abs2 > 0.5);  // should be near equilibrium
    std::printf("  PASS: test_gl_preserve_ic (|psi|^2=%.4f)\n", abs2);
}

void test_gl_with_field() {
    // Test extended solver with applied field
    const int nx = 16, ny = 16;
    int N = nx * ny;
    std::vector<double> psi_r(N, 0.0), psi_i(N, 0.0);

    double alpha = -1.0, beta = 1.0, kappa = 2.0, dx = 0.5;
    double H = 0.1;

    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                      SUPERMAG_GL_GAUGE, H,
                                      psi_r.data(), psi_i.data());
    assert(rc == SUPERMAG_OK);

    // With field, ψ should still be finite
    bool has_nonzero = false;
    for (int i = 0; i < N; ++i) {
        double abs2 = psi_r[i] * psi_r[i] + psi_i[i] * psi_i[i];
        assert(std::isfinite(abs2));
        if (abs2 > 0.01) has_nonzero = true;
    }
    assert(has_nonzero);
    std::printf("  PASS: test_gl_with_field\n");
}

void test_gl_null() {
    int rc = supermag_gl_minimize(-1.0, 1.0, 1.0, 10, 10, 0.5,
                                  SUPERMAG_GL_SCALAR, 0.0, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_gl_null\n");
}

int main() {
    std::printf("Running Ginzburg-Landau tests...\n");
    test_gl_uniform_equilibrium();
    test_gl_normal_state();
    test_gl_kappa_effect();
    test_gl_preserve_ic();
    test_gl_with_field();
    test_gl_null();
    std::printf("All Ginzburg-Landau tests passed!\n");
    return 0;
}
