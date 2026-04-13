// S/N/F trilayer effective proximity kernel.                       [EQ-13]
//
// Composes the F-layer kernel with N-layer (normal metal) propagation
// using the continued-fraction formula:
//
//   K_SNF = q_N · [ K_F + q_N·tanh(q_N·d_N) ]
//                / [ q_N + K_F·tanh(q_N·d_N) ]
//
// where K_F is the standard bilayer kernel (coth or tanh),
// and q_N = 1/xi_N is real (no exchange splitting in N layer).
//
// Limits:
//   d_N → 0:  K_SNF → K_F  (recovers bilayer)
//   d_N → ∞:  K_SNF → q_N  (S layer sees only normal metal)

#include "supermag/proximity.h"
#include <complex>
#include <cmath>

namespace supermag {
// Forward declarations from kernels.cpp
std::complex<double> kernel_coth(double d_F, double xi_F);
std::complex<double> kernel_tanh(double d_F, double xi_F);
}

// Overflow-safe tanh for real argument (used for N-layer propagation)
static double safe_tanh(double x) {
    if (std::abs(x) > 350.0)
        return (x > 0.0) ? 1.0 : -1.0;
    return std::tanh(x);
}

extern "C" {

int supermag_proximity_kernel_snf(
    double d_F, double xi_F,
    const supermag_trilayer_params_t *tri,
    supermag_phase_t phase,
    double *K_real, double *K_imag)
{
    if (!tri || !K_real || !K_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (d_F <= 0.0 || xi_F <= 0.0 || tri->xi_N <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Compute F-layer kernel (bilayer), already overflow-safe
    std::complex<double> K_F;
    if (phase == SUPERMAG_PHASE_ZERO) {
        K_F = supermag::kernel_coth(d_F, xi_F);
    } else if (phase == SUPERMAG_PHASE_PI) {
        K_F = supermag::kernel_tanh(d_F, xi_F);
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    double d_N = tri->d_N;

    // Edge case: zero N-layer thickness → bilayer
    if (d_N <= 0.0) {
        *K_real = K_F.real();
        *K_imag = K_F.imag();
        return SUPERMAG_OK;
    }

    // N-layer wave vector: real, no exchange splitting
    double q_N = 1.0 / tri->xi_N;
    double tanh_qN_dN = safe_tanh(q_N * d_N);

    // Continued-fraction composition: [EQ-13]
    //   K_SNF = q_N · (K_F + q_N·tanh(q_N·d_N)) / (q_N + K_F·tanh(q_N·d_N))
    std::complex<double> qN(q_N, 0.0);
    std::complex<double> numerator   = K_F + qN * tanh_qN_dN;
    std::complex<double> denominator = qN + K_F * tanh_qN_dN;

    if (std::abs(denominator) < 1e-300) {
        *K_real = 0.0;
        *K_imag = 0.0;
        return SUPERMAG_OK;
    }

    std::complex<double> K_SNF = qN * numerator / denominator;

    *K_real = K_SNF.real();
    *K_imag = K_SNF.imag();
    return SUPERMAG_OK;
}

}
