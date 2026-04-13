// Transfer matrix utility — 2×2 complex matrix operations for
// cascaded proximity layer composition.
//
// See transfer_matrix.h for interface documentation.

#include "transfer_matrix.h"
#include <cmath>

namespace supermag {

// Threshold matching kernels.cpp — beyond this, sinh/cosh overflow.
static constexpr double OVERFLOW_THRESHOLD = 350.0;

Matrix2x2 mat_identity() {
    Matrix2x2 I;
    I.m[0][0] = 1.0; I.m[0][1] = 0.0;
    I.m[1][0] = 0.0; I.m[1][1] = 1.0;
    return I;
}

Matrix2x2 mat_multiply(const Matrix2x2 &A, const Matrix2x2 &B) {
    Matrix2x2 C;
    C.m[0][0] = A.m[0][0]*B.m[0][0] + A.m[0][1]*B.m[1][0];
    C.m[0][1] = A.m[0][0]*B.m[0][1] + A.m[0][1]*B.m[1][1];
    C.m[1][0] = A.m[1][0]*B.m[0][0] + A.m[1][1]*B.m[1][0];
    C.m[1][1] = A.m[1][0]*B.m[0][1] + A.m[1][1]*B.m[1][1];
    return C;
}

std::complex<double> mat_det(const Matrix2x2 &M) {
    return M.m[0][0]*M.m[1][1] - M.m[0][1]*M.m[1][0];
}

Matrix2x2 layer_transfer_matrix(std::complex<double> q, double delta) {
    Matrix2x2 M;

    if (delta <= 0.0 || std::abs(q) < 1e-30) {
        return mat_identity();
    }

    auto qd = q * delta;
    double re = qd.real();

    std::complex<double> sh, ch;
    if (std::abs(re) > OVERFLOW_THRESHOLD) {
        // Asymptotic: for large |Re(z)|,
        //   cosh(z) ≈ sinh(z) ≈ (sgn(Re(z))/2) * exp(|Re(z)|) * exp(i*Im(z))
        // The exponential factor cancels in ratios, so we track only the
        // relative phase.  We normalize so that cosh = 1, and
        // sinh/cosh = tanh → sgn(Re(z)), giving:
        //   M = [[1,       1/q * sgn],
        //        [q * sgn, 1        ]]
        // This is the thick-layer limit where the transfer matrix
        // simply couples impedances.
        double sign = (re > 0.0) ? 1.0 : -1.0;
        std::complex<double> s(sign, 0.0);
        M.m[0][0] = std::complex<double>(1.0, 0.0);
        M.m[0][1] = s / q;
        M.m[1][0] = q * s;
        M.m[1][1] = std::complex<double>(1.0, 0.0);
    } else {
        sh = std::sinh(qd);
        ch = std::cosh(qd);
        M.m[0][0] = ch;
        M.m[0][1] = sh / q;
        M.m[1][0] = q * sh;
        M.m[1][1] = ch;
    }

    return M;
}

std::complex<double> extract_kernel(const Matrix2x2 &M_total) {
    // K_eff = M[1][0] / M[0][0]
    // Guard against zero denominator (pathological case)
    if (std::abs(M_total.m[0][0]) < 1e-300) {
        return std::complex<double>(0.0, 0.0);
    }
    return M_total.m[1][0] / M_total.m[0][0];
}

} // namespace supermag
