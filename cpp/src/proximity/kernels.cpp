// Proximity effect kernels: coth and tanh
// Internal to the proximity solver — not exposed in public C API.
//
// q = (1+i) / xi_F  is the complex inverse decay length in the ferromagnet.
// coth kernel: K = q * coth(q * d_F)  — 0-junction (S/F bilayer with vacuum boundary)
// tanh kernel: K = q * tanh(q * d_F)  — π-junction (semi-infinite/thick-F limit)
//
// Overflow safety: for |Re(q*d_F)| > 350 the raw sinh/cosh would overflow
// double precision.  In that regime coth(z) → sgn(Re(z)) and tanh(z) → sgn(Re(z)),
// so the kernel reduces to q * (±1).

#include <complex>
#include <cmath>

namespace supermag {

// Original overflow threshold (raw sinh/cosh overflow).
static constexpr double OVERFLOW_THRESHOLD = 350.0;

// Tighter safe threshold: coth(z) − 1 ≈ 2·exp(−2Re(z)) < 2e-3 for Re(z) > 3.0.
static constexpr double SAFE_THRESHOLD = 3.0;

static std::complex<double> compute_q(double xi_F) {
    return std::complex<double>(1.0, 1.0) / xi_F;
}

std::complex<double> kernel_coth(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;

    // Overflow-safe: for large |Re(qd)|, coth(z) → sgn(Re(z))
    double re = qd.real();
    if (std::abs(re) > OVERFLOW_THRESHOLD) {
        double sign = (re > 0.0) ? 1.0 : -1.0;
        return q * sign;
    }

    // coth(z) = cosh(z)/sinh(z)
    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (ch / sh);
}

std::complex<double> kernel_tanh(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;

    // Overflow-safe: for large |Re(qd)|, tanh(z) → sgn(Re(z))
    double re = qd.real();
    if (std::abs(re) > OVERFLOW_THRESHOLD) {
        double sign = (re > 0.0) ? 1.0 : -1.0;
        return q * sign;
    }

    // tanh(z) = sinh(z)/cosh(z)
    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (sh / ch);
}

// Safe variants with tighter threshold (Re(qd) > 3.0).
// Uses asymptotic form much earlier, avoiding loss-of-precision
// in the coth/tanh ratio for moderately large arguments.

std::complex<double> kernel_coth_safe(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;

    double re = qd.real();
    if (std::abs(re) > SAFE_THRESHOLD) {
        double sign = (re > 0.0) ? 1.0 : -1.0;
        return q * sign;
    }

    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (ch / sh);
}

std::complex<double> kernel_tanh_safe(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;

    double re = qd.real();
    if (std::abs(re) > SAFE_THRESHOLD) {
        double sign = (re > 0.0) ? 1.0 : -1.0;
        return q * sign;
    }

    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (sh / ch);
}

} // namespace supermag
