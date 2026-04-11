// Proximity effect kernels: coth and tanh
// Internal to the proximity solver — not exposed in public C API.
//
// q = (1+i) / xi_F  is the complex inverse decay length in the ferromagnet.
// coth kernel: K = q * coth(q * d_F)  — 0-junction (S/F bilayer with vacuum boundary)
// tanh kernel: K = q * tanh(q * d_F)  — π-junction (semi-infinite/thick-F limit)

#include <complex>
#include <cmath>

namespace supermag {

static std::complex<double> compute_q(double xi_F) {
    return std::complex<double>(1.0, 1.0) / xi_F;
}

std::complex<double> kernel_coth(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;
    // coth(z) = cosh(z)/sinh(z)
    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (ch / sh);
}

std::complex<double> kernel_tanh(double d_F, double xi_F) {
    auto q = compute_q(xi_F);
    auto qd = q * d_F;
    // tanh(z) = sinh(z)/cosh(z)
    auto sh = std::sinh(qd);
    auto ch = std::cosh(qd);
    return q * (sh / ch);
}

} // namespace supermag
