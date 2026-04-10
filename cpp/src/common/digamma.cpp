#include "digamma.h"
#include <cmath>

namespace supermag {

// Bernoulli numbers B_{2k} / (2k) for k=1..8
// used in the asymptotic expansion of the digamma function:
//   psi(z) ~ ln(z) - 1/(2z) - sum_{k=1}^{N} B_{2k}/(2k * z^{2k})
static const double B2k_over_2k[] = {
    1.0 / 12.0,       // B2/(2)   = 1/6 / 2
    1.0 / 120.0,      // B4/(4)   = -1/30 / 4  ... but sign handled below
    1.0 / 252.0,      // B6/(6)
    1.0 / 240.0,      // B8/(8)
    1.0 / 132.0,      // B10/(10)
    691.0 / 32760.0,   // B12/(12)
    1.0 / 12.0,        // B14/(14)
    3617.0 / 8160.0,   // B16/(16)
};

// Full Bernoulli numbers B_{2k} for k=1..8
// B2=1/6, B4=-1/30, B6=1/42, B8=-1/30, B10=5/66, B12=-691/2730, B14=7/6, B16=-3617/510
// The asymptotic series is: psi(z) = ln(z) - 1/(2z) - sum_{k=1}^N B_{2k}/(2k * z^{2k})
static const double bernoulli_2k[] = {
     1.0 / 6.0,         // B2
    -1.0 / 30.0,        // B4
     1.0 / 42.0,        // B6
    -1.0 / 30.0,        // B8
     5.0 / 66.0,        // B10
    -691.0 / 2730.0,    // B12
     7.0 / 6.0,         // B14
    -3617.0 / 510.0,    // B16
};

std::complex<double> digamma(std::complex<double> z) {
    // Use recurrence: psi(z) = psi(z+1) - 1/z
    // Shift z until |z| > 10 for asymptotic convergence
    std::complex<double> result(0.0, 0.0);

    while (std::abs(z) < 10.0) {
        result -= 1.0 / z;
        z += 1.0;
    }

    // Asymptotic expansion:
    // psi(z) = ln(z) - 1/(2z) - sum_{k=1}^{8} B_{2k} / (2k * z^{2k})
    result += std::log(z) - 0.5 / z;

    std::complex<double> z2 = z * z;
    std::complex<double> z_power = z2; // z^{2k} starts at z^2

    for (int k = 0; k < 8; ++k) {
        result -= bernoulli_2k[k] / (2.0 * (k + 1) * z_power);
        z_power *= z2;
    }

    return result;
}

} // namespace supermag
