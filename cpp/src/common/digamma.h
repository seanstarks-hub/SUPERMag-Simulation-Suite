#ifndef SUPERMAG_DIGAMMA_H
#define SUPERMAG_DIGAMMA_H

#include <complex>

namespace supermag {

// Digamma (psi) function for complex argument z.
// Uses recurrence relation to shift |z| > 10, then asymptotic series.
std::complex<double> digamma(std::complex<double> z);

} // namespace supermag

#endif
