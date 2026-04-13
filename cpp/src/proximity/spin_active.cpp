// Spin-active interface corrections.
// Modifies the effective kernel to account for
// spin-dependent interface scattering (Eschrig PRB 80, 2009).
//
// K_eff = K * cos(θ/2) * (1 - P * sin(θ/2))
//
// where θ is the spin-mixing angle and P is the interface polarization.

#include "supermag/proximity.h"
#include <cmath>
#include <complex>

namespace supermag {

std::complex<double> apply_spin_active(
    const supermag_spin_active_t *sa,
    std::complex<double> K)
{
    if (!sa) return K;

    double theta = sa->mixing_angle;
    double P = sa->polarization;
    if (P < 0.0) P = 0.0;
    if (P > 1.0) P = 1.0;

    // K_eff = K * cos(θ/2) * (1 - P * sin(θ/2))
    double half_theta = theta / 2.0;
    double factor = std::cos(half_theta) * (1.0 - P * std::sin(half_theta));

    return K * factor;
}

} // namespace supermag
