// Magnetic domain effective proximity kernel.                     [EQ-15]
//
// For a ferromagnet with alternating magnetic domains (±E_ex), the
// layer is divided into n_domains equal domains of width d_domain.
// Adjacent domains have opposite magnetization direction, which flips
// the sign of the imaginary part of the wave vector q:
//   domain with +E_ex:  q = (1+i)/xi_F
//   domain with −E_ex:  q = (1−i)/xi_F
//
// Domain walls: if domain_wall > 0, a smooth transition region is
// inserted between each pair of domains. The wall is discretized into
// 10 sub-slices with linearly interpolated exchange energy.
//
// Transfer matrices are cascaded from vacuum to S-side.
//
// Limits:
//   n_domains = 1:  matches standard bilayer kernel
//   domain_wall = 0:  sharp domain boundaries

#include "supermag/proximity.h"
#include "transfer_matrix.h"
#include <complex>
#include <cmath>

// Number of sub-slices used to resolve each domain wall
static constexpr int WALL_SLICES = 10;

extern "C" {

int supermag_proximity_kernel_domains(
    double d_F, double xi_F, double E_ex,
    const supermag_domain_params_t *dom,
    supermag_phase_t phase,
    double *K_real, double *K_imag)
{
    if (!dom || !K_real || !K_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (d_F <= 0.0 || xi_F <= 0.0 || E_ex <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;
    if (dom->domain_width <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Infer number of domains from total thickness and domain width
    int N = static_cast<int>(std::round(d_F / dom->domain_width));
    if (N < 1) N = 1;
    double d_dom = dom->domain_width;
    double d_wall = (dom->domain_wall > 0.0) ? dom->domain_wall : 0.0;

    // Wave vectors for +E_ex and −E_ex domains
    // +E_ex: q = (1+i)/xi_F  (standard)
    // −E_ex: q = (1−i)/xi_F  (flipped imaginary part)
    std::complex<double> q_plus(1.0, 1.0);
    q_plus /= xi_F;
    std::complex<double> q_minus(1.0, -1.0);
    q_minus /= xi_F;

    // Build the total transfer matrix from vacuum (domain N) to S-side (domain 1).
    auto M_total = supermag::mat_identity();

    for (int i = N - 1; i >= 0; --i) {
        // Alternating sign: domain 0 = +E_ex, domain 1 = −E_ex, ...
        bool positive = (i % 2 == 0);
        auto q_dom = positive ? q_plus : q_minus;

        // Effective domain thickness: reduce by half-walls on each side
        double d_eff = d_dom;
        if (d_wall > 0.0 && N > 1) {
            if (i == 0 || i == N - 1)
                d_eff -= d_wall / 2.0;  // edge domain: one wall side
            else
                d_eff -= d_wall;        // interior domain: half-wall each side
            if (d_eff < 0.0) d_eff = 0.0;
        }

        // Domain transfer matrix
        if (d_eff > 0.0) {
            auto M_dom = supermag::layer_transfer_matrix(q_dom, d_eff);
            M_total = supermag::mat_multiply(M_dom, M_total);
        }

        // Domain wall between domain i and domain i-1 (towards S-side)
        // The wall interpolates E_ex linearly from one domain to the next.
        if (d_wall > 0.0 && i > 0) {
            bool next_positive = ((i - 1) % 2 == 0);
            auto q_start = positive ? q_plus : q_minus;
            auto q_end   = next_positive ? q_plus : q_minus;

            double slice_thickness = d_wall / WALL_SLICES;
            for (int s = 0; s < WALL_SLICES; ++s) {
                // Interpolation parameter: 0 at domain i, 1 at domain i-1
                double t = (s + 0.5) / WALL_SLICES;
                // Linear interpolation of q between domains
                auto q_wall = q_start * (1.0 - t) + q_end * t;
                auto M_slice = supermag::layer_transfer_matrix(q_wall, slice_thickness);
                M_total = supermag::mat_multiply(M_slice, M_total);
            }
        }
    }

    // Extract kernel based on phase (boundary condition)
    std::complex<double> K;
    if (phase == SUPERMAG_PHASE_ZERO) {
        // coth kernel: K = M[1][0] / M[0][0]  (0-junction)
        K = supermag::extract_kernel(M_total);
    } else if (phase == SUPERMAG_PHASE_PI) {
        // tanh kernel: K = M[0][0] / M[0][1]  (π-junction)
        if (std::abs(M_total.m[0][1]) < 1e-300) {
            *K_real = 0.0;
            *K_imag = 0.0;
            return SUPERMAG_OK;
        }
        K = M_total.m[0][0] / M_total.m[0][1];
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    *K_real = K.real();
    *K_imag = K.imag();
    return SUPERMAG_OK;
}

}
