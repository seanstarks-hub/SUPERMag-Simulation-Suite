// Graded ferromagnet effective proximity kernel.                   [EQ-14]
//
// For a ferromagnet with position-dependent exchange energy E_ex(x),
// the layer is sliced into n_slices sub-layers. Each slice has a local
// coherence length:
//   xi_F_local = xi_F_ref * sqrt(E_ex_ref / E_ex_local)
//
// Transfer matrices are cascaded from the vacuum boundary inward:
//   M_total = M_1 · M_2 · ... · M_N   (1 = S-side, N = vacuum-side)
//
// Exchange profiles (x ∈ [0, 1] normalized position within F layer):
//   LINEAR:      E(x) = E_left + (E_right - E_left) * x
//   EXPONENTIAL: E(x) = E_left * exp(ln(E_right/E_left) * x)
//   STEP:        E(x) = E_left for x < 0.5, E_right for x >= 0.5
//
// Limits:
//   E_ex_left == E_ex_right → uniform F, should match bilayer kernel
//   n_slices → ∞ → continuous grading (converges)

#include "supermag/proximity.h"
#include "transfer_matrix.h"
#include <complex>
#include <cmath>

// Compute local exchange energy at normalized position t ∈ [0,1]
// (0 = S/F interface, 1 = vacuum/free surface)
static double exchange_profile(double E_surface, double E_bulk,
                               supermag_grade_profile_t profile, double t) {
    switch (profile) {
    case SUPERMAG_GRADE_LINEAR:
        return E_surface + (E_bulk - E_surface) * t;
    case SUPERMAG_GRADE_EXPONENTIAL:
        if (E_surface <= 0.0 || E_bulk <= 0.0)
            return E_surface;  // fallback
        return E_surface * std::exp(std::log(E_bulk / E_surface) * t);
    case SUPERMAG_GRADE_STEP:
        return (t < 0.5) ? E_surface : E_bulk;
    default:
        return E_surface;
    }
}

extern "C" {

int supermag_proximity_kernel_graded(
    double d_F, double xi_F_ref,
    const supermag_graded_params_t *grade,
    supermag_phase_t phase,
    double *K_real, double *K_imag)
{
    if (!grade || !K_real || !K_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (d_F <= 0.0 || xi_F_ref <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;
    if (grade->n_slabs < 1)
        return SUPERMAG_ERR_INVALID_DIM;
    if (grade->E_ex_surface <= 0.0 || grade->E_ex_bulk <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    int N = grade->n_slabs;
    double delta = d_F / N;  // thickness of each slice

    // Build total transfer matrix from vacuum side (slice N) to S side (slice 1).
    // Convention: slice index i goes from N-1 (vacuum) to 0 (S-side).
    // We accumulate M_total = M_0 * M_1 * ... * M_{N-1}
    // where M_0 is the S-side slice and M_{N-1} is the vacuum-side slice.
    //
    // For a free surface (coth kernel / 0-junction), we cascade and then
    // extract K = M[0][0]/M[0][1] = q_coth behavior.
    // For tanh kernel (pi-junction), K = M[1][0]/M[0][0].
    //
    // We compose right-to-left: start from vacuum boundary.

    auto M_total = supermag::mat_identity();
    double E_ref = grade->E_ex_surface;

    for (int i = N - 1; i >= 0; --i) {
        // Normalized position at center of slab
        double t = (static_cast<double>(i) + 0.5) / N;
        double E_local = exchange_profile(grade->E_ex_surface, grade->E_ex_bulk,
                                          grade->profile, t);

        // Local coherence length: xi_F_local = xi_F_ref * sqrt(E_ref / E_local)
        double xi_local = xi_F_ref * std::sqrt(E_ref / E_local);
        auto q_local = std::complex<double>(1.0, 1.0) / xi_local;

        auto M_slice = supermag::layer_transfer_matrix(q_local, delta);
        M_total = supermag::mat_multiply(M_slice, M_total);
    }

    // Extract kernel based on phase (boundary condition)
    std::complex<double> K;
    if (phase == SUPERMAG_PHASE_ZERO) {
        // coth kernel: K = M[0][0] / M[0][1]  (free surface BC)
        if (std::abs(M_total.m[0][1]) < 1e-300) {
            *K_real = 0.0;
            *K_imag = 0.0;
            return SUPERMAG_OK;
        }
        K = M_total.m[0][0] / M_total.m[0][1];
    } else if (phase == SUPERMAG_PHASE_PI) {
        // tanh kernel: K = M[1][0] / M[0][0]  (fixed surface BC)
        K = supermag::extract_kernel(M_total);
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    *K_real = K.real();
    *K_imag = K.imag();
    return SUPERMAG_OK;
}

}
