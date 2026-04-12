// Proximity effect Tc solver — digamma-based self-consistency equations.
//
// Thin-S model:
//   ln(Tc0/T) = Re[ psi(1/2 + gamma*K/(2*pi*T) + Lambda_dep) - psi(1/2) ]
//   where K is the kernel (coth or tanh), gamma is interface transparency.
//   Find highest root F(T) = 0 via Brent's method.
//
// Fominov model (PRB 66, 014507):
//   2x2 matrix M(T) encodes S-layer thickness, gamma, gamma_B, kernel.
//   Tc is where det M(T) = 0.
//   Uses the same scan + Brent strategy on the determinant.

#define _USE_MATH_DEFINES
#include "supermag/proximity.h"
#include "../common/digamma.h"
#include "../solvers/root_scalar.h"
#include <cmath>
#include <complex>
#include <cstring>

namespace supermag {
// Forward declarations from kernels.cpp
std::complex<double> kernel_coth(double d_F, double xi_F);
std::complex<double> kernel_tanh(double d_F, double xi_F);
}

// Context struct passed to the root-finding callback
struct ThinSContext {
    double Tc0;
    double gamma;
    double lambda_dep;
    std::complex<double> K;  // kernel value at this d_F
};

// F(T) = ln(Tc0/T) - Re[ psi(1/2 + A(T)) - psi(1/2) ]
// where A(T) = gamma * K / (2*pi*T / Tc0) + lambda_dep
// We use reduced temperature t = T/Tc0, so the Matsubara argument is scaled.
//
// More precisely, the self-consistency condition is:
//   ln(Tc0/T) = Re[ psi(1/2 + Gamma_N * K_F / (2*pi*kB*T)) - psi(1/2) ]
//
// In the thin-S limit with dimensionless coupling:
//   A = gamma * xi_F * K / (2 * pi * T / Tc0)  but we absorb scales into gamma.
//
// Simplified thin-S model:
//   F(T) = ln(Tc0/T) - Re[ psi(1/2 + gamma*K_norm + lambda_dep) - psi(1/2) ]
// where K_norm = K * xi_F (dimensionless kernel), and gamma absorbs coupling strength.
//
// For the Nb/CuNi system (Fominov fit), the relevant formula is:
//   ln(Tc0/Tc) = Re[ psi(1/2 + alpha(Tc)) - psi(1/2) ]
//   alpha(Tc) = gamma / (gamma_B + K_F^{-1} * q_S * cot(q_S * d_S))
// where q_S = sqrt(2*pi*T/(D_S)), K_F = kernel.

static double thin_s_equation(double T, void *ctx) {
    auto *c = static_cast<ThinSContext*>(ctx);
    if (T <= 0.0) return -1.0;

    double log_ratio = std::log(c->Tc0 / T);

    // Complex pair-breaking parameter: A = gamma * K * Tc0 / (2*pi*T)  [EQ-4]
    // No eta = xi_S/d_S prefactor — gamma absorbs coupling strength.
    std::complex<double> half(0.5, 0.0);
    std::complex<double> A_complex = c->gamma * c->K * c->Tc0 / (2.0 * M_PI * T);
    A_complex += c->lambda_dep;

    auto psi_arg = supermag::digamma(half + A_complex);
    auto psi_half = supermag::digamma(half);

    return log_ratio - (psi_arg - psi_half).real();
}

// Fominov model context
struct FominovContext {
    double Tc0;
    double d_S;
    double xi_S;
    double gamma;
    double gamma_B;
    double lambda_dep;
    std::complex<double> K;  // kernel at this d_F
};

// det M(T) for the Fominov model.
// M is a 2x2 matrix encoding S and F layer self-consistency:
//
// The Fominov (PRB 66 014507) Tc equation for an S/F bilayer:
//   ln(Tc0/Tc) = Re[ psi(1/2 + alpha) - psi(1/2) ]
// where
//   alpha = gamma * K / (1 + gamma_B * K + Omega_S(T))
//   Omega_S(T) = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S/xi_S)
//
// The S-layer thermal impedance Omega_S(T) captures the T-dependent
// pair-breaking from finite S-layer thickness. It reduces to:
//   - Omega_S → 0 as T → 0 (recovers simplified formula)
//   - Omega_S → xi_S/d_S in the thin-S limit
//   - Omega_S → sqrt(T/Tc0) in the thick-S limit (d_S >> xi_S)
static double fominov_determinant(double T, void *ctx) {
    auto *c = static_cast<FominovContext*>(ctx);
    if (T <= 0.0) return -1.0;

    double log_ratio = std::log(c->Tc0 / T);

    // S-layer thermal impedance  [EQ-5]
    // Omega_S(T) = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S / xi_S)
    double sqrt_t_ratio = std::sqrt(T / c->Tc0);
    double qS_dS = sqrt_t_ratio * c->d_S / c->xi_S;
    double Omega_S = 0.0;
    if (c->xi_S > 0.0 && c->d_S > 0.0) {
        // coth(x) = cosh(x)/sinh(x); for large x, coth → 1
        if (qS_dS > 15.0) {
            Omega_S = sqrt_t_ratio;  // thick-S limit
        } else if (qS_dS < 1e-10) {
            Omega_S = 0.0;  // T → 0 limit
        } else {
            Omega_S = sqrt_t_ratio * std::cosh(qS_dS) / std::sinh(qS_dS);
        }
    }

    // Effective pair-breaking including interface barrier and S-layer impedance:
    // alpha = gamma * K / (1 + gamma_B * K + Omega_S(T))
    // No eta = xi_S/d_S prefactor — gamma absorbs coupling strength.
    std::complex<double> one(1.0, 0.0);
    std::complex<double> alpha = c->gamma * c->K / (one + c->gamma_B * c->K + Omega_S);

    // Scale by temperature: the Matsubara frequency normalization
    alpha *= c->Tc0 / (2.0 * M_PI * T);
    alpha += c->lambda_dep;

    std::complex<double> half(0.5, 0.0);
    auto psi_arg = supermag::digamma(half + alpha);
    auto psi_half = supermag::digamma(half);

    return log_ratio - (psi_arg - psi_half).real();
}

extern "C" {

int supermag_proximity_solve_tc(
    const supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double *tc_out)
{
    if (!params || !tc_out)
        return SUPERMAG_ERR_NULL_PTR;

    double d_F = params->d_F;
    double xi_F = params->xi_F;
    double Tc0 = params->Tc0;

    if (Tc0 <= 0.0 || xi_F <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Compute kernel  [EQ-2, EQ-3]
    // 0-junction: coth kernel (S/F bilayer, PHASE_ZERO)
    // π-junction: tanh kernel (semi-infinite F, PHASE_PI)
    std::complex<double> K;
    if (params->phase == SUPERMAG_PHASE_ZERO) {
        K = supermag::kernel_coth(d_F, xi_F);
    } else if (params->phase == SUPERMAG_PHASE_PI) {
        K = supermag::kernel_tanh(d_F, xi_F);
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    double lambda_dep = supermag_depairing_total(depairing);
    double T_min = 0.01;
    double T_max = Tc0;

    if (params->model == SUPERMAG_MODEL_THIN_S) {
        ThinSContext ctx;
        ctx.Tc0 = Tc0;
        ctx.gamma = params->gamma;
        ctx.lambda_dep = lambda_dep;
        ctx.K = K;

        double root = supermag::root_scalar_solve(thin_s_equation, &ctx, T_min, T_max, 1e-9);
        if (std::isnan(root)) {
            *tc_out = 0.0;
        } else {
            *tc_out = root;
        }
    } else if (params->model == SUPERMAG_MODEL_FOMINOV) {
        FominovContext ctx;
        ctx.Tc0 = Tc0;
        ctx.d_S = params->d_S;
        ctx.xi_S = params->xi_S;
        ctx.gamma = params->gamma;
        ctx.gamma_B = params->gamma_B;
        ctx.lambda_dep = lambda_dep;
        ctx.K = K;

        double root = supermag::root_scalar_solve(fominov_determinant, &ctx, T_min, T_max, 1e-9);
        if (std::isnan(root)) {
            *tc_out = 0.0;
        } else {
            *tc_out = root;
        }
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    return SUPERMAG_OK;
}

int supermag_proximity_solve_tc_batch(
    const supermag_proximity_params_t *params,
    const double *d_F_array,
    int n_dF,
    const supermag_depairing_t *depairing,
    double *tc_out)
{
    if (!params || !d_F_array || !tc_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_dF < 1)
        return SUPERMAG_ERR_INVALID_DIM;

    // Make a mutable copy of params to vary d_F per point
    supermag_proximity_params_t local_params;
    std::memcpy(&local_params, params, sizeof(local_params));

    for (int i = 0; i < n_dF; ++i) {
        local_params.d_F = d_F_array[i];
        int rc = supermag_proximity_solve_tc(&local_params, depairing, &tc_out[i]);
        if (rc != SUPERMAG_OK)
            return rc;
    }

    return SUPERMAG_OK;
}

} // extern "C"
