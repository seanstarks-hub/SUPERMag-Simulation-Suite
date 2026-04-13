// Proximity effect Tc solver — digamma-based self-consistency equations.
//
// Thin-S model:
//   ln(Tc0/T) = Re[ psi(1/2 + alpha) - psi(1/2) ]
//   alpha = gamma / (gamma_B + K) * Tc0/(2*pi*T) + lambda_dep
//   K in the denominator: as d_F->0, K->inf, alpha->0, Tc->Tc0.
//
// Fominov model (PRB 66, 014507):
//   Same form with alpha = gamma / (gamma_B + K + Omega_S(T)) * Tc0/(2*pi*T)
//   Omega_S = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S/xi_S)
//   Reduces to thin-S when Omega_S = 0.

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
// Forward declaration from spin_active.cpp
std::complex<double> apply_spin_active(const supermag_spin_active_t *sa, std::complex<double> K);
}

// Context struct passed to the root-finding callback
struct ThinSContext {
    double Tc0;
    double gamma;
    double gamma_B;
    double lambda_dep;
    std::complex<double> K;  // kernel value at this d_F
};

// F(T) = ln(Tc0/T) - Re[ psi(1/2 + A(T)) - psi(1/2) ]
// where A(T) = gamma / (gamma_B + K) * Tc0 / (2*pi*T) + lambda_dep
//
// The pair-breaking parameter alpha = gamma / (gamma_B + K) puts K in the
// denominator: as d_F -> 0, the coth kernel K -> infinity, so alpha -> 0
// and Tc -> Tc0 (physical thin-F limit). gamma is the coupling strength
// (numerator), gamma_B is the interface barrier (denominator).
//
// This is the thin-S limit of the Fominov model (Omega_S = 0).

static double thin_s_equation(double T, void *ctx) {
    auto *c = static_cast<ThinSContext*>(ctx);
    if (T <= 0.0) return -1.0;

    double log_ratio = std::log(c->Tc0 / T);

    // Complex pair-breaking parameter: A = gamma / (gamma_B + K) * Tc0 / (2*pi*T)  [EQ-4]
    std::complex<double> half(0.5, 0.0);
    std::complex<double> alpha = c->gamma / (c->gamma_B + c->K);
    std::complex<double> A_complex = alpha * c->Tc0 / (2.0 * M_PI * T);
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
//   alpha = gamma / (gamma_B + K + Omega_S(T))  * Tc0 / (2*pi*T)
//   Omega_S(T) = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S/xi_S)
//
// K is in the denominator: as d_F -> 0, K (coth kernel) diverges,
// alpha -> 0, and Tc -> Tc0 (physical thin-F limit).
// This reduces to the thin-S equation (EQ-4) when Omega_S = 0.
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

    // Effective pair-breaking with K in denominator:
    // alpha = gamma / (gamma_B + K + Omega_S(T))
    // As d_F -> 0, K (coth kernel) diverges, so alpha -> 0 and Tc -> Tc0.
    std::complex<double> alpha = c->gamma / (c->gamma_B + c->K + Omega_S);

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
    // Standard bilayer: 0-junction coth, π-junction tanh
    // Extended geometries use Phase 1 kernel functions (EQ-13, EQ-14, EQ-15)
    std::complex<double> K;
    if (params->geometry == SUPERMAG_GEOM_TRILAYER) {
        // S/N/F trilayer kernel [EQ-13]
        if (!params->geom_params)
            return SUPERMAG_ERR_NULL_PTR;
        auto *tri = static_cast<const supermag_trilayer_params_t*>(params->geom_params);
        double Kr, Ki;
        int krc = supermag_proximity_kernel_snf(d_F, xi_F, tri, params->phase, &Kr, &Ki);
        if (krc != SUPERMAG_OK) return krc;
        K = std::complex<double>(Kr, Ki);
    } else if (params->geometry == SUPERMAG_GEOM_GRADED) {
        // Graded ferromagnet kernel [EQ-14]
        if (!params->geom_params)
            return SUPERMAG_ERR_NULL_PTR;
        auto *grade = static_cast<const supermag_graded_params_t*>(params->geom_params);
        double Kr, Ki;
        int krc = supermag_proximity_kernel_graded(d_F, xi_F, grade, params->phase, &Kr, &Ki);
        if (krc != SUPERMAG_OK) return krc;
        K = std::complex<double>(Kr, Ki);
    } else if (params->geometry == SUPERMAG_GEOM_DOMAINS) {
        // Magnetic domain kernel [EQ-15]
        if (!params->geom_params)
            return SUPERMAG_ERR_NULL_PTR;
        auto *dom = static_cast<const supermag_domain_params_t*>(params->geom_params);
        double Kr, Ki;
        int krc = supermag_proximity_kernel_domains(d_F, xi_F, params->E_ex, dom,
                                                     params->phase, &Kr, &Ki);
        if (krc != SUPERMAG_OK) return krc;
        K = std::complex<double>(Kr, Ki);
    } else if (params->phase == SUPERMAG_PHASE_ZERO) {
        K = supermag::kernel_coth(d_F, xi_F);
    } else if (params->phase == SUPERMAG_PHASE_PI) {
        K = supermag::kernel_tanh(d_F, xi_F);
    } else {
        return SUPERMAG_ERR_INVALID_MODEL;
    }

    // Apply spin-active interface correction to kernel
    if (params->spin_active) {
        K = supermag::apply_spin_active(params->spin_active, K);
    }

    double lambda_dep = supermag_depairing_total(depairing);
    double T_min = 0.01;
    double T_max = Tc0;

    // Self-consistency equation selection:
    // Non-bilayer geometries use thin_s equation (kernel is the only difference)
    int eq_model = params->model;
    // Fominov multimode falls back to single-mode Fominov for now
    if (eq_model == SUPERMAG_MODEL_FOMINOV_MULTI)
        eq_model = SUPERMAG_MODEL_FOMINOV;

    if (eq_model == SUPERMAG_MODEL_THIN_S) {
        ThinSContext ctx;
        ctx.Tc0 = Tc0;
        ctx.gamma = params->gamma;
        ctx.gamma_B = params->gamma_B;
        ctx.lambda_dep = lambda_dep;
        ctx.K = K;

        double root = supermag::root_scalar_solve(thin_s_equation, &ctx, T_min, T_max, 1e-9);
        if (std::isnan(root)) {
            *tc_out = 0.0;
        } else {
            *tc_out = root;
        }
    } else if (eq_model == SUPERMAG_MODEL_FOMINOV) {
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
