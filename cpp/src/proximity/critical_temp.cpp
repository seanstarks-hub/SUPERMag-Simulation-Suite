// Proximity effect Tc solver — digamma-based self-consistency equations.
//
// Thin-S model:
//   ln(Tc0/T) = Re[ psi(1/2 + alpha) - psi(1/2) ]
//   alpha = gamma / (gamma_B + K) * Tc0 / (2*pi*T) + Lambda_dep
//   K is the kernel (coth for 0-junction), gamma is coupling, gamma_B is barrier.
//   K in the denominator: as d_F→0, K (coth)→∞, so alpha→0 and Tc→Tc0.
//   Find highest root F(T) = 0 via Brent's method.
//
// Fominov model (PRB 66, 014507):
//   Same digamma equation with T-dependent S-layer impedance:
//   alpha = gamma / (gamma_B + K + Omega_S(T)) * Tc0 / (2*pi*T) + Lambda_dep
//   Omega_S(T) = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S / xi_S)
//   Reduces to thin-S when d_S→∞ (Omega_S→0) and gamma_B→0.

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
//
// Self-consistency equation for S/F bilayer Tc:
//   ln(Tc0/Tc) = Re[ psi(1/2 + alpha) - psi(1/2) ]
//   alpha = gamma / (gamma_B + K) * Tc0 / (2*pi*T) + lambda_dep
//
// K = q*coth(q*d_F) for 0-junction (diverges as d_F→0 → alpha→0 → Tc→Tc0).
// gamma is the coupling strength (numerator), gamma_B is the interface barrier
// (additive with K in denominator).
//
// Reduces to Fominov (EQ-5) by adding Omega_S(T) to the denominator.

static double thin_s_equation(double T, void *ctx) {
    auto *c = static_cast<ThinSContext*>(ctx);
    if (T <= 0.0) return -1.0;

    double log_ratio = std::log(c->Tc0 / T);

    // Pair-breaking parameter: A = gamma / (gamma_B + K) * Tc0 / (2*pi*T)  [EQ-4]
    // K is in the denominator: larger K (thicker F) → smaller A → less suppression.
    std::complex<double> half(0.5, 0.0);
    std::complex<double> A_complex = c->gamma / (c->gamma_B + c->K) * c->Tc0 / (2.0 * M_PI * T);
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
//
// The Fominov (PRB 66 014507) Tc equation for an S/F bilayer:
//   ln(Tc0/Tc) = Re[ psi(1/2 + alpha) - psi(1/2) ]
// where
//   alpha = gamma / (gamma_B + K + Omega_S(T)) * Tc0 / (2*pi*T) + lambda_dep
//   Omega_S(T) = sqrt(T/Tc0) * coth(sqrt(T/Tc0) * d_S/xi_S)
//
// K is in the denominator alongside Omega_S and gamma_B.
// As d_F->0, K (coth) diverges, so alpha->0 and Tc->Tc0.
// As d_S->inf, Omega_S->0, recovering the thin-S formula.
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

    // Effective pair-breaking with K in the denominator:
    // alpha = gamma / (gamma_B + K + Omega_S(T))  [EQ-5]
    // K (coth) diverges as d_F->0, driving alpha->0 and Tc->Tc0.
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
    // K appears in the denominator of α, so the coth divergence at d_F→0
    // drives α→0 and Tc→Tc0 (correct physics).
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
