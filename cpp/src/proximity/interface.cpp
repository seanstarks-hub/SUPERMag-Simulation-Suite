// Interface boundary condition and roughness corrections.
//
// Kupriyanov-Lukichev boundary:
//   gamma_KL = (rho_S * xi_S) / (rho_F * xi_F) * (1 / R_B)
//
// Roughness (Debye-Waller-like suppression):
//   K_rough = K * exp(-(delta/xi_F)^2 / 2)

#include "supermag/interface.h"
#include <cmath>

extern "C" {

int supermag_kl_boundary(
    double rho_S, double xi_S,
    double rho_F, double xi_F,
    double R_B,
    double *gamma_out)
{
    if (!gamma_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (rho_F <= 0.0 || xi_F <= 0.0 || R_B <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    *gamma_out = (rho_S * xi_S) / (rho_F * xi_F * R_B);
    return SUPERMAG_OK;
}

int supermag_apply_roughness(
    double *K_real, double *K_imag,
    const supermag_roughness_t *roughness,
    double xi_F)
{
    if (!K_real || !K_imag)
        return SUPERMAG_ERR_NULL_PTR;
    if (!roughness)
        return SUPERMAG_OK;  // no roughness to apply
    if (xi_F <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    double ratio = roughness->delta_roughness / xi_F;
    double factor = std::exp(-ratio * ratio / 2.0);

    *K_real *= factor;
    *K_imag *= factor;
    return SUPERMAG_OK;
}

} // extern "C"
