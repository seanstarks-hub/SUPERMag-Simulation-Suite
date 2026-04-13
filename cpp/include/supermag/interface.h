#ifndef SUPERMAG_INTERFACE_H
#define SUPERMAG_INTERFACE_H

#include "error.h"
#include "proximity.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Kupriyanov-Lukichev boundary condition coefficient.
 * Returns gamma_KL = (rho_S * xi_S) / (rho_F * xi_F) * (1 / R_B)
 * for given material and interface parameters.
 *
 * @param rho_S     S-layer resistivity (Ohm·nm)
 * @param xi_S      S-layer coherence length (nm)
 * @param rho_F     F-layer resistivity (Ohm·nm)
 * @param xi_F      F-layer coherence length (nm)
 * @param R_B       Interface barrier resistance (Ohm·nm^2), must be > 0
 * @param gamma_out Output: KL boundary parameter
 */
int supermag_kl_boundary(
    double rho_S, double xi_S,
    double rho_F, double xi_F,
    double R_B,
    double *gamma_out
);

/**
 * Apply roughness correction to the effective kernel.
 * Roughness suppresses the pair amplitude at the interface via
 * a Debye-Waller-like factor:
 *   K_rough = K * exp(-(delta/xi_F)^2 / 2)
 *
 * @param K_real        Real part of kernel (modified in-place)
 * @param K_imag        Imaginary part of kernel (modified in-place)
 * @param roughness     Roughness parameters
 * @param xi_F          Ferromagnet coherence length (nm)
 */
int supermag_apply_roughness(
    double *K_real, double *K_imag,
    const supermag_roughness_t *roughness,
    double xi_F
);

#ifdef __cplusplus
}
#endif

#endif
