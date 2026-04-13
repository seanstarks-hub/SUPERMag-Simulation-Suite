#include "supermag/proximity.h"
#include "supermag/constants.h"
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern "C" {

double supermag_depairing_total(const supermag_depairing_t *dp) {
    if (!dp) return 0.0;
    return dp->ag + dp->zeeman + dp->orbital + dp->spin_orbit;
}

int supermag_depairing_compute(
    const supermag_depairing_input_t *input,
    supermag_depairing_t *output)
{
    if (!input || !output)
        return SUPERMAG_ERR_NULL_PTR;
    if (input->Tc0 <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // CODATA constants — all in SI
    const double hbar = supermag_const_hbar();    // J·s
    const double kB   = supermag_const_kB();      // J/K
    const double muB  = supermag_const_mu_B();    // J/T
    const double e    = supermag_const_e();       // C

    // Convert input units to SI for internal computation:
    //   Gamma_s, Gamma_so:  meV → J  (× 1.602176634e-22)
    //   D:  nm^2/ps → m^2/s  (× 1e-6)
    //   thickness:  nm → m  (× 1e-9)
    //   H: already Tesla
    //   Tc0: already Kelvin
    const double meV_to_J  = 1.602176634e-22;
    const double nm2ps_to_m2s = 1.0e-6;
    const double nm_to_m   = 1.0e-9;

    double Gamma_s_SI  = input->Gamma_s  * meV_to_J;   // J
    double Gamma_so_SI = input->Gamma_so * meV_to_J;   // J
    double D_SI        = input->D        * nm2ps_to_m2s; // m^2/s
    double thick_SI    = input->thickness * nm_to_m;    // m
    double kBTc0       = kB * input->Tc0;               // J

    // EQ-7A: Abrikosov-Gorkov  λ_AG = Γ_s / (2·kB·Tc0)
    output->ag = Gamma_s_SI / (2.0 * kBTc0);

    // EQ-7B: Zeeman  λ_Z = (μ_B·H)² / (2π·kB·Tc0)²
    double muBH = muB * input->H;
    double denom_z = 2.0 * M_PI * kBTc0;
    output->zeeman = (muBH * muBH) / (denom_z * denom_z);

    // EQ-7C: Orbital  λ_orb = D·(e·H)²·d² / (3·ℏ²·2π·kB·Tc0)
    double eH = e * input->H;
    output->orbital = D_SI * eH * eH * thick_SI * thick_SI
                    / (3.0 * hbar * hbar * 2.0 * M_PI * kBTc0);

    // EQ-7D: Spin-orbit  λ_SO = Γ_so / (2·kB·Tc0)
    output->spin_orbit = Gamma_so_SI / (2.0 * kBTc0);

    return SUPERMAG_OK;
}

}
