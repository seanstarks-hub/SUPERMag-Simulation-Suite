// Individual depairing channel computations.
// Each function computes one dimensionless pair-breaking parameter
// from laboratory-measurable physical inputs.
// All functions use T_kelvin as the temperature denominator.

#include "supermag/depairing.h"
#include "supermag/proximity.h"
#include "supermag/constants.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern "C" {

// AG pair-breaking from spin-flip scattering
//   lambda_AG = gamma_s / (2 * kB * T)
double supermag_depairing_ag(double gamma_s_meV, double T_kelvin) {
    if (T_kelvin <= 0.0) return 0.0;
    double kB = supermag_const_kB();  // meV/K
    return gamma_s_meV / (2.0 * kB * T_kelvin);
}

// Zeeman pair-breaking
//   lambda_Z = (mu_B * H)^2 / (2*pi*kB*T)^2
double supermag_depairing_zeeman(double H_tesla, double T_kelvin) {
    if (T_kelvin <= 0.0) return 0.0;
    double kB = supermag_const_kB();
    double mu_B = supermag_const_mu_B();
    double denom = 2.0 * M_PI * kB * T_kelvin;
    return (mu_B * H_tesla) * (mu_B * H_tesla) / (denom * denom);
}

// Orbital pair-breaking (perpendicular field)
//   lambda_orb_perp = D*(e*H)^2*d^2 / (3*hbar^2*(2*pi*kB*T))
double supermag_depairing_orbital_perp(double D_nm2ps, double H_tesla,
                                       double thickness_nm, double T_kelvin) {
    if (T_kelvin <= 0.0) return 0.0;
    double kB = supermag_const_kB();
    double hbar = supermag_const_hbar();
    double e = supermag_const_e();
    double denom = 3.0 * hbar * hbar * 2.0 * M_PI * kB * T_kelvin;
    if (denom == 0.0) return 0.0;
    return D_nm2ps * (e * H_tesla) * (e * H_tesla) * thickness_nm * thickness_nm / denom;
}

// Orbital pair-breaking (parallel field)
//   lambda_orb_par = D*(e*H)^2*d^2 / (12*hbar^2*(2*pi*kB*T))
double supermag_depairing_orbital_par(double D_nm2ps, double H_tesla,
                                      double thickness_nm, double T_kelvin) {
    if (T_kelvin <= 0.0) return 0.0;
    double kB = supermag_const_kB();
    double hbar = supermag_const_hbar();
    double e = supermag_const_e();
    double denom = 12.0 * hbar * hbar * 2.0 * M_PI * kB * T_kelvin;
    if (denom == 0.0) return 0.0;
    return D_nm2ps * (e * H_tesla) * (e * H_tesla) * thickness_nm * thickness_nm / denom;
}

// Spin-orbit coupling pair-breaking
//   lambda_SOC = Gamma_so / (2 * kB * T)
double supermag_depairing_soc(double Gamma_so_meV, double T_kelvin) {
    if (T_kelvin <= 0.0) return 0.0;
    double kB = supermag_const_kB();
    return Gamma_so_meV / (2.0 * kB * T_kelvin);
}

// Compute all depairing channels from physical inputs
int supermag_depairing_from_physical(
    double gamma_s_meV, double H_tesla,
    double D_nm2ps, double thickness_nm,
    double Gamma_so_meV, double T_kelvin,
    supermag_depairing_t *output)
{
    if (!output)
        return SUPERMAG_ERR_NULL_PTR;
    if (T_kelvin <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    output->ag = supermag_depairing_ag(gamma_s_meV, T_kelvin);
    output->zeeman = supermag_depairing_zeeman(H_tesla, T_kelvin);
    output->orbital = supermag_depairing_orbital_perp(D_nm2ps, H_tesla, thickness_nm, T_kelvin);
    output->spin_orbit = supermag_depairing_soc(Gamma_so_meV, T_kelvin);

    return SUPERMAG_OK;
}

} // extern "C"
