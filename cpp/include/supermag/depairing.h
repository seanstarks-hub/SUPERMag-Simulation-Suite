#ifndef SUPERMAG_DEPAIRING_H
#define SUPERMAG_DEPAIRING_H

#include "error.h"
#include "proximity.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Individual depairing channel computations.
 * Each returns the dimensionless pair-breaking parameter for one channel.
 * All use T_kelvin as the temperature denominator (not Tc0).
 */

/* Abrikosov-Gor'kov (AG) pair-breaking from spin-flip scattering.
 * lambda_AG = gamma_s_meV / (2 * kB * T_kelvin) */
double supermag_depairing_ag(double gamma_s_meV, double T_kelvin);

/* Zeeman pair-breaking.
 * lambda_Z = (mu_B * H_tesla)^2 / (2*pi*kB*T_kelvin)^2 */
double supermag_depairing_zeeman(double H_tesla, double T_kelvin);

/* Orbital pair-breaking (perpendicular field).
 * lambda_orb_perp = D*(e*H)^2*d^2 / (3*hbar^2*(2*pi*kB*T)) */
double supermag_depairing_orbital_perp(double D_nm2ps, double H_tesla,
                                       double thickness_nm, double T_kelvin);

/* Orbital pair-breaking (parallel field).
 * lambda_orb_par = D*(e*H)^2*d^2 / (12*hbar^2*(2*pi*kB*T)) */
double supermag_depairing_orbital_par(double D_nm2ps, double H_tesla,
                                      double thickness_nm, double T_kelvin);

/* Spin-orbit coupling pair-breaking.
 * lambda_SOC = Gamma_so / (2 * kB * T_kelvin) */
double supermag_depairing_soc(double Gamma_so_meV, double T_kelvin);

/**
 * Compute all depairing channels from physical laboratory inputs.
 * Fills the output depairing_t struct with all four channels.
 *
 * @param gamma_s_meV    Spin-flip scattering rate (meV)
 * @param H_tesla        Applied magnetic field (T)
 * @param D_nm2ps        Diffusion coefficient (nm^2/ps)
 * @param thickness_nm   Film thickness (nm)
 * @param Gamma_so_meV   Spin-orbit scattering rate (meV)
 * @param T_kelvin       Temperature (K), must be > 0
 * @param output         Output: filled depairing struct
 */
int supermag_depairing_from_physical(
    double gamma_s_meV, double H_tesla,
    double D_nm2ps, double thickness_nm,
    double Gamma_so_meV, double T_kelvin,
    supermag_depairing_t *output
);

#ifdef __cplusplus
}
#endif

#endif
