#ifndef SUPERMAG_JOSEPHSON_H
#define SUPERMAG_JOSEPHSON_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Josephson current-phase relation (CPR) for S/F/S junctions.
 * 0-pi transition detection as function of F-layer thickness and temperature. */
int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T,
    int n_phases, double* phase_arr, double* current_out
);

#ifdef __cplusplus
}
#endif

#endif
