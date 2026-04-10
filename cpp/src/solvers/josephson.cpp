// STUB — Josephson current-phase relation solver
// TODO: Compute CPR for S/F/S junctions, detect 0-pi transitions.

#include "supermag/josephson.h"

extern "C" {

int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T,
    int n_phases, double* phase_arr, double* current_out)
{
    (void)d_F; (void)xi_F; (void)E_ex; (void)T;
    (void)n_phases; (void)phase_arr; (void)current_out;
    // TODO: Implement Josephson CPR calculation
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
