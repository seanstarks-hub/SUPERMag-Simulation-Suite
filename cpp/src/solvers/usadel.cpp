// STUB — Usadel diffusive-limit solver
// TODO: Implement self-consistency loop for Delta(x) with Kupriyanov-Lukichev
//       boundary conditions at S/F interfaces. Uses tridiagonal solver for
//       1D finite-difference discretization.

#include "supermag/usadel.h"

extern "C" {

int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    int n_grid, double* Delta_out, double* x_out)
{
    (void)Tc0; (void)d_S; (void)d_F;
    (void)xi_S; (void)xi_F; (void)E_ex;
    (void)n_grid; (void)Delta_out; (void)x_out;
    // TODO: Implement Usadel equation solver
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
