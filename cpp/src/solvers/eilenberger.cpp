// STUB — Eilenberger clean-limit solver
// TODO: Implement Riccati parameterization for quasiclassical Green's function.

#include "supermag/eilenberger.h"

extern "C" {

int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    int n_grid, double* f_out, double* x_out)
{
    (void)Tc0; (void)d_S; (void)d_F;
    (void)xi_S; (void)E_ex;
    (void)n_grid; (void)f_out; (void)x_out;
    // TODO: Implement Eilenberger equation solver
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
