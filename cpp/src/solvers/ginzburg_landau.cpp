// STUB — Ginzburg-Landau free energy minimization
// TODO: Implement iterative minimization of GL functional on 2D grid.

#include "supermag/ginzburg_landau.h"

extern "C" {

int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    double* psi_real, double* psi_imag)
{
    (void)alpha; (void)beta; (void)kappa;
    (void)nx; (void)ny; (void)dx;
    (void)psi_real; (void)psi_imag;
    // TODO: Implement Ginzburg-Landau solver
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
