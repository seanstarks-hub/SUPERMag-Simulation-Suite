// STUB — Spin-triplet superconductivity solver
// TODO: Long-range triplet correlations from inhomogeneous magnetization.

#include "supermag/triplet.h"

extern "C" {

int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    int n_grid, double* f_triplet_out, double* x_out)
{
    (void)n_layers; (void)thicknesses; (void)magnetization_angles;
    (void)n_grid; (void)f_triplet_out; (void)x_out;
    // TODO: Implement spin-triplet solver
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
