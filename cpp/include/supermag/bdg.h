#ifndef SUPERMAG_BDG_H
#define SUPERMAG_BDG_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* BdG tight-binding Hamiltonian diagonalization
 * For ultra-thin layers, strong spin-orbit coupling, interface roughness.
 *
 * mu: chemical potential (meV).  0.0 reproduces legacy behaviour.
 * eigenvectors_out: if non-NULL, filled with eigenvector matrix (dim×dim, row-major).
 *                   Must be preallocated to (2*n_sites)^2 doubles.
 *                   Pass NULL to skip eigenvector computation. */
int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues,
    double* eigenvectors_out
);

#ifdef __cplusplus
}
#endif

#endif
