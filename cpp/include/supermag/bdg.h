#ifndef SUPERMAG_BDG_H
#define SUPERMAG_BDG_H

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: BdG tight-binding Hamiltonian diagonalization
 * For ultra-thin layers, strong spin-orbit coupling, interface roughness. */
int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex,
    double* eigenvalues_out, int* n_eigenvalues
);

#ifdef __cplusplus
}
#endif

#endif
