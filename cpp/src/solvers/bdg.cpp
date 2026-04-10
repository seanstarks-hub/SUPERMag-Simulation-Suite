// STUB — BdG tight-binding Hamiltonian solver
// TODO: Construct and diagonalize Bogoliubov-de Gennes Hamiltonian.

#include "supermag/bdg.h"

extern "C" {

int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex,
    double* eigenvalues_out, int* n_eigenvalues)
{
    (void)n_sites; (void)t_hop; (void)Delta; (void)E_ex;
    (void)eigenvalues_out; (void)n_eigenvalues;
    // TODO: Implement BdG solver
    return SUPERMAG_ERR_NO_CONVERGE;
}

}
