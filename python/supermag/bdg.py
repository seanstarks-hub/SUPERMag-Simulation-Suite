"""
STUB — BdG (Bogoliubov-de Gennes) tight-binding solver Python interface.

Constructs and diagonalizes the BdG Hamiltonian for microscopic treatment
of S/F heterostructures. Needed when quasiclassical approximations break
down (ultra-thin layers, strong spin-orbit coupling, interface roughness).

See docs/theory/bdg_discretization.md for details.
"""

import numpy as np


def solve(n_sites, t_hop, Delta, E_ex):
    """
    Diagonalize the BdG Hamiltonian on a 1D tight-binding lattice.

    Parameters
    ----------
    n_sites : int
        Number of lattice sites.
    t_hop : float
        Nearest-neighbor hopping energy (eV).
    Delta : float
        Superconducting pairing potential (meV).
    E_ex : float
        Exchange splitting in ferromagnet region (meV).

    Returns
    -------
    eigenvalues : numpy.ndarray
        BdG eigenvalues (meV).

    Raises
    ------
    NotImplementedError
        This solver is not yet implemented.
    """
    raise NotImplementedError(
        "BdG solver not yet implemented. See docs/theory/bdg_discretization.md"
    )
