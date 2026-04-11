"""
BdG (Bogoliubov-de Gennes) tight-binding solver — Python interface.

Constructs and diagonalizes the BdG Hamiltonian for microscopic treatment
of S/F heterostructures. Needed when quasiclassical approximations break
down (ultra-thin layers, strong spin-orbit coupling, interface roughness).

The 2N×2N Nambu-doubled Hamiltonian is:

    H_BdG = [[ H_e↑,    Δ·I  ],
             [ Δ*·I,  -H_e↓* ]]

where H_eσ = -t·tridiag(1,0,1) + (-μ ± E_ex)·I.

See docs/theory/bdg_discretization.ipynb for details.

References:
    de Gennes, P.G. (1966). Superconductivity of Metals and Alloys.
"""

import numpy as np


_USE_NATIVE = False
try:
    from supermag._native import _bdg_solve as _native_bdg_solve
    _USE_NATIVE = True
except ImportError:
    pass


def solve(n_sites, t_hop, Delta, E_ex, mu=0.0):
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
    mu : float, optional
        Chemical potential (meV). Default: 0.0.

    Returns
    -------
    eigenvalues : numpy.ndarray
        BdG eigenvalues (meV), sorted. Length 2*n_sites.
    """
    if _USE_NATIVE and mu == 0.0:
        return _native_bdg_solve(n_sites, t_hop, Delta, E_ex)

    # Pure Python fallback
    N = int(n_sites)
    Delta_eV = Delta * 1e-3  # meV → eV
    E_ex_eV = E_ex * 1e-3    # meV → eV
    mu_eV = mu * 1e-3         # meV → eV

    H = np.zeros((2 * N, 2 * N), dtype=complex)

    # On-site terms
    for i in range(N):
        H[i, i] = -mu_eV + E_ex_eV            # electron (spin-up): -µ + h
        H[N + i, N + i] = mu_eV + E_ex_eV      # hole (spin-down): +µ + h
        # Pairing
        H[i, N + i] = Delta_eV
        H[N + i, i] = np.conj(Delta_eV)

    # Hopping (nearest-neighbor)
    for i in range(N - 1):
        H[i, i + 1] = -t_hop            # electron hopping
        H[i + 1, i] = -t_hop
        H[N + i, N + i + 1] = t_hop     # hole hopping (sign flip)
        H[N + i + 1, N + i] = t_hop

    eigenvalues = np.linalg.eigvalsh(H)
    return eigenvalues * 1e3  # eV → meV
