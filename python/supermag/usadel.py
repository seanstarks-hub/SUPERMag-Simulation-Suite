"""
STUB — Usadel diffusive-limit solver Python interface.

Solves the Usadel equation for quasiclassical Green's functions in the
diffusive limit with Kupriyanov-Lukichev boundary conditions at S/F interfaces.

See docs/theory/usadel_equations.md for full mathematical treatment.
"""

import numpy as np


def solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex, n_grid=200):
    """
    Solve the Usadel equation for an S/F bilayer.

    Computes the self-consistent superconducting order parameter Delta(x)
    across the S/F heterostructure.

    Parameters
    ----------
    Tc0 : float
        Bulk superconductor critical temperature (K).
    d_S : float
        Superconductor layer thickness (nm).
    d_F : float
        Ferromagnet layer thickness (nm).
    xi_S : float
        Superconductor coherence length (nm).
    xi_F : float
        Ferromagnet coherence length (nm).
    E_ex : float
        Exchange energy in ferromagnet (meV).
    n_grid : int, optional
        Number of spatial grid points. Default: 200.

    Returns
    -------
    x : numpy.ndarray
        Position array (nm).
    Delta : numpy.ndarray
        Order parameter profile (meV).

    Raises
    ------
    NotImplementedError
        This solver is not yet implemented.
    """
    raise NotImplementedError(
        "Usadel solver not yet implemented. See docs/theory/usadel_equations.md"
    )
