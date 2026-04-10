"""
STUB — Eilenberger clean-limit solver Python interface.

Solves the Eilenberger equation using Riccati parameterization
for quasiclassical Green's functions in the clean limit.

See docs/theory/usadel_equations.ipynb for the diffusive counterpart.
"""

import numpy as np


def solve(Tc0, d_S, d_F, xi_S, E_ex, n_grid=200):
    """
    Solve the Eilenberger equation for an S/F bilayer.

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
    E_ex : float
        Exchange energy in ferromagnet (meV).
    n_grid : int, optional
        Number of spatial grid points. Default: 200.

    Returns
    -------
    x : numpy.ndarray
        Position array (nm).
    f : numpy.ndarray
        Anomalous Green's function profile.

    Raises
    ------
    NotImplementedError
        This solver is not yet implemented.
    """
    raise NotImplementedError(
        "Eilenberger solver not yet implemented. See docs/theory/usadel_equations.ipynb"
    )
