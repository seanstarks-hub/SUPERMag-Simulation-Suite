"""
STUB — Josephson junction CPR solver Python interface.

Computes current-phase relations for S/F/S, S/F/I/F/S junctions.
Detects 0-pi transitions as function of F-layer thickness and temperature.

See docs/theory/josephson_cpr.md for details.
"""

import numpy as np


def current_phase_relation(d_F, xi_F, E_ex, T, n_phases=100):
    """
    Compute Josephson current-phase relation I(phi).

    Parameters
    ----------
    d_F : float
        Ferromagnet thickness (nm).
    xi_F : float
        Ferromagnet coherence length (nm).
    E_ex : float
        Exchange energy (meV).
    T : float
        Temperature (K).
    n_phases : int, optional
        Number of phase points from 0 to 2*pi.

    Returns
    -------
    phi : numpy.ndarray
        Phase difference array (rad).
    I : numpy.ndarray
        Supercurrent array (normalized).

    Raises
    ------
    NotImplementedError
        This solver is not yet implemented.
    """
    raise NotImplementedError(
        "Josephson CPR solver not yet implemented. See docs/theory/josephson_cpr.md"
    )
