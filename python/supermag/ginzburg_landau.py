"""
STUB — Ginzburg-Landau free energy solver Python interface.

Minimizes the GL functional for vortex states, mixed-state configurations,
and domain structures near Tc.

See docs/theory/ginzburg_landau.md for details.
"""

import numpy as np


def minimize(alpha, beta, kappa, nx, ny, dx):
    """
    Minimize the Ginzburg-Landau free energy on a 2D grid.

    Parameters
    ----------
    alpha : float
        GL alpha parameter (temperature dependent).
    beta : float
        GL beta parameter.
    kappa : float
        GL kappa = lambda_L / xi.
    nx, ny : int
        Grid dimensions.
    dx : float
        Grid spacing (nm).

    Returns
    -------
    psi : numpy.ndarray
        Complex order parameter on grid, shape (ny, nx).

    Raises
    ------
    NotImplementedError
        This solver is not yet implemented.
    """
    raise NotImplementedError(
        "GL solver not yet implemented. See docs/theory/ginzburg_landau.md"
    )
