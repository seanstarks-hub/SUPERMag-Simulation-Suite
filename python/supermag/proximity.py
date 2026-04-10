"""
S/F Proximity Effect solver — Python interface.

Computes:
- Pair amplitude F(x) in the ferromagnet layer (decaying oscillation)
- Critical temperature Tc(d_F) for S/F bilayer (oscillatory suppression)

Physics:
    In a superconductor/ferromagnet bilayer, Cooper pairs leak into the
    ferromagnet. The exchange field causes the pair amplitude to oscillate
    with a characteristic wavelength ~xi_F (FFLO-like oscillation).
    This leads to non-monotonic Tc(d_F) behavior.

References:
    Buzdin, A.I. (1982). JETP Lett. 35, 178.
    Radović, Z. et al. (1991). Phys. Rev. B 44, 759.
"""

import numpy as np


# Try to import the native C++ extension; fall back to pure Python
_USE_NATIVE = False
try:
    from supermag._native import (
        _pair_amplitude as _native_pair_amplitude,
        _critical_temp as _native_critical_temp,
    )
    _USE_NATIVE = True
except ImportError:
    pass


def pair_amplitude(d_F, xi_F, F0=1.0, n_points=500):
    """
    Compute the pair amplitude F(x) in the ferromagnet layer of an S/F bilayer.

    The pair amplitude exhibits FFLO-like decaying oscillation:
        F(x) = F0 * exp(-x/xi_F) * cos(x/xi_F)

    Parameters
    ----------
    d_F : float
        Ferromagnet layer thickness (nm).
    xi_F : float
        Ferromagnet coherence length (nm). Typically 0.5–5 nm.
    F0 : float, optional
        Pair amplitude at the S/F interface (dimensionless). Default: 1.0.
    n_points : int, optional
        Number of spatial grid points. Default: 500.

    Returns
    -------
    x : numpy.ndarray
        Position array from 0 to d_F (nm), shape (n_points,).
    F : numpy.ndarray
        Pair amplitude F(x) (dimensionless), shape (n_points,).

    Examples
    --------
    >>> import supermag
    >>> x, F = supermag.pair_amplitude(d_F=10.0, xi_F=2.0)
    >>> print(f"F at interface: {F[0]:.3f}")
    F at interface: 1.000
    """
    if _USE_NATIVE:
        return _native_pair_amplitude(F0, xi_F, d_F, n_points)

    # Pure Python fallback
    x = np.linspace(0, d_F, n_points)
    F = F0 * np.exp(-x / xi_F) * np.cos(x / xi_F)
    return x, F


def critical_temperature(Tc0, d_S, d_F_array, E_ex, xi_S, xi_F):
    """
    Compute critical temperature Tc as a function of ferromagnet thickness d_F.

    Uses single-mode Buzdin approximation. Tc(d_F) shows oscillatory
    non-monotonic behavior due to FFLO-like pair amplitude oscillation
    in the ferromagnet layer.

    Parameters
    ----------
    Tc0 : float
        Bulk superconductor critical temperature (K). E.g., 9.2 K for Nb.
    d_S : float
        Superconductor layer thickness (nm).
    d_F_array : array_like
        Array of ferromagnet thicknesses (nm).
    E_ex : float
        Exchange energy in ferromagnet (meV). E.g., 256 meV for Fe.
    xi_S : float
        Superconductor coherence length (nm). E.g., 38 nm for Nb.
    xi_F : float
        Ferromagnet coherence length (nm). E.g., 0.7 nm for Fe.

    Returns
    -------
    Tc : numpy.ndarray
        Critical temperature for each d_F value (K), shape matching d_F_array.

    Examples
    --------
    >>> import numpy as np
    >>> import supermag
    >>> d_F = np.linspace(0.5, 20.0, 50)
    >>> Tc = supermag.critical_temperature(Tc0=9.2, d_S=50.0,
    ...     d_F_array=d_F, E_ex=256.0, xi_S=38.0, xi_F=0.7)
    """
    d_F_array = np.asarray(d_F_array, dtype=np.float64)

    if _USE_NATIVE:
        return _native_critical_temp(Tc0, d_S, xi_S, xi_F, E_ex, d_F_array)

    # Pure Python fallback: single-mode Buzdin approximation
    inv_xi_F = 1.0 / xi_F
    Tc_out = np.empty_like(d_F_array)

    for i, d_F in enumerate(d_F_array):
        if d_F < 0:
            Tc_out[i] = Tc0
            continue

        arg = d_F * inv_xi_F
        exp_decay = np.exp(-2.0 * arg)
        real_exp = exp_decay * np.cos(2.0 * arg)
        imag_exp = -exp_decay * np.sin(2.0 * arg)

        real_num = 1.0 - real_exp
        imag_num = -imag_exp

        kappa_real = xi_S / xi_F
        kappa_imag = xi_S / xi_F
        denom = kappa_real**2 + kappa_imag**2
        ratio_real = (real_num * kappa_real + imag_num * kappa_imag) / denom

        eta = (xi_S / d_S) * ratio_real
        Tc = Tc0 * (1.0 - eta)
        Tc_out[i] = np.clip(Tc, 0.0, Tc0)

    return Tc_out
