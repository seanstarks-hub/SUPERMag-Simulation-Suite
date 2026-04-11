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

Models:
    thin_s   — Thin-S limit with digamma self-consistency equation
    fominov  — PRB 66, 014507 (2002): includes gamma_B interface barrier

References:
    Buzdin, A.I. (1982). JETP Lett. 35, 178.
    Radović, Z. et al. (1991). Phys. Rev. B 44, 759.
    Fominov, Ya.V. et al. (2002). Phys. Rev. B 66, 014507.
"""

import numpy as np


# Try to import the native C++ extension; fall back to pure Python
_USE_NATIVE = False
try:
    from supermag._native import (
        _pair_amplitude as _native_pair_amplitude,
        _solve_tc_batch as _native_solve_tc_batch,
    )
    _USE_NATIVE = True
except ImportError:
    pass


def pair_amplitude(d_F, xi_F, phase="zero", n_points=500):
    """
    Compute the pair amplitude F(x) in the ferromagnet layer of an S/F bilayer.

    For phase="zero" (0-junction):
        F(x) = exp(-x/xi_F) * cos(x/xi_F)
    For phase="pi" (pi-junction):
        F(x) = exp(-x/xi_F) * sin(x/xi_F)

    Parameters
    ----------
    d_F : float
        Ferromagnet layer thickness (nm).
    xi_F : float
        Ferromagnet coherence length (nm). Typically 0.5–10 nm.
    phase : str, optional
        "zero" for 0-junction (coth kernel) or "pi" for pi-junction (tanh kernel).
        Default: "zero".
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
    if d_F <= 0:
        raise ValueError(f"d_F must be positive, got {d_F}")
    if xi_F <= 0:
        raise ValueError(f"xi_F must be positive, got {xi_F}")
    if phase not in ("zero", "pi"):
        raise ValueError(f"phase must be 'zero' or 'pi', got {phase!r}")
    if n_points < 2:
        raise ValueError(f"n_points must be >= 2, got {n_points}")

    if _USE_NATIVE:
        phase_int = 1 if phase == "pi" else 0
        return _native_pair_amplitude(d_F, xi_F, phase_int, n_points)

    # Pure Python fallback
    x = np.linspace(0, d_F, n_points)
    envelope = np.exp(-x / xi_F)
    if phase == "pi":
        F = envelope * np.sin(x / xi_F)
    else:
        F = envelope * np.cos(x / xi_F)
    return x, F


def critical_temperature(Tc0, d_S, d_F_array, E_ex, xi_S, xi_F,
                         gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
                         model="thin_s", phase="zero",
                         depairing=None):
    """
    Compute critical temperature Tc as a function of ferromagnet thickness d_F.

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
    gamma : float, optional
        Interface transparency parameter (dimensionless). Default: 0.3.
    gamma_B : float, optional
        Interface barrier parameter (Fominov model). Default: 0.0.
    D_F : float, optional
        Diffusion coefficient in ferromagnet (m^2/s). Default: 2.5e-4.
    model : str, optional
        Equation model: "thin_s" or "fominov". Default: "thin_s".
    phase : str, optional
        "zero" or "pi". Default: "zero".
    depairing : dict, optional
        Depairing channels: {"ag": float, "zeeman": float, "orbital": float,
        "spin_orbit": float}. All default to 0.0.

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
    if Tc0 <= 0:
        raise ValueError(f"Tc0 must be positive, got {Tc0}")
    if d_S <= 0:
        raise ValueError(f"d_S must be positive, got {d_S}")
    if xi_S <= 0:
        raise ValueError(f"xi_S must be positive, got {xi_S}")
    if xi_F <= 0:
        raise ValueError(f"xi_F must be positive, got {xi_F}")
    if model not in ("thin_s", "fominov"):
        raise ValueError(f"model must be 'thin_s' or 'fominov', got {model!r}")
    if phase not in ("zero", "pi"):
        raise ValueError(f"phase must be 'zero' or 'pi', got {phase!r}")

    d_F_array = np.asarray(d_F_array, dtype=np.float64)

    if _USE_NATIVE:
        model_int = 1 if model == "fominov" else 0
        phase_int = 1 if phase == "pi" else 0
        dp = depairing or {}
        return _native_solve_tc_batch(
            Tc0, d_S, xi_S, xi_F, gamma, gamma_B, E_ex, D_F,
            model_int, phase_int,
            dp.get("ag", 0.0), dp.get("zeeman", 0.0),
            dp.get("orbital", 0.0), dp.get("spin_orbit", 0.0),
            d_F_array)

    # Pure Python fallback: digamma-based self-consistency
    from scipy.special import digamma as _digamma

    inv_xi_F = 1.0 / xi_F
    Tc_out = np.empty_like(d_F_array)

    lambda_dep = 0.0
    if depairing:
        lambda_dep = sum(depairing.get(k, 0.0)
                         for k in ("ag", "zeeman", "orbital", "spin_orbit"))

    for i, d_F in enumerate(d_F_array):
        if d_F <= 0:
            Tc_out[i] = Tc0
            continue

        # Compute kernel magnitude (simplified scalar form)  [EQ-2, EQ-3]
        # 0-junction: K = q*coth(q*d_F) — cosh/sinh
        # π-junction: K = q*tanh(q*d_F) — sinh/cosh
        q = (1.0 + 1.0j) / xi_F
        qd = q * d_F
        if phase == "zero":
            K = q * np.cosh(qd) / np.sinh(qd)
        elif phase == "pi":
            K = q * np.sinh(qd) / np.cosh(qd)

        # Effective coupling  [EQ-4, EQ-5]
        # No eta = xi_S/d_S prefactor — gamma absorbs coupling strength.
        if model == "fominov" and gamma_B > 0:
            alpha_K = gamma * K / (1.0 + gamma_B * K)
        else:
            alpha_K = gamma * K

        # Brent-like scan for highest root of
        # F(T) = ln(Tc0/T) - Re[psi(0.5 + alpha*Tc0/(2*pi*T) + lambda_dep) - psi(0.5)]
        best_tc = 0.0
        psi_half = float(_digamma(0.5))
        T_vals = np.linspace(0.01, Tc0, 1000)
        for j in range(len(T_vals) - 1):
            T_a, T_b = T_vals[j], T_vals[j + 1]

            def _f(T):
                A = alpha_K * Tc0 / (2.0 * np.pi * T) + lambda_dep
                psi_val = _digamma(complex(0.5 + A))
                return np.log(Tc0 / T) - (psi_val - psi_half).real

            fa, fb = _f(T_a), _f(T_b)
            if fa * fb < 0:
                # Bisection refinement
                for _ in range(60):
                    T_m = 0.5 * (T_a + T_b)
                    fm = _f(T_m)
                    if fm * fa < 0:
                        T_b = T_m
                    else:
                        T_a = T_m
                        fa = fm
                root = 0.5 * (T_a + T_b)
                if root > best_tc:
                    best_tc = root

        Tc_out[i] = np.clip(best_tc, 0.0, Tc0)

    return Tc_out
