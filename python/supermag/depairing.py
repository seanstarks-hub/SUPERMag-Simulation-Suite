"""
Depairing channel computations and optimizer utilities.

Computes dimensionless pair-breaking parameters from physical laboratory
inputs (scattering rates, magnetic fields, diffusion coefficients),
and provides optimization/fitting tools for matching Tc(d_F) data.

Physics:
    Pair-breaking (depairing) channels suppress the superconducting order
    parameter. The total depairing parameter enters the self-consistency
    equation as an additive shift:
        λ_dep = λ_AG + λ_Z + λ_orb + λ_SO       (EQ-7)

    Individual channels:
        λ_AG  = Γ_s / (2·kB·T)                   (EQ-7A, spin-flip)
        λ_Z   = (μ_B·H)² / (2π·kB·T)²           (EQ-7B, Zeeman)
        λ_orb = D·(eH)²·d² / (3ℏ²·2π·kB·T)      (EQ-7C, orbital)
        λ_SO  = Γ_so / (2·kB·T)                  (EQ-7D, spin-orbit)

References:
    Abrikosov, A.A. & Gor'kov, L.P. (1961). JETP 12, 1243.
    Maki, K. (1966). Prog. Theor. Phys. 36, 1.
"""

import numpy as np


# Try to import the native C++ extension; fall back to pure Python
_USE_NATIVE = False
try:
    from supermag._native import (
        _depairing_ag as _native_depairing_ag,
        _depairing_zeeman as _native_depairing_zeeman,
        _depairing_orbital_perp as _native_depairing_orbital_perp,
        _depairing_orbital_par as _native_depairing_orbital_par,
        _depairing_soc as _native_depairing_soc,
        _depairing_from_physical as _native_depairing_from_physical,
        _optimize_tc as _native_optimize_tc,
        _inverse_tc as _native_inverse_tc,
        _fit_tc as _native_fit_tc,
    )
    _USE_NATIVE = True
except ImportError:
    pass


# ── Physical constants (CODATA 2018) ────────────────────────────

_HBAR = 1.054571817e-34    # J·s
_KB   = 1.380649e-23       # J/K
_MU_B = 9.2740100783e-24   # J/T
_E    = 1.602176634e-19    # C
_MEV_TO_J = 1.602176634e-22
_NM2PS_TO_M2S = 1.0e-6
_NM_TO_M = 1.0e-9


# ── Individual depairing channels ───────────────────────────────

def depairing_ag(gamma_s_meV, T_kelvin):
    """
    Abrikosov-Gor'kov pair-breaking from spin-flip scattering.

    λ_AG = Γ_s / (2·kB·T)

    Parameters
    ----------
    gamma_s_meV : float
        Spin-flip scattering rate (meV).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    float
        Dimensionless AG depairing parameter.
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        return _native_depairing_ag(gamma_s_meV, T_kelvin)
    Gamma_s_SI = gamma_s_meV * _MEV_TO_J
    return Gamma_s_SI / (2.0 * _KB * T_kelvin)


def depairing_zeeman(H_tesla, T_kelvin):
    """
    Zeeman (Pauli paramagnetic) pair-breaking.

    λ_Z = (μ_B·H)² / (2π·kB·T)²

    Parameters
    ----------
    H_tesla : float
        Applied magnetic field (T).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    float
        Dimensionless Zeeman depairing parameter.
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        return _native_depairing_zeeman(H_tesla, T_kelvin)
    muBH = _MU_B * H_tesla
    denom = 2.0 * np.pi * _KB * T_kelvin
    return (muBH * muBH) / (denom * denom)


def depairing_orbital_perp(D_nm2ps, H_tesla, thickness_nm, T_kelvin):
    """
    Orbital pair-breaking for perpendicular applied field.

    λ_orb⊥ = D·(eH)²·d² / (3ℏ²·2π·kB·T)

    Parameters
    ----------
    D_nm2ps : float
        Diffusion coefficient (nm²/ps).
    H_tesla : float
        Applied magnetic field (T).
    thickness_nm : float
        Film thickness (nm).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    float
        Dimensionless orbital depairing parameter (perp. field).
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        return _native_depairing_orbital_perp(D_nm2ps, H_tesla,
                                               thickness_nm, T_kelvin)
    D_SI = D_nm2ps * _NM2PS_TO_M2S
    d_SI = thickness_nm * _NM_TO_M
    eH = _E * H_tesla
    return D_SI * eH * eH * d_SI * d_SI / (
        3.0 * _HBAR * _HBAR * 2.0 * np.pi * _KB * T_kelvin)


def depairing_orbital_par(D_nm2ps, H_tesla, thickness_nm, T_kelvin):
    """
    Orbital pair-breaking for parallel applied field.

    λ_orb∥ = D·(eH)²·d² / (12ℏ²·2π·kB·T)

    Parameters
    ----------
    D_nm2ps : float
        Diffusion coefficient (nm²/ps).
    H_tesla : float
        Applied magnetic field (T).
    thickness_nm : float
        Film thickness (nm).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    float
        Dimensionless orbital depairing parameter (parallel field).
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        return _native_depairing_orbital_par(D_nm2ps, H_tesla,
                                              thickness_nm, T_kelvin)
    D_SI = D_nm2ps * _NM2PS_TO_M2S
    d_SI = thickness_nm * _NM_TO_M
    eH = _E * H_tesla
    return D_SI * eH * eH * d_SI * d_SI / (
        12.0 * _HBAR * _HBAR * 2.0 * np.pi * _KB * T_kelvin)


def depairing_soc(Gamma_so_meV, T_kelvin):
    """
    Spin-orbit coupling pair-breaking.

    λ_SO = Γ_so / (2·kB·T)

    Parameters
    ----------
    Gamma_so_meV : float
        Spin-orbit scattering rate (meV).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    float
        Dimensionless spin-orbit depairing parameter.
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        return _native_depairing_soc(Gamma_so_meV, T_kelvin)
    Gamma_so_SI = Gamma_so_meV * _MEV_TO_J
    return Gamma_so_SI / (2.0 * _KB * T_kelvin)


def depairing_from_physical(gamma_s_meV, H_tesla, D_nm2ps,
                             thickness_nm, Gamma_so_meV, T_kelvin):
    """
    Compute all depairing channels from physical laboratory inputs.

    Parameters
    ----------
    gamma_s_meV : float
        Spin-flip scattering rate (meV).
    H_tesla : float
        Applied magnetic field (T).
    D_nm2ps : float
        Diffusion coefficient (nm²/ps).
    thickness_nm : float
        Film thickness (nm).
    Gamma_so_meV : float
        Spin-orbit scattering rate (meV).
    T_kelvin : float
        Temperature (K). Must be > 0.

    Returns
    -------
    dict
        ``{"ag": float, "zeeman": float, "orbital": float, "spin_orbit": float}``
    """
    if T_kelvin <= 0:
        raise ValueError(f"T_kelvin must be positive, got {T_kelvin}")
    if _USE_NATIVE:
        ag, zeeman, orbital, spin_orbit = _native_depairing_from_physical(
            gamma_s_meV, H_tesla, D_nm2ps, thickness_nm,
            Gamma_so_meV, T_kelvin)
        return {"ag": ag, "zeeman": zeeman,
                "orbital": orbital, "spin_orbit": spin_orbit}

    return {
        "ag": depairing_ag(gamma_s_meV, T_kelvin),
        "zeeman": depairing_zeeman(H_tesla, T_kelvin),
        "orbital": depairing_orbital_perp(D_nm2ps, H_tesla,
                                           thickness_nm, T_kelvin),
        "spin_orbit": depairing_soc(Gamma_so_meV, T_kelvin),
    }


# ── Optimizer utilities ─────────────────────────────────────────

def optimize_tc(Tc0, d_S, xi_S, xi_F, E_ex, gamma=0.3, gamma_B=0.0,
                D_F=2.5e-4, model="thin_s", phase="zero",
                d_F_lo=0.5, d_F_hi=50.0, Tc_target=4.2,
                depairing=None):
    """
    Find the F-layer thickness that produces a target Tc.

    Uses golden-section search to minimize |Tc(d_F) - Tc_target|.

    Parameters
    ----------
    Tc0 : float
        Bulk superconductor Tc (K).
    d_S : float
        Superconductor thickness (nm).
    xi_S, xi_F : float
        Coherence lengths (nm).
    E_ex : float
        Exchange energy (meV).
    gamma, gamma_B : float
        Interface parameters.
    D_F : float
        Diffusion coefficient (m²/s).
    model, phase : str
        Equation model and junction phase.
    d_F_lo, d_F_hi : float
        Search bounds for d_F (nm).
    Tc_target : float
        Target critical temperature (K).
    depairing : dict, optional
        Depairing channels.

    Returns
    -------
    float
        Optimal d_F (nm) that gives Tc closest to Tc_target.
    """
    if Tc_target <= 0:
        raise ValueError(f"Tc_target must be positive, got {Tc_target}")
    if d_F_lo >= d_F_hi:
        raise ValueError(f"d_F_lo must be < d_F_hi, got {d_F_lo} >= {d_F_hi}")

    if _USE_NATIVE:
        model_int = 1 if model == "fominov" else 0
        phase_int = 1 if phase == "pi" else 0
        dp = depairing or {}
        return _native_optimize_tc(
            Tc0, d_S, xi_S, xi_F, gamma, gamma_B, E_ex, D_F,
            model_int, phase_int,
            dp.get("ag", 0.0), dp.get("zeeman", 0.0),
            dp.get("orbital", 0.0), dp.get("spin_orbit", 0.0),
            d_F_lo, d_F_hi, Tc_target)

    # Pure Python fallback: golden-section search
    from supermag.proximity import critical_temperature

    gr = (np.sqrt(5) + 1) / 2
    a, b = d_F_lo, d_F_hi

    dp = depairing or {}

    def _tc_at(d_F_val):
        return critical_temperature(
            Tc0=Tc0, d_S=d_S, d_F_array=np.array([d_F_val]),
            E_ex=E_ex, xi_S=xi_S, xi_F=xi_F,
            gamma=gamma, gamma_B=gamma_B, D_F=D_F,
            model=model, phase=phase, depairing=dp)[0]

    for _ in range(100):
        if abs(b - a) < 1e-10:
            break
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        fc = abs(_tc_at(c) - Tc_target)
        fd = abs(_tc_at(d) - Tc_target)
        if fc < fd:
            b = d
        else:
            a = c

    return (a + b) / 2


def inverse_tc(Tc0, d_S, xi_S, xi_F, E_ex, gamma=0.3, gamma_B=0.0,
               D_F=2.5e-4, model="thin_s", phase="zero",
               Tc_target=4.2, d_F_lo=0.5, d_F_hi=50.0,
               depairing=None):
    """
    Find d_F that produces exactly Tc_target via Brent's method.

    Unlike ``optimize_tc`` (which minimizes |Tc - target|), this finds
    the exact root of Tc(d_F) - Tc_target = 0.

    Parameters
    ----------
    Tc0, d_S, xi_S, xi_F, E_ex, gamma, gamma_B, D_F, model, phase
        Same as ``optimize_tc``.
    Tc_target : float
        Target Tc (K).
    d_F_lo, d_F_hi : float
        Bracket for Brent search (nm).
    depairing : dict, optional
        Depairing channels.

    Returns
    -------
    float
        d_F (nm) that gives Tc = Tc_target, or NaN if no root in bracket.
    """
    if Tc_target <= 0:
        raise ValueError(f"Tc_target must be positive, got {Tc_target}")
    if d_F_lo >= d_F_hi:
        raise ValueError(f"d_F_lo must be < d_F_hi, got {d_F_lo} >= {d_F_hi}")

    if _USE_NATIVE:
        model_int = 1 if model == "fominov" else 0
        phase_int = 1 if phase == "pi" else 0
        dp_dict = depairing or {}
        return _native_inverse_tc(
            Tc0, d_S, xi_S, xi_F, gamma, gamma_B, E_ex, D_F,
            model_int, phase_int,
            dp_dict.get("ag", 0.0), dp_dict.get("zeeman", 0.0),
            dp_dict.get("orbital", 0.0), dp_dict.get("spin_orbit", 0.0),
            Tc_target, d_F_lo, d_F_hi)

    # Pure Python fallback: bisection
    from supermag.proximity import critical_temperature
    dp = depairing or {}

    def _f(d_F_val):
        tc = critical_temperature(
            Tc0=Tc0, d_S=d_S, d_F_array=np.array([d_F_val]),
            E_ex=E_ex, xi_S=xi_S, xi_F=xi_F,
            gamma=gamma, gamma_B=gamma_B, D_F=D_F,
            model=model, phase=phase, depairing=dp)[0]
        return tc - Tc_target

    a, b = d_F_lo, d_F_hi
    fa, fb = _f(a), _f(b)
    if fa * fb > 0:
        return (a if abs(fa) < abs(fb) else b)

    for _ in range(100):
        c = (a + b) / 2
        fc = _f(c)
        if abs(fc) < 1e-12 or (b - a) < 1e-12:
            break
        if fa * fc < 0:
            b = c
        else:
            a, fa = c, fc

    return (a + b) / 2


def fit_tc(Tc0, d_S, xi_S, xi_F, E_ex, gamma=0.3, gamma_B=0.0,
           D_F=2.5e-4, model="thin_s", phase="zero",
           d_F_data=None, Tc_data=None,
           fit_gamma=True, fit_gamma_B=False,
           fit_E_ex=False, fit_xi_F=False,
           depairing=None):
    """
    Fit proximity parameters to experimental Tc(d_F) data.

    Uses Nelder-Mead to minimize chi² between computed and measured Tc.

    Parameters
    ----------
    Tc0, d_S, xi_S, xi_F, E_ex, gamma, gamma_B, D_F, model, phase
        Initial parameter values. Fitted parameters are overwritten.
    d_F_data : array_like
        Experimental F-layer thicknesses (nm), length N >= 2.
    Tc_data : array_like
        Experimental Tc values (K), length N.
    fit_gamma : bool
        Include gamma in fit. Default: True.
    fit_gamma_B : bool
        Include gamma_B in fit. Default: False.
    fit_E_ex : bool
        Include E_ex in fit. Default: False.
    fit_xi_F : bool
        Include xi_F in fit. Default: False.
    depairing : dict, optional
        Depairing channels.

    Returns
    -------
    dict
        ``{"gamma": float, "gamma_B": float, "E_ex": float,
          "xi_F": float, "chi2": float}`` — fitted values + final chi².
    """
    if d_F_data is None or Tc_data is None:
        raise ValueError("d_F_data and Tc_data are required")
    d_F_data = np.asarray(d_F_data, dtype=np.float64)
    Tc_data = np.asarray(Tc_data, dtype=np.float64)
    if len(d_F_data) < 2:
        raise ValueError(f"Need at least 2 data points, got {len(d_F_data)}")
    if len(d_F_data) != len(Tc_data):
        raise ValueError("d_F_data and Tc_data must have same length")

    if _USE_NATIVE:
        model_int = 1 if model == "fominov" else 0
        phase_int = 1 if phase == "pi" else 0
        dp = depairing or {}
        result = _native_fit_tc(
            Tc0, d_S, xi_S, xi_F, gamma, gamma_B, E_ex, D_F,
            model_int, phase_int,
            dp.get("ag", 0.0), dp.get("zeeman", 0.0),
            dp.get("orbital", 0.0), dp.get("spin_orbit", 0.0),
            d_F_data, Tc_data,
            int(fit_gamma), int(fit_gamma_B),
            int(fit_E_ex), int(fit_xi_F))
        return result

    # Pure Python fallback: scipy Nelder-Mead
    from scipy.optimize import minimize as _minimize
    from supermag.proximity import critical_temperature

    dp = depairing or {}

    # Build initial parameter vector
    param_names = []
    x0 = []
    if fit_gamma:
        param_names.append("gamma")
        x0.append(gamma)
    if fit_gamma_B:
        param_names.append("gamma_B")
        x0.append(gamma_B)
    if fit_E_ex:
        param_names.append("E_ex")
        x0.append(E_ex)
    if fit_xi_F:
        param_names.append("xi_F")
        x0.append(xi_F)

    if not param_names:
        raise ValueError("At least one fit flag must be True")

    base = {"Tc0": Tc0, "d_S": d_S, "xi_S": xi_S, "xi_F": xi_F,
            "E_ex": E_ex, "gamma": gamma, "gamma_B": gamma_B,
            "D_F": D_F, "model": model, "phase": phase, "depairing": dp}

    def _chi2(params):
        kw = dict(base)
        for name, val in zip(param_names, params):
            if val <= 0:
                return 1e20
            kw[name] = val
        Tc_calc = critical_temperature(d_F_array=d_F_data, **kw)
        return float(np.sum((Tc_calc - Tc_data) ** 2))

    res = _minimize(_chi2, x0, method="Nelder-Mead",
                    options={"maxiter": 2000, "xatol": 1e-8, "fatol": 1e-10})

    fitted = dict(base)
    for name, val in zip(param_names, res.x):
        fitted[name] = float(val)

    return {
        "gamma": fitted["gamma"],
        "gamma_B": fitted["gamma_B"],
        "E_ex": fitted["E_ex"],
        "xi_F": fitted["xi_F"],
        "chi2": float(res.fun),
    }
