"""
Eilenberger clean-limit solver — Python interface.

Solves the Eilenberger equation using Riccati parameterization
for quasiclassical Green's functions in the clean limit.

Eilenberger equation (1D, Matsubara):
    ℏv_{F,x} ∂f/∂x = −2(ω_n + iE_ex)f + 2Δ(x)g

Riccati parameterization (Schopohl & Maki, 1995):
    f = 2a/(1 + aā), g = (1 − aā)/(1 + aā)

Riccati ODE:
    ℏv_{F,x} ∂a/∂x = −2(ω_n + iE_ex)a + Δ(x)(1 − a²)

See docs/theory/usadel_equations.ipynb for the diffusive counterpart.

References:
    Eilenberger, G. (1968). Z. Phys. 214, 195.
    Schopohl, N. & Maki, K. (1995). Phys. Rev. B 52, 490.
"""

import numpy as np


_USE_NATIVE = False
try:
    from supermag._native import _eilenberger_solve as _native_eilenberger_solve
    _USE_NATIVE = True
except ImportError:
    pass


def solve(Tc0, d_S, d_F, xi_S, E_ex, n_grid=200):
    """
    Solve the Eilenberger equation for an S/F bilayer (clean limit).

    Integrates the Riccati ODE for the anomalous Green's function
    through S and F layers with Fermi-surface averaging.

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
        Position array (nm), spanning [−d_S, d_F].
    f : numpy.ndarray
        Anomalous Green's function |f(x)| (Fermi-surface averaged).
    """
    if _USE_NATIVE:
        return _native_eilenberger_solve(Tc0, d_S, d_F, xi_S, E_ex, n_grid)

    # Pure Python fallback — Riccati integration with angular averaging
    kB_meV = 8.617333262e-2  # meV/K

    # BCS gap
    Delta_0 = 1.764 * kB_meV * Tc0
    T = 0.5 * Tc0
    Delta_T = Delta_0 * np.sqrt(max(1.0 - T / Tc0, 0.0))

    # Fermi velocity scale: ℏv_F / kB = hbar_vF_nm_meV
    # ξ_S = ℏv_F / (2πkBTc) → ℏv_F = 2π·kB·Tc · ξ_S
    hbar_vF = 2.0 * np.pi * kB_meV * Tc0 * xi_S  # nm·meV

    # Spatial grid spanning S and F layers
    n_S = n_grid // 2
    n_F = n_grid - n_S
    x_S = np.linspace(-d_S, 0, n_S, endpoint=False)
    x_F = np.linspace(0, d_F, n_F)
    x = np.concatenate([x_S, x_F])
    dx_val = (d_S + d_F) / max(n_grid - 1, 1)

    # Matsubara frequency (first, dominant)
    omega_1 = np.pi * kB_meV * T  # lowest Matsubara freq (meV)

    # BCS bulk value of Riccati parameter
    a_BCS = Delta_T / (omega_1 + np.sqrt(omega_1**2 + Delta_T**2))

    # Gauss-Legendre quadrature for angular averaging
    n_angles = 16
    from numpy.polynomial.legendre import leggauss
    nodes, weights = leggauss(n_angles)
    # Map from [-1, 1] to cos(θ) ∈ [-1, 1]
    cos_theta = nodes
    # Only right-movers (cos_θ > 0) for forward integration
    # and left-movers for backward — average over all

    # Profiles
    def delta_profile(xp):
        return Delta_T if xp < 0 else 0.0

    def eex_profile(xp):
        return 0.0 if xp < 0 else E_ex

    # Integrate Riccati ODE for each angle
    f_avg = np.zeros(n_grid)

    for k in range(n_angles):
        cos_k = cos_theta[k]
        if abs(cos_k) < 1e-10:
            continue

        vx = hbar_vF * cos_k  # ℏv_{F,x} in nm·meV

        # Riccati ODE: vx · da/dx = -2(ω + iE_ex)a + Δ(1 - a²)
        a = np.zeros(n_grid, dtype=complex)

        if cos_k > 0:
            # Right-mover: integrate left → right
            a[0] = a_BCS + 0j
            for j in range(n_grid - 1):
                xj = x[j]
                D = delta_profile(xj)
                E = eex_profile(xj)
                da = (-2.0 * (omega_1 + 1j * E) * a[j]
                      + D * (1.0 - a[j]**2)) / vx
                a[j + 1] = a[j] + da * dx_val
                # Clamp for stability
                if abs(a[j + 1]) > 2.0:
                    a[j + 1] *= 2.0 / abs(a[j + 1])
        else:
            # Left-mover: integrate right → left
            a[-1] = 0.0 + 0j  # normal-metal boundary
            for j in range(n_grid - 1, 0, -1):
                xj = x[j]
                D = delta_profile(xj)
                E = eex_profile(xj)
                da = (-2.0 * (omega_1 + 1j * E) * a[j]
                      + D * (1.0 - a[j]**2)) / vx
                a[j - 1] = a[j] - da * dx_val
                if abs(a[j - 1]) > 2.0:
                    a[j - 1] *= 2.0 / abs(a[j - 1])

        # f = 2a / (1 + |a|²)
        f_k = 2.0 * a / (1.0 + np.abs(a)**2)
        f_avg += weights[k] * np.abs(f_k)

    # Normalize by total weight (integral of sin(θ)dθ over [0,π] = 2)
    f_avg /= np.sum(weights)

    return x, f_avg
