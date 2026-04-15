"""
Usadel diffusive-limit solver — Python interface.

Solves the linearized Usadel equation for quasiclassical Green's functions
in the diffusive limit with Kupriyanov–Lukichev boundary conditions at
S/F interfaces.

Theta parameterization: g = cos(θ), f = sin(θ), with normalization g² + ff̄ = 1.

Usadel equation in F layer:
    D_F θ'' + 2i(ω_n + iE_ex) sin(θ) = 0

Linearized (θ ≪ 1):
    D_F θ'' + 2i(ω_n + iE_ex) θ = 0

Solution: θ(x) ∝ exp(−κx), κ = √(2(ω_n + iE_ex)/D_F)

See docs/theory/usadel_equations.ipynb for full mathematical treatment.

References:
    Usadel, K.D. (1970). Phys. Rev. Lett. 25, 507.
    Kupriyanov, M.Yu. & Lukichev, V.F. (1988). Sov. Phys. JETP 67, 1163.
"""

import numpy as np


_USE_NATIVE = False
try:
    from supermag._native import _usadel_solve as _native_usadel_solve
    _USE_NATIVE = True
except ImportError:
    pass


def solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex, n_grid=200, T=None):
    """
    Solve the Usadel equation for an S/F bilayer.

    Computes the self-consistent superconducting order parameter Δ(x)
    across the S/F heterostructure using a linearized approach with
    Matsubara frequency summation.

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
    T : float, optional
        Temperature (K). If None, defaults to 0.5 * Tc0.

    Returns
    -------
    x : numpy.ndarray
        Position array (nm), spanning [−d_S, d_F].
    Delta : numpy.ndarray
        Order parameter profile (meV), real.
    """
    if _USE_NATIVE:
        T_val = T if T is not None else 0.5 * Tc0
        return _native_usadel_solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex, T_val, n_grid)

    # Pure Python fallback — linearized Usadel with Matsubara sum
    kB_meV = 8.617333262e-2  # meV/K

    # BCS gap at T=0
    Delta_0 = 1.764 * kB_meV * Tc0  # meV

    # Diffusion constants from coherence lengths
    # ξ = sqrt(ℏD / 2πkBTc) → D = ξ² · 2π·kB·Tc / ℏ
    # For the solver we work in reduced units: D_S = ξ_S², D_F = ξ_F² · E_ex
    # (setting ℏ=1 and absorbing constants into κ)

    # Spatial grid: S region [-d_S, 0] and F region [0, d_F]
    n_S = n_grid // 2
    n_F = n_grid - n_S
    x_S = np.linspace(-d_S, 0, n_S, endpoint=False)
    x_F = np.linspace(0, d_F, n_F)
    x = np.concatenate([x_S, x_F])
    h_S = d_S / max(n_S - 1, 1)
    h_F = d_F / max(n_F - 1, 1)

    # Temperature: use specified T or default to 0.5 * Tc0
    T_use = T if T is not None else 0.5 * Tc0
    Delta_T = Delta_0 * np.sqrt(1.0 - T_use / Tc0)

    # Number of Matsubara frequencies
    n_matsubara = 100

    # Initialize order parameter
    Delta = np.zeros(n_grid)
    # Bulk BCS value in S, zero in F
    Delta[:n_S] = Delta_T
    Delta[n_S:] = 0.0

    # Self-consistency iteration
    for _iteration in range(30):
        theta_sum = np.zeros(n_grid, dtype=complex)

        for n_m in range(n_matsubara):
            omega_n = np.pi * kB_meV * T_use * (2 * n_m + 1)  # meV

            # Solve linearized Usadel on S and F grids separately
            # then match at the interface

            # --- S region: D_S θ'' + 2ω_n sin(θ) ≈ 0 for linearized ---
            # κ_S = sqrt(2ω_n) / ξ_S (in units where D_S = ξ_S²·E_scale)
            kappa_S = np.sqrt(2.0 * omega_n / (kB_meV * Tc0)) / xi_S
            # BCS bulk value of θ
            theta_BCS = np.arctan(Delta_T / omega_n)

            # θ_S(x) = θ_BCS + (θ_interface - θ_BCS) · cosh(κ_S(x+d_S)) / cosh(κ_S·d_S)
            # (Neumann BC at x = −d_S, matching at x = 0)

            # --- F region: θ'' + 2(ω_n + iE_ex)/(ξ_F²·E_ex) · θ = 0 ---
            kappa_F = (1.0 / xi_F) * np.sqrt(2.0 * (omega_n + 1j * E_ex) / E_ex)
            # θ_F(x) = θ_interface · exp(−κ_F · x)  (decays into F)

            # Interface matching (simplified KL BC):
            # σ_S ∂θ_S/∂x |₀ = σ_F ∂θ_F/∂x |₀
            # θ_interface ≈ θ_BCS · kappa_S / (kappa_S + kappa_F · ξ_S/ξ_F)
            # Simplified for transparency γ ~ 1:
            theta_int = theta_BCS * kappa_S * xi_S / (kappa_S * xi_S + kappa_F * xi_F)

            # Build θ(x) on the grid
            theta_x = np.zeros(n_grid, dtype=complex)

            # S region
            kS_dS = kappa_S * d_S
            for j in range(n_S):
                xj = x_S[j]
                if kS_dS < 50:
                    ratio = np.cosh(kappa_S * (xj + d_S)) / np.cosh(kS_dS)
                else:
                    ratio = np.exp(kappa_S * (xj + d_S) - kS_dS)
                theta_x[j] = theta_BCS + (theta_int - theta_BCS) * ratio

            # F region
            for j in range(n_F):
                xj = x_F[j]
                kF_x = kappa_F * xj
                if np.real(kF_x) < 50:
                    theta_x[n_S + j] = theta_int * np.exp(-kF_x)
                else:
                    theta_x[n_S + j] = 0.0

            theta_sum += theta_x

        # Self-consistency: Δ(x) = λ · πT · Σ sin(θ)
        # In linearized form: Δ ∝ πT · Σ Re[θ(ω_n, x)]
        # Use BCS coupling constant via Tc relation
        # 1/λ = πT Σ 1/sqrt(ω_n² + Δ²) → λ ≈ 0.3 (typical for Nb)
        lambda_BCS = 0.3
        Delta_new = lambda_BCS * np.pi * kB_meV * T_use * np.real(theta_sum)

        # Only update S region (Δ = 0 in F by definition for non-SC)
        Delta_new[n_S:] = 0.0

        # Check convergence
        change = np.max(np.abs(Delta_new[:n_S] - Delta[:n_S]))
        # Mixing
        Delta = 0.3 * Delta_new + 0.7 * Delta
        Delta[n_S:] = 0.0

        if change < 1e-6:
            break

    return x, np.abs(Delta)
