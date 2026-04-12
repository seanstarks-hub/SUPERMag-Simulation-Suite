"""
Josephson junction CPR solver — Python interface.

Computes current-phase relations for S/F/S junctions.
Detects 0-π transitions as function of F-layer thickness and temperature.

The critical current oscillates as (Buzdin, diffusive limit):

    I_c(d_F) ∝ exp(-d_F/ξ_F) cos(d_F/ξ_F − π/4)

with first 0–π transition at d_F ≈ (3π/4)ξ_F.

The full CPR sums over Matsubara frequencies ω_n = πT(2n+1):

    I(φ) = (2πT) Σ_n Re[ Δ²sin(φ) · K_F(ω_n) / D_n(φ) ]

See docs/theory/josephson_cpr.ipynb for details.

References:
    Buzdin, A.I. (2005). Rev. Mod. Phys. 77, 935.
    Ryazanov, V.V. et al. (2001). Phys. Rev. Lett. 86, 2427.
"""

import numpy as np


_USE_NATIVE = False
try:
    from supermag._native import _josephson_cpr as _native_josephson_cpr
    _USE_NATIVE = True
except ImportError:
    pass


def current_phase_relation(d_F, xi_F, E_ex, T, n_phases=100, Tc0=9.2):
    """
    Compute Josephson current-phase relation I(φ) for an S/F/S junction.

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
        Number of phase points from 0 to 2π.
    Tc0 : float, optional
        Bulk superconductor critical temperature (K). Default: 9.2 (Nb).

    Returns
    -------
    phi : numpy.ndarray
        Phase difference array (rad), shape (n_phases,).
    I : numpy.ndarray
        Supercurrent array (normalized to max |I|=1), shape (n_phases,).
    """
    if _USE_NATIVE:
        return _native_josephson_cpr(d_F, xi_F, E_ex, T, Tc0, n_phases)

    # Pure Python fallback — Matsubara frequency sum for S/F/S CPR  [EQ-9]
    phi = np.linspace(0, 2 * np.pi, n_phases, endpoint=False)

    # BCS gap Δ(T) = 1.764·kB·Tc0·√(1−T/Tc0)
    kB_meV = 8.617333262e-2  # meV/K
    Tc_ref = Tc0
    Delta_0 = 1.764 * kB_meV * Tc_ref
    t_ratio = min(T / Tc_ref, 0.9999)
    Delta = Delta_0 * np.sqrt(max(1.0 - t_ratio, 0.0))
    Delta2 = Delta * Delta

    # Temperature in meV for Matsubara sum
    T_meV = kB_meV * max(T, 0.01)

    # Matsubara frequency sum:
    #   I(φ) = T_meV · Σ_n Re[ P_n · Δ²·sin(φ) /
    #          √((ω_n² + Δ²·sin²(φ/2)) · (ω_n² + Δ²)) ]
    #   P_n = exp(-q_n·d_F) · exp(iπ/4)
    #   q_n = √(2·(ω_n/E_ex + i)) / ξ_F
    #   ω_n = π·T_meV·(2n+1)
    N_MAX = 500
    omega_cut = 20.0 * Delta
    phase_rot = np.exp(1j * np.pi / 4.0)

    sin_phi = np.sin(phi)
    sin_half = np.sin(phi / 2.0)
    sin_half2 = sin_half * sin_half

    I_raw = np.zeros(n_phases)

    for n in range(N_MAX + 1):
        omega_n = np.pi * T_meV * (2.0 * n + 1.0)
        if omega_n > omega_cut and n > 0:
            break

        omega_n2 = omega_n * omega_n

        # Complex wave vector at this Matsubara frequency  [EQ-1 generalized]
        q_n = np.sqrt(2.0 * (omega_n / E_ex + 1.0j)) / xi_F

        # F-layer propagator
        P_n = np.exp(-q_n * d_F) * phase_rot

        # Denominator: √(ω_n² + Δ²·sin²(φ/2)) · √(ω_n² + Δ²)
        denom_phi = np.sqrt(omega_n2 + Delta2 * sin_half2)
        denom_const = np.sqrt(omega_n2 + Delta2)

        # Accumulate: Re[ P_n · Δ²·sin(φ) / (denom_phi · denom_const) ]
        I_raw += np.real(P_n * Delta2 * sin_phi / (denom_phi * denom_const))

    # Multiply by T_meV (Matsubara prefactor)
    I_raw *= T_meV

    # Normalize to max |I| = 1
    max_I = np.max(np.abs(I_raw))
    I_norm = I_raw / max_I if max_I > 0 else I_raw

    return phi, I_norm
