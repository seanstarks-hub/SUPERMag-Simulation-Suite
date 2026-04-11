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

    # Pure Python fallback — Buzdin model for S/F/S CPR
    phi = np.linspace(0, 2 * np.pi, n_phases, endpoint=False)

    # Complex wave vector in F layer: q = (1+i)/ξ_F
    q = (1.0 + 1.0j) / xi_F

    # Critical current kernel (Buzdin, diffusive limit):
    #   I_c ∝ exp(-d_F/ξ_F) cos(d_F/ξ_F − π/4)
    # which equals Re[exp(-q·d_F) · exp(iπ/4)]
    Kc = np.exp(-q * d_F)
    I_c = np.real(Kc * np.exp(1j * np.pi / 4))

    # Temperature factor from BCS gap Δ(T) = 1.764 kB Tc √(1−T/Tc)
    kB_meV = 8.617333262e-2  # meV/K
    Tc_ref = Tc0
    Delta_0 = 1.764 * kB_meV * Tc_ref
    t_ratio = min(T / Tc_ref, 0.9999)
    Delta_T = Delta_0 * np.sqrt(max(1.0 - t_ratio, 0.0))
    T_factor = np.tanh(Delta_T / (2.0 * kB_meV * max(T, 0.01)))

    # First-harmonic CPR: I(φ) = I_c(d_F, T) sin(φ)
    I_raw = I_c * T_factor * np.sin(phi)

    # Normalize to max |I| = 1
    max_I = np.max(np.abs(I_raw))
    I_norm = I_raw / max_I if max_I > 0 else I_raw

    return phi, I_norm
