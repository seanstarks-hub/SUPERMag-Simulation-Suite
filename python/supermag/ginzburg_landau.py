"""
Ginzburg-Landau free energy solver — Python interface.

Minimizes the GL functional on a 2D grid for vortex states, mixed-state
configurations, and domain structures near Tc.

GL free energy density:
    f = α|ψ|² + (β/2)|ψ|⁴ + (ξ²/2)|∇ψ|²

First GL equation (TDGL relaxation):
    ∂ψ/∂t = −δF/δψ* = −αψ − β|ψ|²ψ + ξ²∇²ψ

The equilibrium is found by iterating until convergence.

See docs/theory/ginzburg_landau.ipynb for details.

References:
    Ginzburg, V.L. & Landau, L.D. (1950). Zh. Eksp. Teor. Fiz. 20, 1064.
    Gor'kov, L.P. (1959). Sov. Phys. JETP 9, 1364.
"""

import numpy as np


_USE_NATIVE = False
try:
    from supermag._native import _gl_minimize as _native_gl_minimize
    _USE_NATIVE = True
except ImportError:
    pass


def minimize(alpha, beta, kappa, nx, ny, dx, H_applied=0.0,
             mode="scalar", seed=42):
    """
    Minimize the Ginzburg-Landau free energy on a 2D grid.

    Uses relaxation of the first GL equation:
        ∂ψ/∂t = −αψ − β|ψ|²ψ + ξ²∇²ψ
    where ξ = dx / κ in the simulation units.

    Parameters
    ----------
    alpha : float
        GL alpha parameter (< 0 below Tc).
    beta : float
        GL beta parameter (> 0 for stability).
    kappa : float
        GL kappa = λ_L / ξ.
    nx, ny : int
        Grid dimensions (must be > 0).
    dx : float
        Grid spacing (nm, must be > 0).
    H_applied : float, optional
        External magnetic field (Tesla). Default: 0.0.
    mode : str, optional
        "scalar" (no vector potential) or "gauge" (includes A field).
        Default: "scalar".
    seed : int, optional
        Random seed for initial condition. Default: 42.

    Returns
    -------
    psi : numpy.ndarray
        Complex order parameter on grid, shape (ny, nx).

    Raises
    ------
    ValueError
        If nx <= 0, ny <= 0, dx <= 0, or beta <= 0.
    NotImplementedError
        If mode="gauge" (not yet implemented in Python fallback).
    """
    if nx <= 0 or ny <= 0:
        raise ValueError(f"nx and ny must be > 0, got nx={nx}, ny={ny}")
    if dx <= 0:
        raise ValueError(f"dx must be > 0, got {dx}")
    if beta <= 0:
        raise ValueError(f"beta must be > 0 for stability, got {beta}")

    mode_map = {"scalar": 0, "gauge": 1}
    mode_int = mode_map.get(mode, 0) if isinstance(mode, str) else int(mode)

    if _USE_NATIVE:
        psi_re, psi_im = _native_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                              mode_int, H_applied)
        psi_re = np.asarray(psi_re).reshape(ny, nx)
        psi_im = np.asarray(psi_im).reshape(ny, nx)
        return psi_re + 1j * psi_im

    if mode == "gauge":
        raise NotImplementedError(
            "Gauge mode (with vector potential) is not available in the "
            "Python fallback. Build the C++ extension for gauge support.")

    # Pure Python fallback — TDGL relaxation (no vector potential)
    # Equilibrium |ψ|² = −α/β for uniform state (α < 0)
    psi_eq = np.sqrt(max(-alpha / beta, 0.0)) if beta > 0 else 0.0

    # Coherence length in grid units
    xi = dx  # ξ_GL is effectively dx in our discretization
    xi2 = xi * xi

    # Initialize with bulk value + small random perturbation
    rng = np.random.default_rng(seed)
    psi = psi_eq * (1.0 + 0.05 * (rng.standard_normal((ny, nx))
                                    + 1j * rng.standard_normal((ny, nx))))

    # Time step for stability: dt < dx²/(4·ξ²) for 2D diffusion
    dt = 0.1 * dx * dx / (4.0 * xi2 + 1e-15)
    n_steps = 2000

    for _step in range(n_steps):
        # Laplacian with periodic boundary conditions
        lap = (np.roll(psi, 1, axis=0) + np.roll(psi, -1, axis=0)
               + np.roll(psi, 1, axis=1) + np.roll(psi, -1, axis=1)
               - 4.0 * psi) / (dx * dx)

        # TDGL: dψ/dt = −αψ − β|ψ|²ψ + ξ²∇²ψ
        dpsi = -alpha * psi - beta * np.abs(psi)**2 * psi + xi2 * lap

        psi += dt * dpsi

        # Check convergence every 100 steps
        if _step % 100 == 99:
            residual = np.max(np.abs(dpsi))
            if residual < 1e-8:
                break

    return psi
