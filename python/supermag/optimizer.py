"""
Multi-parameter design optimization for S/F heterostructures.

Provides Nelder-Mead optimizers that compose :func:`critical_temperature` to
find layer thicknesses and interface parameters meeting a target Tc, with
optional sensitivity analysis and robustness penalties.
"""

import numpy as np
from scipy.optimize import minimize as _minimize

from supermag.proximity import critical_temperature


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _tc_scalar(d_S, d_F, gamma, base_kw):
    """Evaluate Tc at a single (d_S, d_F, gamma) point."""
    Tc_arr = critical_temperature(
        d_S=d_S, d_F_array=np.array([d_F]), gamma=gamma, **base_kw
    )
    return float(Tc_arr[0])


# ---------------------------------------------------------------------------
# sensitivity_at
# ---------------------------------------------------------------------------

def sensitivity_at(Tc0, d_S, d_F, E_ex, xi_S, xi_F,
                   gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
                   model="thin_s", phase="zero", depairing=None,
                   delta=0.01):
    """
    Compute partial derivatives of Tc with respect to d_F, d_S, and gamma.

    Uses central finite differences with fractional step *delta*.

    Parameters
    ----------
    Tc0, d_S, d_F, E_ex, xi_S, xi_F, gamma, gamma_B, D_F, model, phase, depairing
        Operating point.  *d_F* is a scalar thickness (nm).
    delta : float, optional
        Fractional step size for finite differences. Default 0.01 (1 %).

    Returns
    -------
    dict
        ``{"dTc_ddF": float, "dTc_ddS": float, "dTc_dgamma": float,
          "Tc": float}``
    """
    base = dict(Tc0=Tc0, E_ex=E_ex, xi_S=xi_S, xi_F=xi_F,
                gamma_B=gamma_B, D_F=D_F, model=model, phase=phase,
                depairing=depairing or {})

    Tc0_val = _tc_scalar(d_S, d_F, gamma, base)

    def _deriv(param, val):
        h = max(abs(val) * delta, 1e-6)
        if param == "d_F":
            tp = _tc_scalar(d_S, val + h, gamma, base)
            tm = _tc_scalar(d_S, val - h, gamma, base)
        elif param == "d_S":
            tp = _tc_scalar(val + h, d_F, gamma, base)
            tm = _tc_scalar(val - h, d_F, gamma, base)
        else:  # gamma
            tp = _tc_scalar(d_S, d_F, val + h, base)
            tm = _tc_scalar(d_S, d_F, val - h, base)
        return (tp - tm) / (2 * h)

    return {
        "dTc_ddF": _deriv("d_F", d_F),
        "dTc_ddS": _deriv("d_S", d_S),
        "dTc_dgamma": _deriv("gamma", gamma),
        "Tc": Tc0_val,
    }


# ---------------------------------------------------------------------------
# optimize_design
# ---------------------------------------------------------------------------

def optimize_design(Tc0, d_S, xi_S, xi_F, E_ex,
                    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
                    model="thin_s", phase="zero", depairing=None,
                    Tc_target=4.2,
                    vary_d_S=None, vary_d_F=None, vary_gamma=None,
                    max_total_thickness=None):
    """
    Find layer dimensions and interface parameter to achieve a target Tc.

    Uses Nelder-Mead optimisation over a user-selected subset of ``(d_S,
    d_F, gamma)``.  Each ``vary_*`` argument is either ``None`` (fixed) or a
    ``(lo, hi)`` bounds tuple.

    Parameters
    ----------
    Tc0, d_S, xi_S, xi_F, E_ex, gamma, gamma_B, D_F, model, phase, depairing
        Baseline material / geometry parameters.
    Tc_target : float
        Desired critical temperature (K).
    vary_d_S : tuple or None
        ``(lo, hi)`` bounds for d_S optimisation (nm), or *None* to fix.
    vary_d_F : tuple or None
        ``(lo, hi)`` bounds for d_F optimisation (nm), or *None* to fix.
    vary_gamma : tuple or None
        ``(lo, hi)`` bounds for gamma optimisation, or *None* to fix.
    max_total_thickness : float or None
        If set, d_S + d_F is clamped to this value during optimisation.

    Returns
    -------
    dict
        ``{"d_S": float, "d_F": float, "gamma": float,
          "Tc": float, "evaluations": int}``
    """
    base = dict(Tc0=Tc0, E_ex=E_ex, xi_S=xi_S, xi_F=xi_F,
                gamma_B=gamma_B, D_F=D_F, model=model, phase=phase,
                depairing=depairing or {})

    # Build optimisation vector
    param_names = []
    x0 = []
    bounds_lo = []
    bounds_hi = []

    if vary_d_S is not None:
        param_names.append("d_S")
        x0.append(np.clip(d_S, vary_d_S[0], vary_d_S[1]))
        bounds_lo.append(vary_d_S[0])
        bounds_hi.append(vary_d_S[1])
    if vary_d_F is not None:
        param_names.append("d_F")
        x0.append(np.clip(np.mean(vary_d_F), vary_d_F[0], vary_d_F[1]))
        bounds_lo.append(vary_d_F[0])
        bounds_hi.append(vary_d_F[1])
    if vary_gamma is not None:
        param_names.append("gamma")
        x0.append(np.clip(gamma, vary_gamma[0], vary_gamma[1]))
        bounds_lo.append(vary_gamma[0])
        bounds_hi.append(vary_gamma[1])

    if not param_names:
        raise ValueError("At least one of vary_d_S, vary_d_F, vary_gamma "
                         "must be set")

    bounds_lo = np.array(bounds_lo)
    bounds_hi = np.array(bounds_hi)
    evals = [0]

    def _unpack_vec(params):
        p = np.clip(params, bounds_lo, bounds_hi)
        cur = dict(zip(param_names, p))
        ds = cur.get("d_S", d_S)
        df = cur.get("d_F", d_S)  # fallback unused
        if max_total_thickness is not None and ds + df > max_total_thickness:
            excess = ds + df - max_total_thickness
            if "d_F" in param_names and "d_S" in param_names:
                p[param_names.index("d_F")] -= excess / 2
                p[param_names.index("d_S")] -= excess / 2
            elif "d_F" in param_names:
                p[param_names.index("d_F")] -= excess
            elif "d_S" in param_names:
                p[param_names.index("d_S")] -= excess
            p = np.clip(p, bounds_lo, bounds_hi)
        cur_d_S, cur_d_F, cur_gamma = d_S, 1.0, gamma
        for name, val in zip(param_names, p):
            if name == "d_S":
                cur_d_S = float(val)
            elif name == "d_F":
                cur_d_F = float(val)
            else:
                cur_gamma = float(val)
        return cur_d_S, cur_d_F, cur_gamma, p

    def _objective(params):
        evals[0] += 1
        cur_d_S, cur_d_F, cur_gamma, _ = _unpack_vec(np.atleast_1d(params))
        tc = _tc_scalar(cur_d_S, cur_d_F, cur_gamma, base)
        return (tc - Tc_target) ** 2

    # --- Grid pre-search for a good starting point ---
    n_grid = 20
    grid = np.linspace(bounds_lo, bounds_hi, n_grid)
    if len(param_names) == 1:
        grid = grid.reshape(-1, 1)
    else:
        mesh = np.meshgrid(*[np.linspace(lo, hi, n_grid)
                             for lo, hi in zip(bounds_lo, bounds_hi)])
        grid = np.column_stack([m.ravel() for m in mesh])
    best_cost = np.inf
    for row in grid:
        c = _objective(row)
        if c < best_cost:
            best_cost = c
            x0 = list(row)
    evals[0] = 0  # reset counter after pre-search

    # --- Optimise (Nelder-Mead from best grid point) ---
    res = _minimize(_objective, x0, method="Nelder-Mead",
                    options={"maxiter": 2000, "xatol": 1e-6, "fatol": 1e-10})

    final_d_S, final_d_F, final_gamma, _ = _unpack_vec(res.x)

    tc_final = _tc_scalar(final_d_S, final_d_F, final_gamma, base)

    return {
        "d_S": final_d_S,
        "d_F": final_d_F,
        "gamma": final_gamma,
        "Tc": tc_final,
        "evaluations": evals[0],
    }


# ---------------------------------------------------------------------------
# robust_optimize
# ---------------------------------------------------------------------------

def robust_optimize(Tc0, d_S, xi_S, xi_F, E_ex,
                    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
                    model="thin_s", phase="zero", depairing=None,
                    Tc_target=4.2,
                    vary_d_S=None, vary_d_F=None, vary_gamma=None,
                    max_total_thickness=None,
                    tolerances=None, robustness_weight=0.5):
    """
    Optimise for a target Tc while penalising sensitivity to fabrication error.

    Augmented objective::

        L = (Tc - Tc_target)^2 + lambda * sum(|dTc/dp_i| * sigma_i)

    where *sigma_i* comes from *tolerances* and *lambda* is
    *robustness_weight*.

    Parameters
    ----------
    tolerances : dict, optional
        Expected fabrication uncertainties, e.g.
        ``{"d_F": 0.5, "d_S": 1.0, "gamma": 0.02}``.
        Keys not present default to 0.
    robustness_weight : float
        Trade-off parameter lambda. Default 0.5.

    Returns
    -------
    dict
        Same as :func:`optimize_design` plus
        ``{"sensitivity": {"dTc_ddF", "dTc_ddS", "dTc_dgamma", "Tc"}}``.
    """
    tol = tolerances or {}
    base = dict(Tc0=Tc0, E_ex=E_ex, xi_S=xi_S, xi_F=xi_F,
                gamma_B=gamma_B, D_F=D_F, model=model, phase=phase,
                depairing=depairing or {})

    # Build optimisation vector (same as optimize_design)
    param_names = []
    x0 = []
    bounds_lo = []
    bounds_hi = []

    if vary_d_S is not None:
        param_names.append("d_S")
        x0.append(np.clip(d_S, vary_d_S[0], vary_d_S[1]))
        bounds_lo.append(vary_d_S[0])
        bounds_hi.append(vary_d_S[1])
    if vary_d_F is not None:
        param_names.append("d_F")
        x0.append(np.clip(np.mean(vary_d_F), vary_d_F[0], vary_d_F[1]))
        bounds_lo.append(vary_d_F[0])
        bounds_hi.append(vary_d_F[1])
    if vary_gamma is not None:
        param_names.append("gamma")
        x0.append(np.clip(gamma, vary_gamma[0], vary_gamma[1]))
        bounds_lo.append(vary_gamma[0])
        bounds_hi.append(vary_gamma[1])

    if not param_names:
        raise ValueError("At least one of vary_d_S, vary_d_F, vary_gamma "
                         "must be set")

    bounds_lo = np.array(bounds_lo)
    bounds_hi = np.array(bounds_hi)
    evals = [0]

    def _unpack(params):
        p = np.clip(np.atleast_1d(params), bounds_lo, bounds_hi)
        cur = dict(zip(param_names, p))
        ds = cur.get("d_S", d_S)
        df = cur.get("d_F", d_S)
        if max_total_thickness is not None and ds + df > max_total_thickness:
            excess = ds + df - max_total_thickness
            if "d_F" in param_names and "d_S" in param_names:
                p[param_names.index("d_F")] -= excess / 2
                p[param_names.index("d_S")] -= excess / 2
            elif "d_F" in param_names:
                p[param_names.index("d_F")] -= excess
            elif "d_S" in param_names:
                p[param_names.index("d_S")] -= excess
            p = np.clip(p, bounds_lo, bounds_hi)
        cur_d_S, cur_d_F, cur_gamma = d_S, 1.0, gamma
        for name, val in zip(param_names, p):
            if name == "d_S":
                cur_d_S = float(val)
            elif name == "d_F":
                cur_d_F = float(val)
            else:
                cur_gamma = float(val)
        return cur_d_S, cur_d_F, cur_gamma

    def _objective(params):
        evals[0] += 1
        cur_d_S, cur_d_F, cur_gamma = _unpack(params)
        tc = _tc_scalar(cur_d_S, cur_d_F, cur_gamma, base)
        cost = (tc - Tc_target) ** 2

        # sensitivity penalty
        sens = sensitivity_at(
            Tc0=Tc0, d_S=cur_d_S, d_F=cur_d_F, E_ex=E_ex,
            xi_S=xi_S, xi_F=xi_F, gamma=cur_gamma, gamma_B=gamma_B,
            D_F=D_F, model=model, phase=phase, depairing=depairing,
        )
        penalty = (abs(sens["dTc_ddF"]) * tol.get("d_F", 0)
                   + abs(sens["dTc_ddS"]) * tol.get("d_S", 0)
                   + abs(sens["dTc_dgamma"]) * tol.get("gamma", 0))
        return cost + robustness_weight * penalty

    # --- Grid pre-search for a good starting point ---
    n_grid = 20
    if len(param_names) == 1:
        grid = np.linspace(bounds_lo, bounds_hi, n_grid).reshape(-1, 1)
    else:
        mesh = np.meshgrid(*[np.linspace(lo, hi, n_grid)
                             for lo, hi in zip(bounds_lo, bounds_hi)])
        grid = np.column_stack([m.ravel() for m in mesh])
    best_cost = np.inf
    for row in grid:
        c = _objective(row)
        if c < best_cost:
            best_cost = c
            x0 = list(row)
    evals[0] = 0  # reset counter after pre-search

    # --- Optimise (Nelder-Mead from best grid point) ---
    res = _minimize(_objective, x0, method="Nelder-Mead",
                    options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-10})

    final_d_S, final_d_F, final_gamma = _unpack(res.x)
    tc_final = _tc_scalar(final_d_S, final_d_F, final_gamma, base)

    sens_final = sensitivity_at(
        Tc0=Tc0, d_S=final_d_S, d_F=final_d_F, E_ex=E_ex,
        xi_S=xi_S, xi_F=xi_F, gamma=final_gamma, gamma_B=gamma_B,
        D_F=D_F, model=model, phase=phase, depairing=depairing,
    )

    return {
        "d_S": final_d_S,
        "d_F": final_d_F,
        "gamma": final_gamma,
        "Tc": tc_final,
        "evaluations": evals[0],
        "sensitivity": sens_final,
    }
