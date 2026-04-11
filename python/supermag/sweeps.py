"""
Generic Tc parameter sweep utilities.

Provides functions for sweeping any single parameter or building 2D
phase diagrams from the proximity-effect critical temperature solver.
"""

import numpy as np
import scipy.constants as const

from supermag.proximity import critical_temperature


def _estimate_xi_F(D_F, E_ex_meV):
    """Compute xi_F (nm) from D_F (m^2/s) and E_ex (meV)."""
    E_ex_J = E_ex_meV * const.eV * 1e-3
    return np.sqrt(const.hbar * D_F / E_ex_J) * 1e9


# Parameters accepted by critical_temperature() and their roles
_CT_SCALAR_PARAMS = {
    "Tc0", "d_S", "E_ex", "xi_S", "xi_F",
    "gamma", "gamma_B", "D_F", "model", "phase",
}


def tc_parameter_sweep(sweep_var, sweep_values, d_F_array=None, **fixed_params):
    """
    Sweep a single parameter and compute Tc for each value.

    Parameters
    ----------
    sweep_var : str
        Name of the parameter to sweep.  Supported:
        ``"d_F"``, ``"d_S"``, ``"gamma"``, ``"gamma_B"``,
        ``"E_ex"``, ``"xi_F"``, ``"Tc0"``.
    sweep_values : array_like
        Values of *sweep_var* to iterate over.
    d_F_array : array_like, optional
        Ferromagnet thickness array (nm).  Required when *sweep_var*
        is **not** ``"d_F"``.  When *sweep_var* ``== "d_F"`` this is
        ignored (``sweep_values`` is used instead).
    **fixed_params
        All remaining keyword arguments are forwarded to
        :func:`supermag.critical_temperature`.  At minimum you need
        ``Tc0``, ``d_S``, ``E_ex``, ``xi_S``, ``xi_F`` (unless they
        are the sweep variable).

    Returns
    -------
    result : dict
        ``{"sweep_var": sweep_var,
          "sweep_values": np.ndarray,
          "d_F": np.ndarray,
          "Tc": np.ndarray}``

        * When *sweep_var* is ``"d_F"``, ``Tc`` has shape
          ``(len(sweep_values),)`` — one value per thickness.
        * Otherwise ``Tc`` has shape
          ``(len(sweep_values), len(d_F_array))``.

    Notes
    -----
    When ``sweep_var="E_ex"`` and ``"xi_F"`` is **not** in
    *fixed_params*, ``xi_F`` is automatically recomputed at each step
    via ``xi_F = sqrt(hbar * D_F / E_ex)``.  Supply ``D_F`` (default
    2.5e-4 m²/s) for this to give physical values.

    Examples
    --------
    >>> import numpy as np
    >>> from supermag.sweeps import tc_parameter_sweep
    >>> result = tc_parameter_sweep(
    ...     "gamma", np.linspace(0.05, 0.5, 10),
    ...     d_F_array=np.linspace(0.5, 20, 50),
    ...     Tc0=9.2, d_S=50.0, E_ex=256.0, xi_S=38.0, xi_F=0.7)
    >>> result["Tc"].shape
    (10, 50)
    """
    sweep_values = np.asarray(sweep_values, dtype=np.float64)

    if sweep_var == "d_F":
        # Each sweep value is a single thickness → 1-D output
        Tc_list = []
        for val in sweep_values:
            kw = dict(fixed_params, d_F_array=np.array([val]))
            Tc_list.append(critical_temperature(**kw)[0])
        return {
            "sweep_var": sweep_var,
            "sweep_values": sweep_values,
            "d_F": sweep_values,
            "Tc": np.array(Tc_list),
        }

    if d_F_array is None:
        raise ValueError(
            "d_F_array is required when sweep_var is not 'd_F'")
    d_F_array = np.asarray(d_F_array, dtype=np.float64)

    # Detect whether xi_F should be auto-recomputed from E_ex
    auto_xi_F = (sweep_var == "E_ex" and "xi_F" not in fixed_params)

    Tc_grid = np.empty((len(sweep_values), len(d_F_array)))

    for i, val in enumerate(sweep_values):
        kw = dict(fixed_params)
        kw[sweep_var] = val
        kw["d_F_array"] = d_F_array
        if auto_xi_F:
            D_F = kw.pop("D_F", 2.5e-4)
            kw["xi_F"] = _estimate_xi_F(D_F, val)
            kw["D_F"] = D_F
        Tc_grid[i] = critical_temperature(**kw)

    return {
        "sweep_var": sweep_var,
        "sweep_values": sweep_values,
        "d_F": d_F_array,
        "Tc": Tc_grid,
    }


def tc_phase_diagram(var1, vals1, var2, vals2, d_F_value=None,
                     **fixed_params):
    """
    Compute a 2-D grid of Tc values over two swept parameters.

    One axis is always ``d_F`` unless *d_F_value* is supplied (in which
    case both axes are non-d_F parameters evaluated at a single d_F).

    Parameters
    ----------
    var1, var2 : str
        Parameter names for the two axes.
    vals1, vals2 : array_like
        Values along each axis.
    d_F_value : float, optional
        If both *var1* and *var2* differ from ``"d_F"``, supply a
        single thickness in nm.
    **fixed_params
        Forwarded to :func:`supermag.critical_temperature`.

    Returns
    -------
    result : dict
        ``{var1: vals1, var2: vals2, "Tc": np.ndarray}``
        where ``Tc`` has shape ``(len(vals1), len(vals2))``.
    """
    vals1 = np.asarray(vals1, dtype=np.float64)
    vals2 = np.asarray(vals2, dtype=np.float64)

    has_df = (var1 == "d_F" or var2 == "d_F")

    if not has_df and d_F_value is None:
        raise ValueError(
            "d_F_value is required when neither axis is 'd_F'")

    Tc_grid = np.empty((len(vals1), len(vals2)))

    for i, v1 in enumerate(vals1):
        for j, v2 in enumerate(vals2):
            kw = dict(fixed_params)
            kw[var1] = v1
            kw[var2] = v2

            # Figure out d_F_array for the solver call
            if var1 == "d_F":
                kw["d_F_array"] = np.array([v1])
            elif var2 == "d_F":
                kw["d_F_array"] = np.array([v2])
            else:
                kw["d_F_array"] = np.array([d_F_value])

            # Auto-compute xi_F when sweeping E_ex
            if ("E_ex" in (var1, var2) and "xi_F" not in fixed_params
                    and var1 != "xi_F" and var2 != "xi_F"):
                D_F = kw.pop("D_F", 2.5e-4)
                kw["xi_F"] = _estimate_xi_F(D_F, kw["E_ex"])
                kw["D_F"] = D_F

            # Remove the var names from kw if they collide with
            # the d_F_array keyword
            kw.pop("d_F", None)

            Tc_grid[i, j] = critical_temperature(**kw)[0]

    return {
        var1: vals1,
        var2: vals2,
        "Tc": Tc_grid,
    }
