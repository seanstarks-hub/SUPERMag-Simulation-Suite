"""
Material parameter database for superconductor/ferromagnet heterostructures.

All values are physically reasonable order-of-magnitude estimates.

Units:
    Tc      — Critical temperature (K)
    xi_S    — Superconductor coherence length (nm)
    lambda_L — London penetration depth (nm)
    Delta_0 — Superconducting gap at T=0 (meV)
    E_ex    — Exchange energy (meV)
    xi_F    — Ferromagnet coherence length (nm)
    D_F     — Diffusion coefficient in ferromagnet (m^2/s)
"""

MATERIALS = {
    # --- Superconductors ---
    "Nb": {
        "type": "superconductor",
        "Tc": 9.2,           # K
        "xi_S": 38.0,        # nm — dirty-limit coherence length
        "lambda_L": 39.0,    # nm — London penetration depth
        "Delta_0": 1.55,     # meV — BCS gap at T=0
    },
    "Pb": {
        "type": "superconductor",
        "Tc": 7.2,           # K
        "xi_S": 83.0,        # nm
        "lambda_L": 37.0,    # nm
        "Delta_0": 1.35,     # meV
    },
    "Al": {
        "type": "superconductor",
        "Tc": 1.2,           # K
        "xi_S": 1600.0,      # nm — very long in clean Al
        "lambda_L": 16.0,    # nm
        "Delta_0": 0.18,     # meV
    },
    # --- Ferromagnets ---
    "Fe": {
        "type": "ferromagnet",
        "E_ex": 256.0,       # meV
        "xi_F": 0.7,         # nm
        "D_F": 2.5e-4,       # m^2/s
    },
    "Co": {
        "type": "ferromagnet",
        "E_ex": 309.0,       # meV
        "xi_F": 0.5,         # nm
        "D_F": 1.8e-4,       # m^2/s
    },
    "Ni": {
        "type": "ferromagnet",
        "E_ex": 75.0,        # meV
        "xi_F": 2.3,         # nm
        "D_F": 5.0e-4,       # m^2/s
    },
    "Py": {
        "type": "ferromagnet",
        "E_ex": 20.0,        # meV — Permalloy (Ni80Fe20), weak ferromagnet
        "xi_F": 5.0,         # nm
        "D_F": 3.0e-4,       # m^2/s
    },
    "CuNi": {
        "type": "ferromagnet",
        "E_ex": 5.0,         # meV — dilute ferromagnetic alloy
        "xi_F": 10.0,        # nm
        "D_F": 4.0e-4,       # m^2/s
    },
}


def get_material(name):
    """
    Retrieve material parameters by name.

    Parameters
    ----------
    name : str
        Material name (e.g., "Nb", "Fe", "Py").

    Returns
    -------
    dict
        Dictionary of material parameters with units documented in module docstring.

    Raises
    ------
    KeyError
        If material name is not in the database.

    Examples
    --------
    >>> from supermag.materials import get_material
    >>> nb = get_material("Nb")
    >>> print(f"Nb Tc = {nb['Tc']} K")
    Nb Tc = 9.2 K
    """
    if name not in MATERIALS:
        available = ", ".join(sorted(MATERIALS.keys()))
        raise KeyError(f"Material '{name}' not found. Available: {available}")
    return dict(MATERIALS[name])  # return a copy


def list_materials():
    """
    List all available material names, grouped by type.

    Returns
    -------
    dict
        {"superconductor": [...], "ferromagnet": [...]}
    """
    result = {"superconductor": [], "ferromagnet": []}
    for name, params in sorted(MATERIALS.items()):
        mat_type = params.get("type", "unknown")
        if mat_type in result:
            result[mat_type].append(name)
    return result
