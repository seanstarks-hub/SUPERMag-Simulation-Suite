"""
Material parameter database for superconductor/ferromagnet heterostructures.

All values are physically reasonable order-of-magnitude estimates.

Units:
    Tc      — Critical temperature (K)
    xi_S    — Superconductor coherence length (nm)
    lambda_L — London penetration depth (nm)
    Delta_0 — Superconducting gap at T=0 (meV)
    rho     — Normal-state resistivity (μΩ·cm)
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
        "rho": 15.0,         # μΩ·cm — dirty Nb thin film
    },
    "Pb": {
        "type": "superconductor",
        "Tc": 7.2,           # K
        "xi_S": 83.0,        # nm
        "lambda_L": 37.0,    # nm
        "Delta_0": 1.35,     # meV
        "rho": 22.0,         # μΩ·cm
    },
    "Al": {
        "type": "superconductor",
        "Tc": 1.2,           # K
        "xi_S": 1600.0,      # nm — very long in clean Al
        "lambda_L": 16.0,    # nm
        "Delta_0": 0.18,     # meV
        "rho": 2.7,          # μΩ·cm — clean Al
    },
    # --- Ferromagnets ---
    "Fe": {
        "type": "ferromagnet",
        "E_ex": 256.0,       # meV
        "xi_F": 0.7,         # nm
        "D_F": 2.5e-4,       # m^2/s
        "rho": 10.0,         # μΩ·cm
    },
    "Co": {
        "type": "ferromagnet",
        "E_ex": 309.0,       # meV
        "xi_F": 0.5,         # nm
        "D_F": 1.8e-4,       # m^2/s
        "rho": 6.3,          # μΩ·cm
    },
    "Ni": {
        "type": "ferromagnet",
        "E_ex": 75.0,        # meV
        "xi_F": 2.3,         # nm
        "D_F": 5.0e-4,       # m^2/s
        "rho": 6.9,          # μΩ·cm
    },
    "Py": {
        "type": "ferromagnet",
        "E_ex": 20.0,        # meV — Permalloy (Ni80Fe20), weak ferromagnet
        "xi_F": 5.0,         # nm
        "D_F": 3.0e-4,       # m^2/s
        "rho": 40.0,         # μΩ·cm — high-resistivity alloy
    },
    "CuNi": {
        "type": "ferromagnet",
        "E_ex": 5.0,         # meV — dilute ferromagnetic alloy
        "xi_F": 10.0,        # nm
        "D_F": 4.0e-4,       # m^2/s
        "rho": 35.0,         # μΩ·cm
    },
    "Cu0.43Ni0.57": {
        "type": "ferromagnet",
        "E_ex": 11.2,        # meV — h~130K from Fominov fit (kB*130K ≈ 11.2 meV)
        "xi_F": 4.2,         # nm — sqrt(hbar*D_F/(2*E_ex))
        "D_F": 4.0e-4,       # m^2/s
        "rho": 50.0,         # μΩ·cm
    },
}


# ── S/F Interface Catalogue ──────────────────────────────────────────────────
#
# Phenomenological (gamma, gamma_B) values for each SC×FM pair.
# gamma is the effective pair-breaking coupling; gamma_B is the interface
# barrier parameter.  These are NOT computable from bulk resistivities alone —
# they are fitted to experimental Tc(d_F) data or estimated from validated
# pairs.  "validated" entries are anchored to the validation suite;
# "estimated" entries are extrapolated.
#
# gamma_B encodes both intrinsic interface properties (alloy interdiffusion)
# and fabrication quality.  Use FABRICATION_TIERS to adjust for deposition
# conditions.

INTERFACES = {
    # ── Validated (anchored to validation suite) ────
    ("Nb", "Fe"):           {"gamma": 0.30, "gamma_B": 0.00,
                             "source": "validated", "ref": "Buzdin, JETP Lett. 35, 178 (1982)"},
    ("Nb", "Ni"):           {"gamma": 0.30, "gamma_B": 0.00,
                             "source": "validated", "ref": "Radovic, PRB 44, 759 (1991)"},
    ("Nb", "Cu0.43Ni0.57"): {"gamma": 0.15, "gamma_B": 0.30,
                             "source": "validated", "ref": "Fominov, PRB 66, 014507 (2002)"},
    # ── Estimated Nb pairs ──────────────────────────
    ("Nb", "Co"):           {"gamma": 0.28, "gamma_B": 0.00,
                             "source": "estimated", "ref": "scaled from Nb/Fe"},
    ("Nb", "Py"):           {"gamma": 0.12, "gamma_B": 0.05,
                             "source": "estimated", "ref": "Rusanov (2004) range"},
    ("Nb", "CuNi"):         {"gamma": 0.10, "gamma_B": 0.25,
                             "source": "estimated", "ref": "scaled from Nb/Cu0.43Ni0.57"},
    # ── Estimated Pb pairs (stronger coupling than Nb) ──
    ("Pb", "Fe"):           {"gamma": 0.45, "gamma_B": 0.00,
                             "source": "estimated", "ref": "scaled from Nb/Fe, Pb mismatch"},
    ("Pb", "Co"):           {"gamma": 0.42, "gamma_B": 0.00,
                             "source": "estimated", "ref": "scaled from Pb/Fe"},
    ("Pb", "Ni"):           {"gamma": 0.40, "gamma_B": 0.00,
                             "source": "estimated", "ref": "scaled from Pb/Fe"},
    ("Pb", "Py"):           {"gamma": 0.18, "gamma_B": 0.05,
                             "source": "estimated", "ref": "scaled from Nb/Py, Pb mismatch"},
    ("Pb", "CuNi"):         {"gamma": 0.15, "gamma_B": 0.25,
                             "source": "estimated", "ref": "scaled from Nb/CuNi"},
    ("Pb", "Cu0.43Ni0.57"): {"gamma": 0.22, "gamma_B": 0.30,
                             "source": "estimated", "ref": "scaled from Nb/Cu0.43Ni0.57"},
    # ── Estimated Al pairs (weak-coupling, low Tc0) ─
    ("Al", "Fe"):           {"gamma": 0.08, "gamma_B": 0.00,
                             "source": "estimated", "ref": "Al weak coupling"},
    ("Al", "Co"):           {"gamma": 0.07, "gamma_B": 0.00,
                             "source": "estimated", "ref": "Al weak coupling"},
    ("Al", "Ni"):           {"gamma": 0.08, "gamma_B": 0.00,
                             "source": "estimated", "ref": "Al weak coupling"},
    ("Al", "Py"):           {"gamma": 0.03, "gamma_B": 0.05,
                             "source": "estimated", "ref": "Al weak coupling"},
    ("Al", "CuNi"):         {"gamma": 0.02, "gamma_B": 0.25,
                             "source": "estimated", "ref": "Al weak coupling"},
    ("Al", "Cu0.43Ni0.57"): {"gamma": 0.03, "gamma_B": 0.30,
                             "source": "estimated", "ref": "Al weak coupling"},
}


# ── Fabrication Quality Tiers ────────────────────────────────────────────────
#
# gamma_B accounts for both intrinsic (alloy interdiffusion) and extrinsic
# (oxide, roughness) interface barriers.  The INTERFACES catalogue gives base
# values for clean interfaces.  These tiers add a correction to gamma_B that
# reflects typical laboratory deposition conditions.

FABRICATION_TIERS = {
    "clean":     {"gamma_B_add": 0.00,
                  "desc": "In-situ MBE / e-beam, no vacuum break, base < 1e-9 Torr"},
    "sputtered": {"gamma_B_add": 0.05,
                  "desc": "DC/RF magnetron sputtering, base < 1e-7 Torr"},
    "oxidized":  {"gamma_B_add": 0.30,
                  "desc": "Ex-situ transfer with native oxide at S/F interface"},
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


def get_interface(sc_name, fm_name, tier="sputtered"):
    """Look up interface parameters for an S/F pair.

    Parameters
    ----------
    sc_name : str
        Superconductor name (e.g., "Nb").
    fm_name : str
        Ferromagnet name (e.g., "Fe").
    tier : str, optional
        Fabrication quality: ``"clean"``, ``"sputtered"``, or ``"oxidized"``.
        Adjusts gamma_B.  Default: ``"sputtered"``.

    Returns
    -------
    dict
        ``{"gamma": float, "gamma_B": float, "source": str, "ref": str, "tier": str}``

    Raises
    ------
    KeyError
        If the (sc, fm) pair is not in INTERFACES.
    ValueError
        If tier is not recognized.
    """
    if tier not in FABRICATION_TIERS:
        raise ValueError(
            f"Unknown fabrication tier '{tier}'. "
            f"Available: {', '.join(sorted(FABRICATION_TIERS))}")
    key = (sc_name, fm_name)
    if key not in INTERFACES:
        available = [f"{s}/{f}" for s, f in sorted(INTERFACES)]
        raise KeyError(
            f"Interface '{sc_name}/{fm_name}' not in catalogue. "
            f"Available: {', '.join(available)}")
    entry = dict(INTERFACES[key])
    entry["gamma_B"] += FABRICATION_TIERS[tier]["gamma_B_add"]
    entry["tier"] = tier
    return entry
