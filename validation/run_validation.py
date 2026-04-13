#!/usr/bin/env python3
"""
Run all validation cases and report pass/fail.

Each validation directory contains:
- params.json: Physical parameters
- expected/*.csv: Expected output data

This script runs the corresponding SUPERMag solver and compares
against expected data. Fails if max absolute deviation exceeds tolerance.
"""

import json
import os
import sys
import numpy as np

from supermag.proximity import pair_amplitude, critical_temperature
from supermag.triplet import solve as triplet_solve


def run_buzdin_1982():
    """Validate against Buzdin (1982) S/F proximity effect."""
    base = os.path.join(os.path.dirname(__file__), "buzdin_1982")
    with open(os.path.join(base, "params.json")) as f:
        params = json.load(f)

    tolerance = params.get("tolerance", 1e-3)
    results = []

    # --- Tc(d_F) validation ---
    tc_file = os.path.join(base, "expected", "tc_vs_df.csv")
    if os.path.exists(tc_file):
        expected = np.genfromtxt(tc_file, delimiter=",", skip_header=1)
        d_F_expected = expected[:, 0]
        Tc_expected = expected[:, 1]

        Tc_computed = critical_temperature(
            Tc0=params["Tc0_K"], d_S=params["d_S_nm"],
            d_F_array=d_F_expected,
            E_ex=params["E_ex_meV"],
            xi_S=params["xi_S_nm"], xi_F=params["xi_F_nm"],
            gamma=params.get("gamma", 0.3),
        )

        # Normalize by Tc0 for comparison
        max_dev = np.max(np.abs(Tc_computed - Tc_expected)) / params["Tc0_K"]
        passed = max_dev < tolerance
        results.append(("Buzdin1982/Tc_vs_dF", passed, max_dev, tolerance))

    # --- Pair amplitude validation ---
    pa_file = os.path.join(base, "expected", "pair_amplitude.csv")
    if os.path.exists(pa_file):
        expected = np.genfromtxt(pa_file, delimiter=",", skip_header=1)
        x_expected = expected[:, 0]
        F_expected = expected[:, 1]

        x_computed, F_computed = pair_amplitude(
            d_F=x_expected[-1], xi_F=params["xi_F_nm"],
            n_points=len(x_expected),
        )

        max_dev = np.max(np.abs(F_computed - F_expected))
        passed = max_dev < tolerance
        results.append(("Buzdin1982/pair_amplitude", passed, max_dev, tolerance))

    return results


def run_ryazanov_2003():
    """Validate against Ryazanov (2003) / Fominov (2002) Nb/Cu0.43Ni0.57 data."""
    base = os.path.join(os.path.dirname(__file__), "ryazanov_2003")
    with open(os.path.join(base, "params.json")) as f:
        params = json.load(f)

    tolerance = params.get("tolerance", 0.05)
    results = []

    # --- Tc(d_F) validation using Fominov model ---
    tc_file = os.path.join(base, "expected", "tc_vs_df.csv")
    if os.path.exists(tc_file):
        expected = np.genfromtxt(tc_file, delimiter=",", skip_header=1)
        d_F_expected = expected[:, 0]
        Tc_expected = expected[:, 1]

        # Filter out d_F = 0 since solver needs d_F > 0
        mask = d_F_expected > 0
        d_F_nz = d_F_expected[mask]
        Tc_exp_nz = Tc_expected[mask]

        Tc_computed = critical_temperature(
            Tc0=params["Tc0_K"], d_S=params["d_S_nm"],
            d_F_array=d_F_nz,
            E_ex=params["E_ex_meV"],
            xi_S=params["xi_S_nm"], xi_F=params["xi_F_nm"],
            gamma=params.get("gamma", 0.15),
            gamma_B=params.get("gamma_B", 0.3),
            D_F=params.get("D_F_m2s", 4.0e-4),
            model=params.get("model", "fominov"),
            phase=params.get("phase", "zero"),
        )

        # Normalize by Tc0 for comparison
        max_dev = np.max(np.abs(Tc_computed - Tc_exp_nz)) / params["Tc0_K"]
        passed = max_dev < tolerance
        results.append(("Ryazanov2003/Tc_vs_dF", passed, max_dev, tolerance))

    return results


def run_radovic_1991():
    """Validate against Radović (1991) Nb/Ni superlattice Tc(d_F)."""
    base = os.path.join(os.path.dirname(__file__), "radovic_1991")
    with open(os.path.join(base, "params.json")) as f:
        params = json.load(f)

    tolerance = params.get("tolerance", 1e-3)
    results = []

    tc_file = os.path.join(base, "expected", "tc_vs_df.csv")
    if os.path.exists(tc_file):
        expected = np.genfromtxt(tc_file, delimiter=",", skip_header=1)
        d_F_expected = expected[:, 0]
        Tc_expected = expected[:, 1]

        Tc_computed = critical_temperature(
            Tc0=params["Tc0_K"], d_S=params["d_S_nm"],
            d_F_array=d_F_expected,
            E_ex=params["E_ex_meV"],
            xi_S=params["xi_S_nm"], xi_F=params["xi_F_nm"],
            gamma=params.get("gamma", 0.3),
            model=params.get("model", "thin_s"),
        )

        max_dev = np.max(np.abs(Tc_computed - Tc_expected)) / params["Tc0_K"]
        passed = max_dev < tolerance
        results.append(("Radovic1991/Tc_vs_dF", passed, max_dev, tolerance))

    return results


def run_bergeret_2005():
    """Validate against Bergeret (2005) S/F1/N/F2/S triplet correlations."""
    base = os.path.join(os.path.dirname(__file__), "bergeret_2005")
    with open(os.path.join(base, "params.json")) as f:
        params = json.load(f)

    tolerance = params.get("tolerance", 1e-6)
    results = []

    ta_file = os.path.join(base, "expected", "triplet_amplitude.csv")
    if os.path.exists(ta_file):
        expected = np.genfromtxt(ta_file, delimiter=",", skip_header=1)
        x_expected = expected[:, 0]
        f_expected = expected[:, 1]

        x_computed, f_computed = triplet_solve(
            n_layers=params["n_layers"],
            thicknesses=params["thicknesses_nm"],
            magnetization_angles=params["magnetization_angles_rad"],
            n_grid=params.get("n_grid", 200),
            xi_F=params.get("xi_F_nm", 1.0),
            xi_N=params.get("xi_N_nm", 10.0),
        )

        max_dev = np.max(np.abs(f_computed - f_expected))
        passed = max_dev < tolerance
        results.append(("Bergeret2005/triplet_amplitude", passed, max_dev, tolerance))

    return results


def main():
    all_results = []

    print("=" * 60)
    print("SUPERMag Validation Suite")
    print("=" * 60)

    # Run all validation cases
    try:
        all_results.extend(run_buzdin_1982())
    except Exception as e:
        all_results.append(("Buzdin1982", False, str(e), "N/A"))

    try:
        all_results.extend(run_ryazanov_2003())
    except Exception as e:
        all_results.append(("Ryazanov2003", False, str(e), "N/A"))

    try:
        all_results.extend(run_radovic_1991())
    except Exception as e:
        all_results.append(("Radovic1991", False, str(e), "N/A"))

    try:
        all_results.extend(run_bergeret_2005())
    except Exception as e:
        all_results.append(("Bergeret2005", False, str(e), "N/A"))

    # Report
    print()
    all_pass = True
    for name, passed, deviation, tol in all_results:
        status = "PASS" if passed else "FAIL"
        if passed:
            print(f"  [{status}] {name}: max_dev={deviation:.2e} < tol={tol:.2e}")
        else:
            print(f"  [{status}] {name}: max_dev={deviation} >= tol={tol}")
            all_pass = False

    print()
    if all_pass:
        print("All validation cases PASSED.")
        return 0
    else:
        print("Some validation cases FAILED.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
