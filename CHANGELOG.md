# Changelog

All notable changes to SUPERMag are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.3.0] — 2026-04-14

### Added

- **`device.ml`** — Layer-stack parser and resolver. `parse_stack "Nb:50/Fe:10"` extracts typed geometry; `resolve` maps bilayer/trilayer/graded/domains stacks to `proximity_params`, auto-detecting bilayer (S/F) and trilayer (S/N/F). Graded and domain geometries require explicit `geometry_hint` + `geom_config`.
- **`design.ml`** — Combinatorial bilayer design explorer. `enumerate_bilayers` runs SC×FM Cartesian product (3×6 = 18 combos) with parallel `Domain.spawn`, returning `design_result` records with material names and `tc_result`. `filter` applies `tc_min`/`tc_max`/`max_d_total` constraints. `to_csv`/`to_json` output formatters.
- **Interface catalogue** (`INTERFACES`) — 18 (SC, FM) → (γ, γ_B, source) entries for all built-in superconductor/ferromagnet pairs. 3 validated against published data (Buzdin 1982, Radovic 1991, Fominov 2002); 15 estimated from validated pairs.
- **Fabrication tiers** (`FABRICATION_TIERS`) — Additive γ_B corrections for interface quality: `"clean"` (+0.00), `"sputtered"` (+0.05), `"oxidized"` (+0.30).
- **`get_interface()`** — Retrieves per-pair (γ, γ_B) with fabrication-tier adjustment. Replaces fixed γ=0.3 across all combinations.
- **`rho` field** — Bulk resistivity (μΩ·cm) added to all 9 built-in materials (3 SC + 6 FM).
- **CLI `--stack` flag** — Device stack notation (`--stack "Nb:50/Fe:10"`) on `sweep_driver.ml`. Parses via `Device.parse_stack`, resolves via `Device.resolve`. Takes priority over `--sc`/`--fm`.
- **CLI `--explore` flag** — Combinatorial exploration mode on `sweep_driver.ml`. Enumerates all SC×FM bilayers over `--range` as d_F array, applies constraint filters (`--tc-min`, `--tc-max`, `--max-d-total`), outputs CSV or JSON.
- **CLI `--domain-width`, `--domain-wall`** — Domain geometry configuration flags forwarded to `Device.resolve`.
- **`test_device.ml`** — 14 Alcotest cases for `parse_stack` and `resolve` (bilayer, trilayer, graded, domains, error paths).
- **`test_design.ml`** — 5 Alcotest cases for `enumerate_bilayers` and `filter` (full count, subset, tc_min, empty, no constraints).
- **`optimize.ml`** — Multi-parameter design optimizer. Nelder-Mead in N-dimensional free-parameter subspace (d_S, d_F, γ, γ_B, E_ex) with fabrication constraint projection. Objectives: `Target_tc`, `Minimize_tc`, `Maximize_ic` (via Josephson CPR), `Multi` (weighted combination). `sensitivity_at` computes ∂Tc/∂p via central finite differences. `robust_optimize` adds sensitivity penalty for design-for-manufacturing.
- **CLI `--optimize`** — Optimization mode on `sweep_driver.ml`. `--target-tc`, `--vary name:lo,hi` (repeatable), `--robust-tol` flags.
- **`test_optimize.ml`** — 7 Alcotest cases: Nelder-Mead on quadratic, constraint projection (clips + passthrough), sensitivity finiteness, Target_tc 1D (vs golden-section), no-free-params fallback, free-param count.
- **Tutorial 07** (`07_explorer.ipynb`) — Combinatorial SC×FM explorer: enumerates all 18 bilayer combinations with per-pair γ/γ_B from the interface catalogue, uses Fominov model for finite-thickness S layers, visualizes min-Tc heatmap and normalized suppression matrix, ranks pairs by target temperature with γ/source metadata, demonstrates constraint filtering.
- **`optimizer.py`** — Multi-parameter design optimization for S/F heterostructures. `optimize_design()` finds (d_S, d_F, γ) via Nelder-Mead to hit a target Tc. `sensitivity_at()` computes ∂Tc/∂p via central finite differences. `robust_optimize()` adds sensitivity penalty for fabrication-tolerant designs.
- **Tutorial 08** (`08_design_optimizer.ipynb`) — Design optimizer workflow: Nb/Cu₀.₄₃Ni₀.₅₇ π-junction with Tc = 5 K target. Covers Tc landscape visualization, single-objective and multi-parameter optimization, sensitivity analysis, robust optimization, and design comparison.
- **Tutorial 09** (`09_supercurrent_diode.ipynb`) — Supercurrent diode effect via asymmetric SQUID model: two S/F/S junctions with different d_F in a flux-biased loop. Computes CPR, extracts I_c⁺/I_c⁻, sweeps d_F₂, external flux, and temperature.
- **Tutorial 10** (`10_pi_junction.ipynb`) — π-junction phase diagram: maps 0–π boundary by comparing Tc(phase="zero") vs Tc(phase="pi"), multi-material comparison (CuNi/Py/Ni), sweet-spot thickness finder, interface barrier effects.
- **Tutorial 11** (`11_triplet_correlations.ipynb`) — Spin-triplet correlations in magnetic multilayers: collinear vs non-collinear bilayers, sin(α) angular dependence, domain-wall triplet generation, wall-width scaling, F1/N/F2 spin-valve switching.

### Changed

- **OCaml `params.ml`** — Added `domain_wall` field to `domain_params` (aligns with C `supermag_domain_params_t.domain_wall`).
- **OCaml `stubs.ml`** — Added `solver_options_struct` ctypes binding with 6 fields matching `supermag_solver_options_t`.
- **OCaml `stubs.mli`** — Added `?tc0:float` to triplet signature; exposed `solver_options` type.
- **OCaml `solvers.ml`** — Added optional `?opts` parameter to usadel, eilenberger, josephson, and GL wrappers.
- **OCaml `sweep_driver.ml`** — Refactored to 25 parameters with tri-mode (explore/optimize/sweep). `--param` changed from required to optional (default `d_F`). Added `build_params_from_stack` helper. Added `--optimize`, `--target-tc`, `--vary`, `--robust-tol` flags.
- **`architecture.md`** — Updated §1 layer diagram and §5 file locations with `device.ml`, `design.ml`, `test_device.ml`, `test_design.ml`.

## [0.2.1] — 2026-04-14

### Added

- **`supermag_solver_options_t`** — Composable options struct (`solver_options.h`) with `matsubara_max`, `omega_cut_factor`, `max_steps`, `max_iter`, `conv_tol`, `root_grid_points`. Pass NULL for defaults.
- **Domain wall interpolation** — `kernel_domains` now supports finite-width domain walls via `domain_wall` field in `supermag_domain_params_t`. WALL_SLICES=10 with linear exchange-energy interpolation.
- **`test_plotting.py`** — 10 new tests covering `plot_pair_amplitude()` and `plot_tc_vs_df()` (save_path, axes reuse, Tc0 reference line).

### Changed

- **Triplet Python API** — Added `T`, `Tc0`, `E_ex_per_layer`, `D_per_layer`, `mode` parameters. Input validation: `n_layers >= 2`, array length checks.
- **Josephson Python API** — Added `gamma_B` parameter for interface barrier damping.
- **Ginzburg-Landau Python API** — Added `H_applied`, `mode` ("scalar"/"gauge"), `seed` parameters. Input validation: `nx > 0`, `ny > 0`, `dx > 0`, `beta > 0`.
- **BdG Python API** — Added input validation (`n_sites > 0`, `t_hop > 0`).
- **C API signatures** — Usadel, Eilenberger, Josephson, GL now accept `const supermag_solver_options_t *opts`. Triplet signature reordered with `Tc0` added.
- **OCaml FFI** — 5 new foreign declarations + wrapper functions for updated solver signatures.
- **README.md** — Fixed `gamma=0.15` → `0.10`, corrected "usage" section heading.
- **Default interface transparency** — `gamma` default changed from 0.15 to 0.10 across notebooks and examples.

### Fixed

- **Dead code removal** — Removed unused `_tridiag_solve()` from `usadel.py`.
- **Notebook 02** (`02_parameter_sweep.ipynb`) — Fixed stale cell title.
- **Documentation** — Removed stale "not yet exposed" C++ notes from josephson.md, ginzburg_landau.md, triplet.md. Updated all API docs to match current signatures.

## [0.2.0] — 2026-04-12

### Added

- **Depairing module** (`python/supermag/depairing.py`) — Individual channel computations (`depairing_ag`, `depairing_zeeman`, `depairing_orbital_perp`, `depairing_orbital_par`, `depairing_soc`) and composite `depairing_from_physical()`. Pure-Python fallback with CODATA constants + native pybind11 dispatch.
- **Optimizer utilities** — `optimize_tc()` (golden-section search), `inverse_tc()` (Brent root-finding), `fit_tc()` (Nelder-Mead least-squares) for matching Tc(d_F) to target values or experimental data. Python fallback + native dispatch.
- **pybind11 bindings** for all 9 depairing/optimizer C functions in `_binding.cpp`. All 7 solver modules + 9 depairing/optimizer functions now have native dispatch.
- **OCaml Phase 4** — 8 new ctypes FFI bindings (`depairing_ag`..`soc`, `depairing_from_physical`, `optimize_tc`, `inverse_tc`, `fit_tc`), typed `solvers.ml` wrappers, `chain.ml` Result-monad pipeline with `Domain.spawn` parallelism, `--depairing` CLI flag in `sweep_driver.ml`, 5 new FFI tests.
- **Tutorial 05** — Depairing channels: physics, computation, and Tc suppression.
- **Tutorial 06** — Fitting Tc(d_F) data to extract interface parameters.
- **Equation registry** — EQ-20 (golden-section), EQ-21 (Brent inverse), EQ-22 (Nelder-Mead fit) added to `architecture.md`.
- **Python test suite** — `test_depairing.py` with 22 tests covering all channels, composite, optimizer, inverse, and fit.

### Changed

- **architecture.md** — Updated §1 layer diagram (all pybind11 wraps complete), §2 equation registry (+3 entries), §3 C API signatures (depairing/optimizer match actual headers), §5 file locations (depairing.py, chain.ml).
- **README.md** — Replaced misleading `pip install supermag` with build-from-source instructions; replaced stale "in progress" roadmap with "Implemented Solvers" section listing all 7 complete solvers; added forward-looking roadmap (cibuildwheel, macOS CI).
- **CONTRIBUTING.md** — Added depairing/optimizer to project layout and test table.
- **`__init__.py`** — Version bumped to 0.2.0; exports depairing/optimizer public API.

### Fixed

- **optimizer.cpp** — Removed unused `e_val` variable (compiler warning).
- **OCaml `solvers.ml`** — Fixed `bdg` function missing `?mu ()` unit arg.

## [0.1.0] — 2026-04-10

### Added

- **Proximity-effect solver** — `critical_temperature()` for Tc(d_F) via digamma self-consistency (thin-S and Fominov models, zero/pi junction phases, optional depairing channels).
- **Pair amplitude** — `pair_amplitude()` computes oscillating F(x) in the ferromagnet layer.
- **Material database** — 8 built-in materials: Nb, Pb, Al (superconductors); Fe, Co, Ni, Py, Cu₀.₄₃Ni₀.₅₇ (ferromagnets). Runtime registration supported.
- **Plotting utilities** — `plot_tc_vs_df()` and `plot_pair_amplitude()` for publication-quality figures.
- **Theme system** — 4 matplotlib presets: publication, presentation, draft, dark. Context manager support.
- **Validation suite** — Automated comparison against Buzdin (1982) and Ryazanov (2003) reference data.
- **C++ backend** — AVX2 SIMD-accelerated tridiagonal solver and proximity kernels with pure-Python fallback.
- **All solvers implemented** — Usadel, Eilenberger, BdG, Ginzburg-Landau, Josephson, and spin-triplet with C++ engines and pybind11 native dispatch.
- **CI/CD** — GitHub Actions workflows for C++ tests, Python tests, wheel building, and weekly validation.
