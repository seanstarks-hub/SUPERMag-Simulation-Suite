# Changelog

All notable changes to SUPERMag are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.1.0] — 2026-04-10

### Added

- **Proximity-effect solver** — `critical_temperature()` for Tc(d_F) via digamma self-consistency (thin-S and Fominov models, zero/pi junction phases, optional depairing channels).
- **Pair amplitude** — `pair_amplitude()` computes oscillating F(x) in the ferromagnet layer.
- **Material database** — 8 built-in materials: Nb, Pb, Al (superconductors); Fe, Co, Ni, Py, Cu₀.₄₃Ni₀.₅₇ (ferromagnets). Runtime registration supported.
- **Plotting utilities** — `plot_tc_vs_df()` and `plot_pair_amplitude()` for publication-quality figures.
- **Theme system** — 4 matplotlib presets: publication, presentation, draft, dark. Context manager support.
- **Validation suite** — Automated comparison against Buzdin (1982) and Ryazanov (2003) reference data.
- **C++ backend** — AVX2 SIMD-accelerated tridiagonal solver and proximity kernels with pure-Python fallback.
- **Stub interfaces** — Usadel, Eilenberger, BdG, Ginzburg-Landau, Josephson, and spin-triplet solver signatures defined.
- **CI/CD** — GitHub Actions workflows for C++ tests, Python tests, wheel building, and weekly validation.
