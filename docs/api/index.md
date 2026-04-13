# API Reference

Complete reference for the `supermag` Python package (v0.1.0).

## Solvers

| Module | Function | Description | Equation |
|--------|----------|-------------|----------|
| [proximity](proximity.md) | `critical_temperature()` | $T_c(d_F)$ for S/F bilayer | EQ-4, EQ-5 |
| [proximity](proximity.md) | `pair_amplitude()` | Pair amplitude $F(x)$ in F layer | EQ-6 |
| [usadel](usadel.md) | `solve()` | Diffusive-limit order parameter $\Delta(x)$ | — |
| [eilenberger](eilenberger.md) | `solve()` | Clean-limit anomalous Green's function $f(x)$ | — |
| [bdg](bdg.md) | `solve()` | Tight-binding BdG eigenvalue spectrum | EQ-10 |
| [ginzburg_landau](ginzburg_landau.md) | `minimize()` | GL free energy relaxation on 2D grid | EQ-11 |
| [josephson](josephson.md) | `current_phase_relation()` | Josephson CPR $I(\varphi)$ via Matsubara sum | EQ-9 |
| [triplet](triplet.md) | `solve()` | Equal-spin triplet correlations in multilayers | EQ-12 |

## Sweep Engines

| Module | Function | Description |
|--------|----------|-------------|
| [sweeps](sweeps.md) | `tc_parameter_sweep()` | 1D parameter sweep of $T_c$ |
| [sweeps](sweeps.md) | `tc_phase_diagram()` | 2D phase diagram of $T_c$ |

## Utilities

| Module | Function | Description |
|--------|----------|-------------|
| [materials](materials.md) | `get_material()` | Look up material parameters by name |
| [materials](materials.md) | `list_materials()` | List available materials by type |
| [plotting](plotting.md) | `plot_pair_amplitude()` | Plot $F(x)$ in the F layer |
| [plotting](plotting.md) | `plot_tc_vs_df()` | Plot $T_c$ vs $d_F$ |
| [themes](themes.md) | `apply_theme()` | Apply matplotlib theme globally |
| [themes](themes.md) | `theme_context()` | Temporary theme via context manager |
| [themes](themes.md) | `list_themes()` | List available theme names |
| [themes](themes.md) | `get_theme()` | Get theme rcParams dict |
| [themes](themes.md) | `register_theme()` | Register a custom theme |
