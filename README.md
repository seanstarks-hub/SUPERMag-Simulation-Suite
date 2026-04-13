# SUPERMag

**Computational toolkit for superconductor/ferromagnet heterostructure research**

*SUPERMag Lab, Texas State University*

![License](https://img.shields.io/badge/license-BSD--3--Clause-blue)
![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)
![Install](https://img.shields.io/badge/pip%20install-supermag-brightgreen)

---

## Installation

### Pure-Python (all platforms)

All solvers include pure-Python fallbacks. No compiler needed:

```bash
pip install .
```

```python
import supermag
```

This is the recommended path for most users. It works on Windows, Linux, and macOS (including Apple Silicon) with Python ≥ 3.9.

### With native C++ acceleration

Native C++ dispatch gives 10–100× speedups on numerically intensive solvers.
Building the extension requires a C++17 compiler and CMake ≥ 3.15. The build is
handled automatically by `pip install` — no separate `make build` step is needed.

#### Windows (first-class, CI-tested)

1. Install [Visual Studio 2022 Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022) with the **"Desktop development with C++"** workload.
2. Open a **Developer Command Prompt** (or load `vcvars64.bat`).
3. Run:
   ```bash
   pip install .
   ```

CMake is bundled with Visual Studio. If you need a standalone copy: `pip install cmake`.

#### Linux x86_64 (supported, CI-tested)

```bash
# Ubuntu / Debian
sudo apt install g++ cmake

# Fedora
sudo dnf install gcc-c++ cmake

pip install .
```

AVX2 is used automatically when the CPU supports it; the build succeeds without
it on older hardware.

### Platform support matrix

| | Windows x86_64 | Linux x86_64 | macOS |
|---|---|---|---|
| **Pure-Python** | ✅ | ✅ | ✅ |
| **Native C++** | ✅ CI-tested | ✅ CI-tested | ❌ Not yet supported |
| **Pre-built wheels** | ✅ | ✅ | Roadmap |

## Examples

### Critical temperature vs. ferromagnet thickness

Compute $T_c(d_F)$ for an Nb/Fe bilayer using the proximity-effect solver:

```python
import numpy as np
import supermag
from supermag.plotting import plot_tc_vs_df

nb = supermag.get_material("Nb")   # Tc=9.2 K, xi_S=38 nm
fe = supermag.get_material("Fe")   # E_ex=256 meV, xi_F=0.7 nm

d_F = np.linspace(0.5, 20.0, 100)
Tc = supermag.critical_temperature(
    Tc0=nb["Tc"], d_S=50.0,
    d_F_array=d_F, E_ex=fe["E_ex"],
    xi_S=nb["xi_S"], xi_F=fe["xi_F"],
    gamma=0.15,
)

plot_tc_vs_df(d_F, Tc, Tc0=nb["Tc"], save_path="tc_vs_df.png")
```

### Pair amplitude in the ferromagnet layer

Visualize the oscillating, decaying Cooper-pair amplitude $F(x)$ inside the F layer:

```python
from supermag import pair_amplitude
from supermag.plotting import plot_pair_amplitude

x, F = pair_amplitude(d_F=10.0, xi_F=2.3, phase="zero")
plot_pair_amplitude(x, F, save_path="pair_amplitude.png")
```

### Multi-material comparison

Compare $T_c$ suppression across ferromagnets with different exchange energies:

```python
import numpy as np
import supermag
import matplotlib.pyplot as plt

nb = supermag.get_material("Nb")
d_F = np.linspace(0.5, 20.0, 100)

fig, ax = plt.subplots()
for name in ["Fe", "Ni", "Py"]:
    fm = supermag.get_material(name)
    Tc = supermag.critical_temperature(
        Tc0=nb["Tc"], d_S=50.0, d_F_array=d_F,
        E_ex=fm["E_ex"], xi_S=nb["xi_S"], xi_F=fm["xi_F"],
    )
    ax.plot(d_F, Tc, label=name)

ax.set_xlabel("$d_F$ (nm)")
ax.set_ylabel("$T_c$ (K)")
ax.legend()
plt.show()
```

### Custom material registration

Define your own materials at runtime and use them like built-in entries:

```python
from supermag.materials import MATERIALS

MATERIALS["MgB2"] = {
    "type": "superconductor",
    "Tc": 39.0,
    "xi_S": 5.0,
    "lambda_L": 140.0,
    "Delta_0": 7.1,
}

mgb2 = supermag.get_material("MgB2")
```

### [WIP] Plotting themes

Switch between publication, presentation, draft, and dark themes:

```python
import supermag

# Apply a theme globally
supermag.apply_theme("publication")   # APS/PRB single-column format

# Or use a context manager for a single figure
with supermag.theme_context("presentation"):
    plot_tc_vs_df(d_F, Tc, Tc0=nb["Tc"])
```

Run `supermag.list_themes()` to see all available [WIP] presets.

## Available Materials

| Material | Type | Key Parameters |
|----------|------|----------------|
| Nb | Superconductor | $T_c$ = 9.2 K, $\xi_S$ = 38 nm, $\Delta_0$ = 1.55 meV |
| Pb | Superconductor | $T_c$ = 7.2 K, $\xi_S$ = 83 nm, $\Delta_0$ = 1.35 meV |
| Al | Superconductor | $T_c$ = 1.2 K, $\xi_S$ = 1600 nm, $\Delta_0$ = 0.18 meV |
| Fe | Ferromagnet | $E_{ex}$ = 256 meV, $\xi_F$ = 0.7 nm |
| Co | Ferromagnet | $E_{ex}$ = 309 meV, $\xi_F$ = 0.5 nm |
| Ni | Ferromagnet | $E_{ex}$ = 75 meV, $\xi_F$ = 2.3 nm |
| Py | Ferromagnet | $E_{ex}$ = 20 meV, $\xi_F$ = 5.0 nm (Permalloy) |
| CuNi | Ferromagnet | $E_{ex}$ = 5.0 meV, $\xi_F$ = 10.0 nm |
| Cu₀.₄₃Ni₀.₅₇ | Ferromagnet | $E_{ex}$ = 11.2 meV, $\xi_F$ = 4.2 nm |

Register custom materials at runtime — see the example above.

## Validation

Solver outputs are validated against published reference data:

- **Buzdin (1982)** — Oscillatory $T_c$ suppression in S/F bilayers (JETP Lett. 35, 178)
- **Radovic (1991)** — $T_c(d_F)$ in thin-S limit (Phys. Rev. B 44, 759)
- **Ryazanov (2003)** — Monotonic $T_c(d_F)$ in Nb/Cu₀.₄₃Ni₀.₅₇ (Fominov model)
- **Bergeret (2005)** — Long-range triplet amplitude in non-collinear multilayers (Rev. Mod. Phys. 77, 1321)

Run the validation suite with `make validate`. See `validation/` for details.

## Implemented Solvers

All solver modules have C++ implementations with pybind11 native dispatch and pure-Python fallbacks:

- **Proximity** — $T_c(d_F)$ via digamma self-consistency (thin-S and Fominov models)
- **Usadel** — diffusive-limit quasiclassical Green's functions
- **Eilenberger** — clean-limit Riccati parameterization
- **BdG** — tight-binding Bogoliubov–de Gennes Hamiltonian
- **Ginzburg-Landau** — free energy minimization near $T_c$
- **Josephson** — current-phase relations for S/F/S junctions
- **Spin-Triplet** — long-range triplet correlations in magnetic multilayers

### Roadmap

- macOS native C++ support and CI
- Pre-built wheels for macOS (Intel + Apple Silicon)
- GUI capabilities for least effort useage.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for build instructions, development setup, and contribution guidelines.

## License

BSD-3-Clause — see [LICENSE](LICENSE).
