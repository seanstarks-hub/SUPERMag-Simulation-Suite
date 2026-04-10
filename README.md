# SUPERMag

**Computational toolkit for superconductor/ferromagnet heterostructure research**

*SUPERMag Lab, Texas State University*

---

SUPERMag provides high-performance numerical solvers for modeling thin-film superconductor/ferromagnet (S/F) heterostructures, including proximity effects, critical temperature calculations, Josephson junctions, and spin-triplet superconductivity. The library features C++17 solvers (with AVX2 SIMD acceleration), OCaml-based pipeline orchestration, and a user-friendly Python interface installable via `pip`.

## Installation

```bash
pip install supermag
```

No C++ compiler, OCaml, or Make required — just install the Python wheel and import it.

## Quickstart

```python
import numpy as np
import supermag
from supermag.plotting import plot_tc_vs_df

# Material parameters
nb = supermag.get_material("Nb")   # Tc=9.2K, xi_S=38nm
fe = supermag.get_material("Fe")   # E_ex=256meV, xi_F=0.7nm

# Compute Tc vs ferromagnet thickness
d_F = np.linspace(0.5, 20.0, 100)
Tc = supermag.critical_temperature(
    Tc0=nb["Tc"], d_S=50.0,
    d_F_array=d_F, E_ex=fe["E_ex"],
    xi_S=nb["xi_S"], xi_F=fe["xi_F"],
    gamma=0.15,
)

# Plot
plot_tc_vs_df(d_F, Tc, Tc0=nb["Tc"], save_path="tc_vs_df.png")
```

## Architecture

```
OCaml           ← pipeline orchestration, parameter sweeps, type-safe composition
    ↓ FFI (C calling convention via extern "C" headers)
C++             ← all numerical work: solvers, linear algebra, SIMD intrinsics
    ↑ built by
Makefile        ← builds OCaml, C++, Python wheels, runs tests/benchmarks
    ↓ packages into
Python wheel    ← what end users install and import
```

The C++ layer exposes **C-linkage headers only** (`cpp/include/supermag/`) as the FFI boundary. OCaml calls these via ctypes; Python calls them via pybind11.

## Research Domains

| Domain | Module | Status |
|--------|--------|--------|
| S/F Proximity Effect | `supermag.proximity` | ✅ Implemented |
| Usadel (diffusive) | `supermag.usadel` | 🔲 Stub |
| Eilenberger (clean limit) | `supermag.eilenberger` | 🔲 Stub |
| BdG (tight-binding) | `supermag.bdg` | 🔲 Stub |
| Ginzburg-Landau | `supermag.ginzburg_landau` | 🔲 Stub |
| Josephson Junctions | `supermag.josephson` | 🔲 Stub |
| Spin-Triplet SC | `supermag.triplet` | 🔲 Stub |

## Development Setup

For contributors building from source:

```bash
git clone https://github.com/seanstarks-hub/SUPERMag-Simulation-Suite.git
cd SUPERMag-Simulation-Suite

# Build C++ library
make build

# Run C++ tests
make test

# Run Python tests
pip install numpy matplotlib pytest
make pytest

# Run validation suite
make validate
```

### Requirements
- **C++17 compiler**: MSVC on Windows (`/std:c++17 /arch:AVX2`), GCC/Clang on Linux
- **Python ≥ 3.9** with NumPy ≥ 1.22, Matplotlib ≥ 3.5
- **OCaml ≥ 5.0** with dune (for pipeline development only)

## License

BSD-3-Clause — see [LICENSE](LICENSE).
