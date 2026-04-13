# Contributing to SUPERMag

## Prerequisites

- **Python ≥ 3.9** with NumPy ≥ 1.22, Matplotlib ≥ 3.5, pytest ≥ 7.0
- **C++17 compiler** (optional, for native acceleration): GCC/Clang on Linux, MSVC on Windows
- **CMake ≥ 3.15** (optional, pulled automatically by `pip install`)
- **OCaml ≥ 5.0** with dune (pipeline development only)

## Getting Started

```bash
git clone https://github.com/seanstarks-hub/SUPERMag-Simulation-Suite.git
cd SUPERMag-Simulation-Suite

# Install in development mode (pure-Python, no compiler needed)
pip install -e .

# Or install with native C++ acceleration (requires C++17 compiler + CMake)
pip install .

# Run Python tests
pytest python/tests/

# Run C++ tests (requires CMake)
make test

# Run validation suite against published reference data
make validate
```

All solvers have pure-Python fallbacks, so a C++ compiler is only needed for
native performance. CMake handles the entire C++ build — no manual `make build`
step is required.

## Project Layout

```
cpp/include/supermag/   C-linkage headers (FFI boundary)
cpp/src/                C++17 solver implementations
cpp/test/               C++ unit tests
python/supermag/        Python API (pybind11 bindings + pure-Python fallback)
  proximity.py          Tc(d_F) solver
  depairing.py          Depairing channels + optimizer/fitter utilities
  bdg.py … usadel.py    BdG, Eilenberger, GL, Josephson, triplet, Usadel
  _binding.cpp          pybind11 bridge (all solvers + depairing + optimizer)
python/tests/           Python tests (pytest)
ocaml/                  Pipeline orchestration (dune)
validation/             Comparison against published reference data
docs/theory/            Jupyter notebooks with derivations
docs/tutorials/         Usage tutorial notebooks
benchmarks/             Performance benchmarks
```

## Running Tests

| Command | What it runs |
|---------|-------------|
| `pytest python/tests/` | Python tests (proximity, solvers, sweeps, materials, themes, depairing) |
| `make test` | C++ unit tests via CMake (all 18 test suites) |
| `make validate` | Validation suite against Buzdin (1982), Radovic (1991), Ryazanov (2003), Bergeret (2005) reference data |
| `make bench` | Performance benchmarks |

## Adding a New Solver

1. Implement the C++ solver in `cpp/src/solvers/` (or `cpp/src/proximity/`) with a C-linkage header in `cpp/include/supermag/`.
2. Add pybind11 bindings in `python/supermag/_binding.cpp`.
3. Create or update the corresponding Python module in `python/supermag/` with `_USE_NATIVE` dispatch and a pure-Python fallback (see `depairing.py` for the canonical pattern).
4. Add C++ tests in `cpp/test/` and Python tests in `python/tests/`.
5. Register new public names in `python/supermag/__init__.py`.
6. If applicable, add a validation case in `validation/` with expected reference data.
7. Update `architecture.md` (equation registry, C API table, file listing).

## Code Style

- **C++**: C++17, prefer `const` and value semantics, SIMD kernels isolated in `cpp/src/linalg/`.
- **Python**: Follow existing conventions. NumPy-style docstrings. Type hints welcome.
- **Commits**: Short imperative subject line, blank line, then details if needed.

## Pull Requests

1. Fork the repository and create a feature branch.
2. Ensure `pytest python/tests/` passes (and `make test` if you changed C++ code).
3. Add or update tests for any changed functionality.
4. Open a PR with a clear description of what changed and why.
