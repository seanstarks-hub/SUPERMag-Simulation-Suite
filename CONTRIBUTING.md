# Contributing to SUPERMag

## Prerequisites

- **C++17 compiler**: MSVC on Windows (`/std:c++17 /arch:AVX2`), GCC/Clang on Linux
- **Python ≥ 3.9** with NumPy ≥ 1.22, Matplotlib ≥ 3.5, pytest ≥ 7.0
- **OCaml ≥ 5.0** with dune (pipeline development only)

## Building from Source

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

# Run validation suite against published reference data
make validate
```

## Project Layout

```
cpp/include/supermag/   C-linkage headers (FFI boundary)
cpp/src/                C++17 solver implementations
cpp/test/               C++ unit tests
python/supermag/        Python API (pybind11 bindings + pure-Python fallback)
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
| `make test` | C++ unit tests (proximity, tridiag, digamma, determinant, root_scalar, stubs) |
| `make pytest` | Python tests via pytest |
| `make validate` | Validation suite against Buzdin (1982), Ryazanov (2003) reference data |
| `make bench` | Performance benchmarks |

## Adding a New Solver

1. Implement the C++ solver in `cpp/src/solvers/` with a C-linkage header in `cpp/include/supermag/`.
2. Add pybind11 bindings in `python/supermag/_binding.cpp`.
3. Replace the `NotImplementedError` in the corresponding `python/supermag/` stub module.
4. Add C++ tests in `cpp/test/` and Python tests in `python/tests/`.
5. If applicable, add a validation case in `validation/` with expected reference data.

## Code Style

- **C++**: C++17, prefer `const` and value semantics, SIMD kernels isolated in `cpp/src/linalg/`.
- **Python**: Follow existing conventions. NumPy-style docstrings. Type hints welcome.
- **Commits**: Short imperative subject line, blank line, then details if needed.

## Pull Requests

1. Fork the repository and create a feature branch.
2. Ensure `make test` and `make pytest` pass.
3. Add or update tests for any changed functionality.
4. Open a PR with a clear description of what changed and why.
