# architecture.md — SUPERMag Binding Contract

> **This document is the single source of truth for all code generation,
> modification, and review in this repository.**
>
> Before writing or modifying ANY code in this repo, read this document
> in its entirety. Every function, formula, and interface described here
> already exists. Your job is to work within this system, never around it.

---

## §0 — Cardinal Rules

1. **Do not invent.** Every solver, module, header, and Python wrapper
   described in this document already exists. If you need a function,
   find it in the locations listed below. Do not create alternatives,
   duplicates, or "simplified versions."

2. **C++ is the source of truth for physics.** The C++ engine under
   `cpp/` defines the canonical formulas. Every Python fallback under
   `python/supermag/` MUST reproduce the identical formula — same
   terms, same scaling, same variable names. If they disagree, the
   Python fallback is wrong.

3. **Headers are the API contract.** The C function signatures in
   `cpp/include/supermag/*.h` are frozen. Python wrappers, pybind11
   bindings, and OCaml FFI stubs must match these signatures exactly.
   Do not add parameters, change types, or rename arguments.

4. **Every dependency must be declared.** If a Python module imports
   a package (scipy, numpy, matplotlib, cycler, etc.), that package
   MUST appear in `python/pyproject.toml` under `[project.dependencies]`
   or `[project.optional-dependencies]`.

5. **Validate inputs identically across layers.** If the C++ function
   checks `xi_F <= 0 → SUPERMAG_ERR_INVALID_DIM`, the Python wrapper
   MUST raise `ValueError` for the same condition before calling the
   fallback. Parameter validation tables are in §4.

---

## §1 — Layer Architecture

```
┌─────────────────────────────────────────────────────────┐
│  User Code / Notebooks / CLI                            │
├─────────────────────────────────────────────────────────┤
│  Python API         python/supermag/*.py                │
│    ├─ proximity.py    pair_amplitude, critical_temperature
│    ├─ usadel.py       solve                             │
│    ├─ eilenberger.py  solve                             │
│    ├─ bdg.py          solve                             │
│    ├─ ginzburg_landau.py  minimize                      │
│    ├─ josephson.py    current_phase_relation             │
│    ├─ triplet.py      solve                             │
│    ├─ sweeps.py       tc_parameter_sweep, tc_phase_diagram
│    ├─ materials.py    get_material, list_materials       │
│    ├─ plotting.py     plot_pair_amplitude, plot_tc_vs_df │
│    └─ themes/         apply_theme, theme_context, etc.   │
├─────────────────────────────────────────────────────────┤
│  pybind11 bridge    python/supermag/_binding.cpp         │
│    Module name: _native                                  │
│    Wraps: proximity (pair_amplitude, solve_tc_batch)     │
│    TODO: wrap usadel, eilenberger, bdg, gl, josephson,   │
│          triplet when C++ solvers are production-ready    │
├─────────────────────────────────────────────────────────┤
│  C++ Engine         cpp/src/**/*.cpp                     │
│    Headers:         cpp/include/supermag/*.h              │
│    Static library:  build/supermag.lib (.a on Linux)     │
├─────────────────────────────────────────────────────────┤
│  OCaml Orchestrator ocaml/lib/**/  ocaml/bin/            │
│    FFI stubs:       ocaml/lib/ffi/stubs.ml{i}            │
│    Pipeline:        ocaml/lib/pipeline/sweep.ml          │
│    CLI driver:      ocaml/bin/sweep_driver.ml            │
└─────────────────────────────────────────────────────────┘
```

Data flows DOWN only. Python calls C++ via `_native`. OCaml calls C++
via ctypes FFI. Python never calls OCaml. OCaml never calls Python.
C++ never calls up.

---

## §2 — Equation Registry

Every physics formula that appears in code is registered here with a
stable identifier. When implementing or modifying code, reference the
EQ-* identifier in a comment at the implementation site.

### EQ-1: Complex wave vector in ferromagnet
```
q = (1 + i) / ξ_F
```
- C++ impl: `cpp/src/proximity/kernels.cpp` → `compute_q()`
- Python impl: `python/supermag/proximity.py` → inline in fallback

### EQ-2: Proximity kernel, 0-junction (coth)
```
K₀(d_F, ξ_F) = q · coth(q · d_F) = q · cosh(q·d_F) / sinh(q·d_F)
```
- Phase enum: `SUPERMAG_PHASE_ZERO = 0`, Python `phase="zero"`
- C++ impl: `kernels.cpp` → `kernel_coth()`
- Python impl: `proximity.py` fallback, `phase="zero"` branch

### EQ-3: Proximity kernel, π-junction (tanh)
```
K_π(d_F, ξ_F) = q · tanh(q · d_F) = q · sinh(q·d_F) / cosh(q·d_F)
```
- Phase enum: `SUPERMAG_PHASE_PI = 1`, Python `phase="pi"`
- C++ impl: `kernels.cpp` → `kernel_tanh()`
- Python impl: `proximity.py` fallback, `phase="pi"` branch

### EQ-4: Thin-S self-consistency equation
```
F(T) = ln(Tc0/T) − Re[ψ(1/2 + α) − ψ(1/2)]
α = γ · K · Tc0 / (2π T) + λ_dep
```
- **No** η = ξ_S/d_S prefactor. γ absorbs coupling strength.
- C++ impl: `critical_temp.cpp` → `thin_s_equation()`
- Python impl: `proximity.py` fallback, `model="thin_s"` branch

### EQ-5: Fominov self-consistency equation (PRB 66, 014507)
```
F(T) = ln(Tc0/T) − Re[ψ(1/2 + α) − ψ(1/2)]
α = γ · K / (1 + γ_B · K) · Tc0 / (2π T) + λ_dep
```
- Reduces to EQ-4 when γ_B → 0.
- C++ impl: `critical_temp.cpp` → `fominov_determinant()`
- Python impl: `proximity.py` fallback, `model="fominov"` branch

### EQ-6: Pair amplitude in F layer
```
phase=zero:  F(x) = exp(−x/ξ_F) · cos(x/ξ_F)
phase=pi:    F(x) = exp(−x/ξ_F) · sin(x/ξ_F)
```
- C++ impl: `pair_amplitude.cpp`
- Python impl: `proximity.py` → `pair_amplitude()` fallback

### EQ-7: Depairing parameter (sum of all channels)
```
λ_dep = ag + zeeman + orbital + spin_orbit
```
- NULL depairing pointer → λ_dep = 0.0
- C++ impl: `depairing.cpp` → `supermag_depairing_total()`
- Python impl: `proximity.py` → inline sum in fallback

### EQ-8: Digamma (ψ) function, complex argument
```
Asymptotic: ψ(z) = ln(z) − 1/(2z) − Σ B_{2k}/(2k · z^{2k})
Recurrence: ψ(z) = ψ(z+1) − 1/z  (shift until |z| > 10)
```
- C++ impl: `cpp/src/common/digamma.cpp`
- Python fallback: `scipy.special.digamma` (MUST be in pyproject.toml)

### EQ-9: Josephson CPR (Buzdin, first harmonic)
```
I_c(d_F) ∝ exp(−d_F/ξ_F) · cos(d_F/ξ_F − π/4)
I(φ) = I_c · T_factor · sin(φ),  normalized to max|I|=1
```
- C++ impl: `josephson.cpp`
- Python impl: `josephson.py`

### EQ-10: BdG Nambu Hamiltonian
```
H_BdG = [[ −µ·I + E_ex·I − t·tridiag,    Δ·I        ],
          [ Δ*·I,                  +µ·I + E_ex·I + t·tridiag ]]
```
- Eigenvalues in meV (convert from eV after diagonalization)
- C++ impl: `bdg.cpp` (Jacobi eigensolver, n ≤ 500)
- Python impl: `bdg.py` (`numpy.linalg.eigvalsh`)

### EQ-11: GL TDGL relaxation
```
∂ψ/∂t = −αψ − β|ψ|²ψ + ξ²∇²ψ
Equilibrium: |ψ|² = −α/β  (uniform, α < 0)
```
- C++ impl: `ginzburg_landau.cpp` (double-buffered snapshot Euler)
- Python impl: `ginzburg_landau.py` (np.roll snapshot Euler)
- Both use snapshot Euler (consistent Laplacian reads).

### EQ-12: Triplet amplitude at non-collinear interface
```
f_↑↑(x) ∝ |sin(α)| · exp(−|x − x_int|/ξ_N)
α = magnetization_angles[i+1] − magnetization_angles[i]
```
- Defaults: ξ_F = 1.0 nm, ξ_N = 10.0 nm (configurable via API parameters)
- C++ impl: `triplet.cpp`
- Python impl: `triplet.py`

---

## §3 — C API Signatures

These are the canonical C function signatures from `cpp/include/supermag/`.
All wrappers (pybind11, Python fallback, OCaml FFI) must map to these
exactly. Do NOT add parameters, change return types, or rename.

**v0.2 changes:** `T` added to Usadel/Eilenberger, `Tc0` added to
Josephson, `mu` added to BdG, `xi_F`/`xi_N` added to Triplet.
New parameters use sentinel defaults (`T<=0` → 0.5·Tc0, `Tc0<=0` → 9.2,
`mu=0` → no shift, `xi_F/xi_N<=0` → legacy defaults).

```c
/* error.h */
const char* supermag_error_string(int code);

/* proximity.h */
double supermag_depairing_total(const supermag_depairing_t *dp);

int supermag_proximity_solve_tc(
    const supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double *tc_out);

int supermag_proximity_solve_tc_batch(
    const supermag_proximity_params_t *params,
    const double *d_F_array, int n_dF,
    const supermag_depairing_t *depairing,
    double *tc_out);

int supermag_proximity_pair_amplitude(
    double d_F, double xi_F, supermag_phase_t phase,
    int n_points, double *x_out, double *F_out);

/* usadel.h */
int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    double T,
    int n_grid, double* Delta_out, double* x_out);

/* eilenberger.h */
int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    double T,
    int n_grid, double* f_out, double* x_out);

/* bdg.h */
int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues);

/* ginzburg_landau.h */
int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    double* psi_real, double* psi_imag);

/* josephson.h */
int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
    int n_phases, double* phase_arr, double* current_out);

/* triplet.h */
int supermag_triplet_solve(
    int n_layers, const double* thicknesses,
    const double* magnetization_angles,
    double xi_F, double xi_N,
    int n_grid, double* f_triplet_out, double* x_out);
```

---

## §4 — Input Validation Contract

Every public function must validate inputs. The C++ layer returns error
codes; the Python layer raises exceptions. Both must reject the same
conditions.

| Parameter | Valid | C++ error code | Python exception |
|-----------|-------|---------------|-----------------|
| `d_F` | > 0 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `xi_F` | > 0 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `xi_S` | > 0 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `Tc0` | > 0 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `n_points` | ≥ 2 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `n_grid` | > 4 | `SUPERMAG_ERR_INVALID_DIM` | `ValueError` |
| `phase` | "zero" or "pi" | default ZERO | `ValueError` |
| `model` | "thin_s" or "fominov" | `SUPERMAG_ERR_INVALID_MODEL` | `ValueError` |
| `gamma_B` | ≥ 0 | (unchecked) | `ValueError` |
| null ptr | — | `SUPERMAG_ERR_NULL_PTR` | N/A (Python) |

---

## §5 — File Locations (Exhaustive)

### C++ Engine
```
cpp/include/supermag/
  proximity.h, usadel.h, eilenberger.h, bdg.h,
  ginzburg_landau.h, josephson.h, triplet.h,
  error.h, constants.h

cpp/src/common/
  error.cpp, constants.cpp, digamma.cpp, digamma.h

cpp/src/proximity/
  kernels.cpp, critical_temp.cpp, pair_amplitude.cpp, depairing.cpp

cpp/src/solvers/
  root_scalar.cpp, root_scalar.h, determinant.cpp,
  usadel.cpp, eilenberger.cpp, bdg.cpp,
  ginzburg_landau.cpp, josephson.cpp, triplet.cpp

cpp/src/linalg/
  tridiag.cpp, simd_kernels.cpp

cpp/test/
  test_proximity.cpp, test_tridiag.cpp, test_stubs.cpp,
  test_digamma.cpp, test_root_scalar.cpp, test_determinant.cpp
```

### Python Package
```
python/supermag/
  __init__.py, proximity.py, usadel.py, eilenberger.py,
  bdg.py, ginzburg_landau.py, josephson.py, triplet.py,
  materials.py, sweeps.py, plotting.py, _binding.cpp

python/supermag/themes/
  __init__.py, publication.py, presentation.py, draft.py, dark.py

python/tests/
  conftest.py, test_proximity.py, test_materials.py,
  test_solvers.py, test_sweeps.py, test_themes.py
```

### OCaml Orchestrator
```
ocaml/lib/ffi/       stubs.ml, stubs.mli
ocaml/lib/pipeline/  sweep.ml
ocaml/lib/types/     params.ml (types only)
ocaml/bin/           sweep_driver.ml
ocaml/test/          test_ffi.ml
```

### Build & CI
```
Makefile
config/detect_platform.mk, windows_x86_64.mk, linux_x86_64.mk, compiler_flags.mk
python/pyproject.toml, python/CMakeLists.txt
.github/workflows/ci.yml, wheels.yml, validation.yml
```

### Validation
```
validation/run_validation.py
validation/buzdin_1982/   params.json, expected/tc_vs_df.csv, expected/pair_amplitude.csv
validation/ryazanov_2003/ params.json, expected/tc_vs_df.csv
validation/radovic_1991/  (placeholder)
validation/bergeret_2005/ (placeholder)
```

---

## §6 — Known Bugs and Intentional Limitations

Document these so the agent does not "fix" them in the wrong direction
or silently reintroduce them.

### KNOWN-BUG-1: Python fallback kernel assignment is swapped
**Status:** RESOLVED.
Fixed: `proximity.py` fallback now correctly assigns coth to
phase="zero" and tanh to phase="pi", matching C++ engine.

### KNOWN-BUG-2: Python fallback has extra η = ξ_S/d_S scaling
**Status:** RESOLVED.
Fixed: `proximity.py` fallback no longer multiplies α by
`eta = xi_S / d_S`. Matches C++ (EQ-4, EQ-5).

### KNOWN-BUG-3: scipy not declared in pyproject.toml
**Status:** RESOLVED.
Fixed: `scipy>=1.8` is now listed in `[project.dependencies]`.

### KNOWN-LIMIT-1: Usadel/Eilenberger hardcode T = 0.5·Tc0
**Status:** RESOLVED.
Both C API signatures now accept a `T` parameter. When `T <= 0`,
the solver uses 0.5·Tc0 for backward compatibility.

### KNOWN-LIMIT-2: Josephson CPR is first-harmonic only
sin(φ) only. Higher harmonics needed for φ₀-junction states.

### KNOWN-LIMIT-3: BdG has no chemical potential (µ)
**Status:** RESOLVED.
C API now accepts `mu` parameter. On-site energy is `−µ ± E_ex`.
Legacy callers pass `mu = 0.0` to reproduce previous behaviour.

### KNOWN-LIMIT-4: Triplet hardcodes ξ_F=1, ξ_N=10 nm
**Status:** RESOLVED.
C API now accepts `xi_F` and `xi_N` parameters. When `<= 0`, the
solver uses the legacy defaults (1.0 nm and 10.0 nm).

### KNOWN-LIMIT-5: C++ GL solver uses in-place Gauss-Seidel
**Status:** RESOLVED.
C++ now uses double-buffered snapshot Euler, matching the Python
`np.roll` approach. Both produce consistent results.

### KNOWN-LIMIT-6: C++ Eilenberger/Usadel use stack arrays [2048]
Should be std::vector for thread safety.

### KNOWN-LIMIT-7: Fominov kernel is T-independent
The full Fominov model has T-dependent q_S·cot(q_S·d_S). Current
implementation uses simplified α = γK/(1+γ_B·K). Documented in
comments of critical_temp.cpp but not implemented.

### KNOWN-LIMIT-8: Josephson CPR uses hardcoded Tc_ref when Tc0 not supplied
**Status:** RESOLVED.
C API now accepts `Tc0` parameter. When `Tc0 <= 0`, the solver falls
back to 9.2 K (Nb). Python wrapper always passes through user-supplied
`Tc0` value.

---

## §7 — Prohibited Actions

When generating or modifying code in this repository, the following
are explicitly forbidden:

1. **Do not create new solver files.** All solver modules exist.
   If something is a stub, extend the existing file.

2. **Do not add physics scaling factors** (η, prefactors, unit
   conversions) to Python fallbacks unless they appear in the
   corresponding C++ implementation AND the EQ-* registry above.

3. **Do not change C header signatures.** The `extern "C"` API in
   `cpp/include/supermag/*.h` is frozen for ABI stability.

4. **Do not duplicate solver logic.** If the Python fallback and C++
   engine implement the same formula, they must reference the same
   EQ-* equation. There is exactly one formula per physical quantity.

5. **Do not hardcode material constants** (Tc, ξ, E_ex) inside solver
   functions. Use the parameters passed by the caller. Exception:
   josephson.cpp Tc_ref=9.2 is KNOWN-LIMIT (documented in §6).

6. **Do not use `else` as a catch-all for phase/model selection.**
   Always explicitly match known values and raise/return an error
   for unrecognized input.

7. **Do not add `extern "C"` wrappers around headers that already
   contain them.** All headers in `cpp/include/supermag/` have their
   own `extern "C"` guards.

---

## §8 — How to Use This Document

### When modifying existing code:
1. Find the relevant EQ-* equation(s) in §2
2. Verify your change preserves the formula exactly
3. Check §4 for input validation requirements
4. Check §6 to see if you're touching a known bug or limitation
5. If modifying Python fallback: confirm formula matches C++ impl

### When adding a new feature:
1. Add a new EQ-* entry to §2 with the formula, C++ location,
   and Python location
2. Implement in C++ first (source of truth)
3. Add Python fallback that uses the identical formula
4. Add to §3 if it's a new C API function
5. Add validation entries to §4
6. Register any new dependencies in pyproject.toml

### When reviewing code:
1. For each formula in the diff, find its EQ-* entry
2. Verify C++ and Python implementations are identical
3. Check that phase="zero"→coth and phase="pi"→tanh (EQ-2, EQ-3)
4. Verify no unregistered scaling factors
5. Confirm input validation matches §4