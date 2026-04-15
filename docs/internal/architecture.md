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
│    ├─ depairing.py    depairing_*, optimize_tc, inverse_tc, fit_tc
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
│    Wraps: proximity, bdg, usadel, eilenberger, gl,       │
│           josephson, triplet, depairing (6), optimizer (3)│
├─────────────────────────────────────────────────────────┤
│  C++ Engine         cpp/src/**/*.cpp                     │
│    Headers:         cpp/include/supermag/*.h              │
│    Static library:  build/supermag.lib (.a on Linux)     │
├─────────────────────────────────────────────────────────┤
│  OCaml Orchestrator ocaml/lib/**/  ocaml/bin/            │
│    FFI stubs:       ocaml/lib/ffi/stubs.ml{i} (18 C fn) │
│    Typed solvers:   ocaml/lib/ffi/solvers.ml             │
│    Pipeline:        ocaml/lib/pipeline/sweep.ml,chain.ml │
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
- Diverges as d_F→0 → α = γ/(γ_B + K) → 0 → Tc → Tc0 ✓
- Phase enum: `SUPERMAG_PHASE_ZERO = 0`, Python `phase="zero"`
- C++ impl: `kernels.cpp` → `kernel_coth()`
- Python impl: `proximity.py` fallback, `phase="zero"` branch

### EQ-3: Proximity kernel, π-junction (tanh)
```
K_π(d_F, ξ_F) = q · tanh(q · d_F) = q · sinh(q·d_F) / cosh(q·d_F)
```
- Vanishes as d_F→0 (no F-layer → no pair-breaking)
- Phase enum: `SUPERMAG_PHASE_PI = 1`, Python `phase="pi"`
- C++ impl: `kernels.cpp` → `kernel_tanh()`
- Python impl: `proximity.py` fallback, `phase="pi"` branch

### EQ-4: Thin-S self-consistency equation
```
F(T) = ln(Tc0/T) − Re[ψ(1/2 + α) − ψ(1/2)]
α = γ / (γ_B + K) · Tc0 / (2π T) + λ_dep
```
- K in denominator: as d_F→0, K(coth)→∞, α→0, Tc→Tc0.
- γ is the coupling strength, γ_B is the interface barrier.
- C++ impl: `critical_temp.cpp` → `thin_s_equation()`
- Python impl: `proximity.py` fallback, `model="thin_s"` branch

### EQ-5: Fominov self-consistency equation (PRB 66, 014507)
```
F(T) = ln(Tc0/T) − Re[ψ(1/2 + α) − ψ(1/2)]
α = γ / (γ_B + K + Ω_S(T)) · Tc0 / (2π T) + λ_dep

Ω_S(T) = √(T/Tc0) · coth(√(T/Tc0) · d_S/ξ_S)
```
- T-dependent S-layer impedance Ω_S(T) accounts for finite-thickness
  superconductor. At T→0: Ω_S→0, recovers simplified formula.
  Thick-S limit (d_S≫ξ_S): Ω_S→√(T/Tc0). Thin-S limit: Ω_S→ξ_S/d_S.
- Reduces to EQ-4 when γ_B → 0 and Ω_S → 0.
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

### EQ-7A: Abrikosov-Gorkov depairing (spin-flip scattering)
```
λ_AG = Γ_s / (2 · k_B · T)
```
- Γ_s: spin-flip scattering rate (SI: Joules)
- C++ impl: `depairing.cpp` → `supermag_depairing_compute()`

### EQ-7B: Zeeman depairing (Pauli paramagnetic limit)
```
λ_Z = (μ_B · H)² / (2π · k_B · T)²
```
- H: applied magnetic field (Tesla)
- C++ impl: `depairing.cpp` → `supermag_depairing_compute()`

### EQ-7C: Orbital depairing (thin-film limit)
```
λ_orb = D · (e·H)² · d² / (3 · ℏ² · 2π · k_B · T)
```
- D: diffusion coefficient, d: film thickness
- C++ impl: `depairing.cpp` → `supermag_depairing_compute()`

### EQ-7D: Spin-orbit depairing
```
λ_SO = Γ_so / (2 · k_B · T)
```
- Γ_so: spin-orbit scattering rate (SI: Joules)
- C++ impl: `depairing.cpp` → `supermag_depairing_compute()`

### EQ-8: Digamma (ψ) function, complex argument
```
Asymptotic: ψ(z) = ln(z) − 1/(2z) − Σ B_{2k}/(2k · z^{2k})
Recurrence: ψ(z) = ψ(z+1) − 1/z  (shift until |z| > 10)
```
- C++ impl: `cpp/src/common/digamma.cpp`
- Python fallback: `scipy.special.digamma` (MUST be in pyproject.toml)

### EQ-9: Josephson CPR (Matsubara frequency sum)
```
I(φ) = T · Σ_n Re[ P_n · Δ²·sin(φ) / √((ω_n² + Δ²sin²(φ/2))·(ω_n² + Δ²)) ]

P_n = exp(−q_n·d_F) · exp(iπ/4)           [F-layer propagator]
q_n = √(2(ω_n/E_ex + i)) / ξ_F            [complex wave vector]
ω_n = π·k_B·T·(2n+1)                      [Matsubara frequency]
Δ(T) = 1.764·k_B·Tc0·√(1−T/Tc0)
```
- Cutoff: ω_n > 20Δ or n > 500.
- At T→0, n=0: recovers Buzdin first-harmonic sin(φ).
- sin²(φ/2) in denominator generates higher harmonics.
- Complex P_n phase produces 0-π oscillation with d_F/ξ_F.
- Normalized to max|I|=1.
- C++ impl: `josephson.cpp`
- Python impl: `josephson.py`

### EQ-10: BdG Nambu Hamiltonian
```
H_BdG = [[ −µ·I + E_ex·I − t·tridiag,    Δ·I        ],
          [ Δ*·I,                  +µ·I + E_ex·I + t·tridiag ]]
```
- t_hop converted from eV to meV internally (×1e3); all computation in meV
- Cyclic Jacobi eigensolver (systematic row-major sweeps, n ≤ 2000)
- Optional eigenvector output via `eigenvectors_out` parameter (NULL to skip)
- C++ impl: `bdg.cpp`
- Python impl: `bdg.py` (`numpy.linalg.eigvalsh`)

### EQ-11: GL TDGL relaxation
```
∂ψ/∂t = −αψ − β|ψ|²ψ + ξ²∇²ψ
Equilibrium: |ψ|² = −α/β  (uniform, α < 0)
ξ² = 1/(2κ²)   (coherence length from GL parameter)
```
- C++ impl: `ginzburg_landau.cpp` (double-buffered snapshot Euler)
  Mode enum selects SCALAR (no gauge) or GAUGE (EQ-18) with applied field.
- Python impl: `ginzburg_landau.py` (np.roll snapshot Euler)
- Both use snapshot Euler (consistent Laplacian reads).

### EQ-12: Triplet amplitude at non-collinear interface
```
f_↑↑(x) ∝ |sin(α)| · exp(−|x − x_int|/ξ_N)
α = magnetization_angles[i+1] − magnetization_angles[i]
```
- Defaults: ξ_F = 1.0 nm, ξ_N = 10.0 nm (configurable via API parameters)
- C++ now uses Usadel-based model (EQ-19) with dual triplet contributions
- Temperature dependence via required `T` parameter and `mode` enum
- C++ impl: `triplet.cpp`
- Python impl: `triplet.py`

### EQ-16: Nonlinear Usadel self-consistency (tridiagonal)
```
Δ(x) self-consistent via Newton-linearized finite-difference:
  D·∂²Δ/∂x² − (ω_n + 1/(2τ_sf))·Δ + g·Δ_BCS = 0   (in S)
  D·∂²Δ/∂x² − (ω_n + E_ex/ℏ)·Δ = 0                  (in F)

Discretized: A·δΔ = −F   (tridiag system, Thomas algorithm)
Matsubara sum: ω_n = π·k_B·T·(2n+1)  cut at ω_n > 20Δ or n > 500
Self-consistency: adaptive mixing (0.3 → 0.5), 100 iter, tol 1e-8
```
- Proportional grid splitting: n_S = n·d_S/(d_S+d_F)
- Uses `supermag_tridiag_solve()` for the linearized Newton step
- C++ impl: `usadel.cpp`

### EQ-17: Riccati-parametrized Eilenberger (RK4)
```
da/dx = −2iω_n·a/v_F + Δ(x)/v_F − Δ*(x)/v_F · a²

Riccati variable: a(x) relates to f, g via
  f = 2a/(1+a·ã),  g = (1−a·ã)/(1+a·ã)

RK4 integration: 4th-order Runge-Kutta on Riccati ODE
Multi-frequency: sum over Matsubara frequencies ω_n
```
- Normalization constraint |a| ≤ 1 (physical bound)
- Forward sweep for right-movers, backward sweep for left-movers
- Proportional grid splitting: n_S = n·d_S/(d_S+d_F)
- C++ impl: `eilenberger.cpp`

### EQ-18: GL with Maxwell self-consistent gauge field
```
∂ψ/∂t = −α·ψ − β·|ψ|²·ψ + ξ²·∇²_A ψ

∇²_A: gauge-covariant Laplacian with Peierls phase factors
  ψ_{i±1} → ψ_{i±1}·exp(∓i·A·dx)

Landau gauge: A_y = H·x (uniform applied field)
Self-consistent field update: A_y ← A_y − η·J_y (every 50 steps)
J_y = Im[ψ*·(∂/∂y − iA_y)ψ]
```
- κ defines coherence length: ξ² = 1/(2κ²)
- Adaptive convergence: max 5000 steps, tol 1e-8, 3 consecutive checks
- C++ impl: `ginzburg_landau.cpp`

### EQ-19: Triplet coupled Usadel model
```
Singlet: f_0(x) ∝ Δ(T)/ℏ · exp(−q_F·x_int)
  q_F = (1+i)/ξ_F  (complex, oscillatory decay in F)

Triplet: f_1(x) = |sin(Δα)| · [c_LR·exp(−x/ξ_N) + c_SR·exp(−x/ξ_F)]
  Δα = misalignment angle between adjacent magnetizations
  c_LR: long-range triplet amplitude (ξ_N decay)
  c_SR = 0.3·c_LR: short-range S_z=0 triplet (ξ_F decay)

Temperature: Δ(T) = 1.764·k_B·Tc0·√(1−T/Tc0)
```
- Singlet-to-triplet conversion at magnetically inhomogeneous interfaces
- Long-range + short-range dual triplet contributions
- C++ impl: `triplet.cpp`

### EQ-13: S/N/F trilayer effective kernel
```
K_SNF = q_N · (K_F + q_N·tanh(q_N·d_N)) / (q_N + K_F·tanh(q_N·d_N))

q_N = 1/ξ_N   (real, no exchange in normal metal)
K_F = kernel_coth or kernel_tanh (phase-dependent bilayer kernel)
```
- Continued-fraction composition of N-layer propagation with F-layer kernel.
- Limits: d_N→0 recovers bilayer; d_N→∞ gives K→q_N (S sees only N metal).
- C++ impl: `kernel_snf.cpp` → `supermag_proximity_kernel_snf()`

### EQ-14: Graded ferromagnet effective kernel
```
M_total = ∏_{i=1}^{N} M_i    (transfer matrix cascade, vacuum → S interface)

M_i = [[cosh(q_i·δ),  sinh(q_i·δ)/q_i],
       [q_i·sinh(q_i·δ),  cosh(q_i·δ)  ]]

q_i = (1+i)/ξ_F(x_i),   ξ_F(x) = ξ_F_ref · √(E_ref / E_ex(x))

K_coth = M[0][0] / M[0][1]   (0-junction)
K_tanh = M[1][0] / M[0][0]   (π-junction)
```
- E_ex(x) profiles: LINEAR, EXPONENTIAL, STEP.
- Uniform E_ex recovers standard bilayer kernel.
- C++ impl: `kernel_graded.cpp` → `supermag_proximity_kernel_graded()`

### EQ-15: Magnetic domain effective kernel
```
Domain structure: alternating ±E_ex magnetization.
+E_ex: q = (1+i)/ξ_F     −E_ex: q = (1−i)/ξ_F

Transfer matrix cascade through N domains (+ optional walls).
Same extraction rules as EQ-14.
```
- Sharp walls (domain_wall=0): adjacent slices with flipped q.
- Finite walls: linearly interpolated exchange energy across wall region.
- Single domain recovers standard bilayer kernel.
- C++ impl: `kernel_domains.cpp` → `supermag_proximity_kernel_domains()`

### EQ-20: Golden-section optimizer for d_F
```
Minimize |Tc(d_F) − Tc_target| over [d_F_lo, d_F_hi] using golden-section search.
φ = (√5+1)/2;  c = b−(b−a)/φ;  d = a+(b−a)/φ
Compare |Tc(c)−target| vs |Tc(d)−target| to narrow bracket.
```
- Converges in ~100 iterations to Δd_F < 1e-10 nm.
- C++ impl: `optimizer.cpp` → `supermag_optimize_tc()`
- Python impl: `depairing.py` → `optimize_tc()` fallback

### EQ-21: Brent inverse solver for d_F
```
Find root of F(d_F) = Tc(d_F) − Tc_target = 0 via Brent's method.
Inverse quadratic interpolation + bisection fallback.
```
- Requires sign change: F(d_F_lo)·F(d_F_hi) < 0.
- If no sign change, returns endpoint with smaller |F|.
- C++ impl: `optimizer.cpp` → `supermag_inverse_tc()`
- Python impl: `depairing.py` → `inverse_tc()` fallback (bisection)

### EQ-22: Nelder-Mead least-squares fit
```
χ² = Σᵢ (Tc_calc(d_F_i; θ) − Tc_data_i)²
θ = {γ, γ_B, E_ex, ξ_F}  (any subset flagged for fitting)
```
- Pure simplex optimization (derivative-free).
- 4 fit flags independently toggle which parameters to vary.
- C++ impl: `optimizer.cpp` → `supermag_fit_tc()` (custom Nelder-Mead)
- Python impl: `depairing.py` → `fit_tc()` fallback (`scipy.optimize.minimize`)

---

## §3 — C API Signatures

These are the canonical C function signatures from `cpp/include/supermag/`.
All wrappers (pybind11, Python fallback, OCaml FFI) must map to these
exactly. Do NOT add parameters, change return types, or rename.

**v0.2 changes:** `T` added to Usadel/Eilenberger, `Tc0` added to
Josephson, `mu` added to BdG, `xi_F`/`xi_N` added to Triplet.

**v0.3 changes (corrective):**
- Eliminated all `_ext` variants. Advanced parameters merged into main signatures.
- `model` enum split into orthogonal `model` (equation) + `geometry` (structure) enums.
- Added mode enums: `supermag_usadel_mode_t`, `supermag_gl_mode_t`, `supermag_triplet_mode_t`.
- Sentinel `T <= 0` defaults removed. T is now required > 0 for Usadel, Eilenberger, Triplet.
- Added individual depairing functions and optimizer/inverse/fit utilities.
- Added `spin_active` interface parameter to proximity params.
- Kernel overflow threshold lowered from 350 to 3.0 for earlier asymptotic switch.

**v0.4 changes:**
- Added `supermag_solver_options_t` composable options struct (`solver_options.h`).
- Usadel, Eilenberger, Josephson, GL signatures now accept `const supermag_solver_options_t *opts` (NULL = defaults).
- Triplet signature reordered: `E_ex_per_layer`/`D_per_layer` before `xi_F`/`xi_N`; added `Tc0` parameter.
- Domain wall interpolation (`WALL_SLICES=10`) for `kernel_domains`; `domain_wall` field added to `supermag_domain_params_t`.

```c
/* error.h */
const char* supermag_error_string(int code);

/* proximity.h — enums */
typedef enum { SUPERMAG_MODEL_THIN_S=0, SUPERMAG_MODEL_FOMINOV=1,
               SUPERMAG_MODEL_FOMINOV_MULTI=2 } supermag_model_t;
typedef enum { SUPERMAG_GEOM_BILAYER=0, SUPERMAG_GEOM_TRILAYER=1,
               SUPERMAG_GEOM_GRADED=2, SUPERMAG_GEOM_DOMAINS=3 } supermag_geometry_t;
typedef enum { SUPERMAG_PHASE_ZERO=0, SUPERMAG_PHASE_PI=1 } supermag_phase_t;

/* proximity.h — params struct */
typedef struct {
    double Tc0, d_S, d_F, xi_S, xi_F, gamma, gamma_B, E_ex, D_F, D_S;
    supermag_model_t model;
    supermag_phase_t phase;
    supermag_geometry_t geometry;
    const void *geom_params;
    const supermag_spin_active_t *spin_active;
} supermag_proximity_params_t;

/* proximity.h */
double supermag_depairing_total(const supermag_depairing_t *dp);

int supermag_depairing_compute(
    const supermag_depairing_input_t *input,
    supermag_depairing_t *output);

/* depairing.h — individual channel functions */
double supermag_depairing_ag(double gamma_s_meV, double T_kelvin);
double supermag_depairing_zeeman(double H_tesla, double T_kelvin);
double supermag_depairing_orbital_perp(double D_nm2ps, double H_tesla,
                                       double thickness_nm, double T_kelvin);
double supermag_depairing_orbital_par(double D_nm2ps, double H_tesla,
                                      double thickness_nm, double T_kelvin);
double supermag_depairing_soc(double Gamma_so_meV, double T_kelvin);

int supermag_depairing_from_physical(
    double gamma_s_meV, double H_tesla,
    double D_nm2ps, double thickness_nm,
    double Gamma_so_meV, double T_kelvin,
    supermag_depairing_t *output);

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

int supermag_proximity_kernel_snf(...);
int supermag_proximity_kernel_graded(...);
int supermag_proximity_kernel_domains(...);

/* optimizer.h */
int supermag_optimize_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double d_F_lo, double d_F_hi,
    double Tc_target,
    double *d_F_out);

int supermag_inverse_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double Tc_target,
    double d_F_lo, double d_F_hi,
    double *d_F_out);

int supermag_fit_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    const double *d_F_data, const double *Tc_data, int n_data,
    int fit_gamma, int fit_gamma_B, int fit_E_ex, int fit_xi_F,
    double *chi2_out);

/* solver_options.h */
typedef struct {
    int    matsubara_max;       /* Default: 500   */
    double omega_cut_factor;    /* Default: 20.0  */
    int    max_steps;           /* Default: 5000  */
    int    max_iter;            /* Default: 100   */
    double conv_tol;            /* Default: 1e-8  */
    int    root_grid_points;    /* Default: 1000  */
} supermag_solver_options_t;

supermag_solver_options_t supermag_default_solver_options(void);

/* usadel.h */
typedef enum { SUPERMAG_USADEL_LINEARIZED=0, SUPERMAG_USADEL_NONLINEAR=1
             } supermag_usadel_mode_t;

int supermag_usadel_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double xi_F, double E_ex,
    double T,
    supermag_usadel_mode_t mode,
    const supermag_solver_options_t *opts,
    int n_grid, double* Delta_out, double* x_out);

/* eilenberger.h */
int supermag_eilenberger_solve(
    double Tc0, double d_S, double d_F,
    double xi_S, double E_ex,
    double T,
    const supermag_solver_options_t *opts,
    int n_grid, double* f_out, double* x_out);

/* bdg.h */
int supermag_bdg_solve(
    int n_sites, double t_hop, double Delta, double E_ex, double mu,
    double* eigenvalues_out, int* n_eigenvalues,
    double* eigenvectors_out);

/* ginzburg_landau.h */
typedef enum { SUPERMAG_GL_SCALAR=0, SUPERMAG_GL_GAUGE=1
             } supermag_gl_mode_t;

int supermag_gl_minimize(
    double alpha, double beta, double kappa,
    int nx, int ny, double dx,
    supermag_gl_mode_t mode, double H_applied,
    const supermag_solver_options_t *opts,
    double* psi_real, double* psi_imag);

/* josephson.h */
int supermag_josephson_cpr(
    double d_F, double xi_F, double E_ex, double T, double Tc0,
    double gamma_B,
    int n_phases, const double* phase_arr,
    const supermag_solver_options_t *opts,
    double* current_out,
    double* Ic_out);

/* triplet.h */
typedef enum { SUPERMAG_TRIPLET_PHENOMENOLOGICAL=0,
               SUPERMAG_TRIPLET_USADEL=1 } supermag_triplet_mode_t;

int supermag_triplet_solve(
    int n_layers, const double* thicknesses,
    const double* magnetization_angles,
    const double* E_ex_per_layer, const double* D_per_layer,
    double xi_F, double xi_N, double T,
    double Tc0,
    supermag_triplet_mode_t mode,
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
| `model` | "thin_s", "fominov", or "fominov_multi" | `SUPERMAG_ERR_INVALID_MODEL` | `ValueError` |
| `gamma_B` | ≥ 0 | (unchecked) | `ValueError` |
| null ptr | — | `SUPERMAG_ERR_NULL_PTR` | N/A (Python) |

---

## §5 — File Locations (Exhaustive)

### C++ Engine
```
cpp/include/supermag/
  proximity.h, usadel.h, eilenberger.h, bdg.h,
  ginzburg_landau.h, josephson.h, triplet.h,
  error.h, constants.h, solver_options.h

cpp/src/common/
  error.cpp, constants.cpp, digamma.cpp, digamma.h,
  solver_options.cpp

cpp/src/proximity/
  kernels.cpp, critical_temp.cpp, pair_amplitude.cpp, depairing.cpp,
  transfer_matrix.cpp, transfer_matrix.h,
  kernel_snf.cpp, kernel_graded.cpp, kernel_domains.cpp,
  spin_active.cpp, depairing_models.cpp, optimizer.cpp

cpp/src/solvers/
  root_scalar.cpp, root_scalar.h, determinant.cpp,
  usadel.cpp, eilenberger.cpp, bdg.cpp,
  ginzburg_landau.cpp, josephson.cpp, triplet.cpp

cpp/src/linalg/
  tridiag.cpp, simd_kernels.cpp

cpp/test/
  test_proximity.cpp, test_tridiag.cpp, test_stubs.cpp,
  test_digamma.cpp, test_root_scalar.cpp, test_determinant.cpp,
  test_transfer_matrix.cpp, test_kernel_snf.cpp,
  test_kernel_graded.cpp, test_kernel_domains.cpp,
  test_usadel.cpp, test_eilenberger.cpp, test_bdg.cpp,
  test_ginzburg_landau.cpp, test_josephson.cpp, test_triplet.cpp,
  test_depairing_models.cpp, test_optimizer.cpp
```

### Python Package
```
python/supermag/
  __init__.py, proximity.py, depairing.py, usadel.py, eilenberger.py,
  bdg.py, ginzburg_landau.py, josephson.py, triplet.py,
  materials.py, sweeps.py, plotting.py, _binding.cpp

python/supermag/themes/
  __init__.py, publication.py, presentation.py, draft.py, dark.py

python/tests/
  conftest.py, test_proximity.py, test_materials.py,
  test_solvers.py, test_sweeps.py, test_themes.py, test_depairing.py,
  test_plotting.py
```

### OCaml Orchestrator
```
ocaml/lib/ffi/       stubs.ml, stubs.mli, solvers.ml
ocaml/lib/pipeline/  sweep.ml, chain.ml
ocaml/lib/types/     params.ml, material.ml, result.ml, geometry.ml
ocaml/bin/           sweep_driver.ml
ocaml/test/          test_ffi.ml, test_chain.ml, test_sweep.ml
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
validation/radovic_1991/  params.json, expected/tc_vs_df.csv
validation/bergeret_2005/ params.json, expected/triplet_amplitude.csv
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
Both C API signatures now accept a `T` parameter. T must be > 0;
sentinel defaults have been removed (T <= 0 returns SUPERMAG_ERR_INVALID_DIM).

### KNOWN-LIMIT-2: Josephson CPR is first-harmonic only
**Status:** RESOLVED.
Replaced first-harmonic sin(φ) with Matsubara frequency sum (EQ-9).
Higher harmonics arise naturally from sin²(φ/2) denominator.

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
**Status:** RESOLVED.
Both `usadel.cpp` and `eilenberger.cpp` now use `std::vector`.
Static `const double` GL quadrature tables are compile-time constants.

### KNOWN-LIMIT-7: Fominov kernel is T-independent
**Status:** RESOLVED.
Fominov determinant now includes T-dependent S-layer impedance
Ω_S(T) = √(T/Tc0)·coth(√(T/Tc0)·d_S/ξ_S) in the α denominator
(EQ-5). Both C++ and Python implementations updated.

### KNOWN-LIMIT-8: Josephson CPR uses hardcoded Tc_ref when Tc0 not supplied
**Status:** RESOLVED.
C API now accepts `Tc0` parameter. When `Tc0 <= 0`, the solver falls
back to 9.2 K (Nb). Python wrapper always passes through user-supplied
`Tc0` value.

### KNOWN-BUG-4: GL coherence length depended on grid spacing
**Status:** RESOLVED (Phase 2).
`xi2 = dx*dx` was incorrect — made ξ depend on grid spacing instead
of physics. Fixed to `xi2 = 1.0 / (2.0 * kappa * kappa)`.
This is a physics-affecting change; GL outputs now depend on κ as intended.

### KNOWN-BUG-5: BdG mixed eV and meV units
**Status:** RESOLVED (Phase 2).
Previously t_hop was in eV while other parameters in meV. Now t_hop
is converted from eV to meV internally (×1e3). All computation and
eigenvalue output is in meV.

### KNOWN-LIMIT-9: Usadel uses analytic profiles only
**Status:** RESOLVED (Phase 2).
Replaced analytic cosh/exp profiles with full nonlinear Newton-linearized
tridiagonal finite-difference solver (EQ-16). Uses `supermag_tridiag_solve()`.

### KNOWN-LIMIT-10: Eilenberger uses Forward Euler
**Status:** RESOLVED (Phase 2).
Replaced Forward Euler with 4th-order Runge-Kutta integration (EQ-17).
Added Matsubara frequency summation (was single n=0 only).
Fixed Riccati clamping from |a|≤2 to |a|≤1.

### KNOWN-LIMIT-11: BdG Jacobi sweep is naive find-max
**Status:** RESOLVED (Phase 2).
Replaced naive find-max-off-diagonal with cyclic Jacobi (systematic
row-major sweeps). Cap raised from 500 to 2000 sites.

### KNOWN-LIMIT-12: Triplet uses purely phenomenological model
**Status:** RESOLVED (Phase 2).
Replaced exponential decay with Usadel-based model (EQ-19) with
complex singlet amplitude, dual triplet contributions, and temperature
dependence via `supermag_triplet_solve_ext()`.

---

## §7 — Prohibited Actions

When generating or modifying code in this repository, the following
are explicitly forbidden:

1. **Do not create new solver files.** All solver modules exist.
   If something is a stub, extend the existing file.

2. **Do not add physics scaling factors** (η, prefactors, unit
   conversions) to Python fallbacks unless they appear in the
   corresponding C++ implementation AND the EQ-* registry above.

3. **Extend C header signatures with care.** The `extern "C"` API in
   `cpp/include/supermag/*.h` should remain stable. New parameters
   (e.g., `const supermag_solver_options_t *opts`) may be appended
   when defaults preserve backward compatibility (NULL = legacy
   behaviour). Document all signature changes in the §3 version log.

4. **Do not duplicate solver logic.** If the Python fallback and C++
   engine implement the same formula, they must reference the same
   EQ-* equation. There is exactly one formula per physical quantity.

5. **Do not hardcode material constants** (Tc, ξ, E_ex) inside solver
   functions. Use the parameters passed by the caller.

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