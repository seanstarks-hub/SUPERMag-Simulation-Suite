# supermag.usadel

Usadel diffusive-limit solver for S/F bilayers. Computes the self-consistent
superconducting order parameter $\Delta(x)$ using linearized Usadel equations
with Matsubara frequency summation and interface matching.

**C++ engine:** `usadel.cpp`

---

## `usadel.solve()`

Solve the Usadel equation for an S/F bilayer.

```python
supermag.usadel.solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex, n_grid=200, T=None)
```

> **C++ note:** The C++ engine supports both `LINEARIZED` and `NONLINEAR`
> Usadel modes via a `mode` enum (EQ-16). The Python API currently uses
> linearized mode only.

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `Tc0` | float | *required* | Bulk superconductor critical temperature (K) |
| `d_S` | float | *required* | Superconductor layer thickness (nm) |
| `d_F` | float | *required* | Ferromagnet layer thickness (nm) |
| `xi_S` | float | *required* | Superconductor coherence length (nm) |
| `xi_F` | float | *required* | Ferromagnet coherence length (nm) |
| `E_ex` | float | *required* | Exchange energy in ferromagnet (meV) |
| `n_grid` | int | `200` | Number of spatial grid points |
| `T` | float or None | `None` | Temperature (K). If `None`, defaults to $0.5 \cdot T_{c0}$ |

### Returns

| Type | Description |
|------|-------------|
| `x` : `numpy.ndarray` | Position array (nm), spanning $[-d_S, d_F]$ |
| `Delta` : `numpy.ndarray` | Order parameter profile $\Delta(x)$ (meV) |

### Example

```python
import supermag

x, Delta = supermag.usadel.solve(
    Tc0=9.2, d_S=50.0, d_F=10.0,
    xi_S=38.0, xi_F=0.7, E_ex=256.0, T=4.0
)
# Delta is BCS-like in S region, decays with oscillation into F
```
