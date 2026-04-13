# supermag.eilenberger

Eilenberger clean-limit solver for S/F bilayers. Integrates the Riccati ODE
for the anomalous Green's function with Fermi-surface averaging.

**C++ engine:** `eilenberger.cpp`

---

## `eilenberger.solve()`

Solve the Eilenberger equation for an S/F bilayer in the clean limit.

```python
supermag.eilenberger.solve(Tc0, d_S, d_F, xi_S, E_ex, n_grid=200, T=None)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `Tc0` | float | *required* | Bulk superconductor critical temperature (K) |
| `d_S` | float | *required* | Superconductor layer thickness (nm) |
| `d_F` | float | *required* | Ferromagnet layer thickness (nm) |
| `xi_S` | float | *required* | Superconductor coherence length (nm) |
| `E_ex` | float | *required* | Exchange energy in ferromagnet (meV) |
| `n_grid` | int | `200` | Number of spatial grid points |
| `T` | float or None | `None` | Temperature (K). If `None`, defaults to $0.5 \cdot T_{c0}$ |

> **Note:** No `xi_F` parameter — the Eilenberger formulation operates in the
> clean limit where the mean free path exceeds the coherence length.

### Returns

| Type | Description |
|------|-------------|
| `x` : `numpy.ndarray` | Position array (nm), spanning $[-d_S, d_F]$ |
| `f` : `numpy.ndarray` | Anomalous Green's function $\|f(x)\|$ (Fermi-surface averaged) |

### Example

```python
import supermag

x, f = supermag.eilenberger.solve(
    Tc0=9.2, d_S=50.0, d_F=15.0,
    xi_S=38.0, E_ex=75.0, T=4.0
)
# f ≈ 1 deep in S, oscillatory decay into F
```
