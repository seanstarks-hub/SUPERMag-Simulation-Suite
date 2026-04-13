# supermag.sweeps

Sweep engines for systematic parameter studies. Wraps `critical_temperature()`
to compute $T_c$ over 1D or 2D parameter grids.

---

## `tc_parameter_sweep()`

Sweep a single parameter and compute $T_c$ for each value.

```python
supermag.tc_parameter_sweep(sweep_var, sweep_values, d_F_array=None, **fixed_params)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `sweep_var` | str | *required* | Parameter to sweep (see table below) |
| `sweep_values` | array_like | *required* | Values to iterate over |
| `d_F_array` | array_like or None | `None` | Ferromagnet thickness array (nm). **Required** when `sweep_var` is not `"d_F"` |
| `**fixed_params` | — | — | Forwarded to `critical_temperature()`. At minimum: `Tc0`, `d_S`, `E_ex`, `xi_S`, `xi_F` (unless swept) |

#### Supported sweep variables

| `sweep_var` | Description |
|-------------|-------------|
| `"d_F"` | Ferromagnet thickness (nm) |
| `"d_S"` | Superconductor thickness (nm) |
| `"gamma"` | Interface transparency |
| `"gamma_B"` | Interface barrier (Fominov) |
| `"E_ex"` | Exchange energy (meV) |
| `"xi_F"` | Ferromagnet coherence length (nm) |
| `"Tc0"` | Bulk critical temperature (K) |

### Returns

| Type | Description |
|------|-------------|
| `dict` | `{"sweep_var": str, "sweep_values": ndarray, "d_F": ndarray, "Tc": ndarray}` |

- When `sweep_var="d_F"`: `Tc` has shape `(len(sweep_values),)` — one $T_c$ per thickness.
- Otherwise: `Tc` has shape `(len(sweep_values), len(d_F_array))`.

### Auto $\xi_F$ recompute

When `sweep_var="E_ex"` and `xi_F` is **not** in `fixed_params`, $\xi_F$ is
automatically recomputed at each step via $\xi_F = \sqrt{\hbar D_F / E_\text{ex}}$.
Supply `D_F` (default 2.5×10⁻⁴ m²/s) for physical values.

### Example

```python
import numpy as np
import supermag

result = supermag.tc_parameter_sweep(
    "gamma", np.linspace(0.05, 0.5, 10),
    d_F_array=np.linspace(0.5, 20, 50),
    Tc0=9.2, d_S=50.0, E_ex=256.0, xi_S=38.0, xi_F=0.7
)
print(result["Tc"].shape)  # (10, 50)
```

---

## `tc_phase_diagram()`

Compute a 2D grid of $T_c$ values over two swept parameters.

```python
supermag.tc_phase_diagram(var1, vals1, var2, vals2, d_F_value=None, **fixed_params)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `var1` | str | *required* | First axis parameter name |
| `vals1` | array_like | *required* | Values along first axis |
| `var2` | str | *required* | Second axis parameter name |
| `vals2` | array_like | *required* | Values along second axis |
| `d_F_value` | float or None | `None` | Fixed $d_F$ (nm). **Required** when neither axis is `"d_F"` |
| `**fixed_params` | — | — | Forwarded to `critical_temperature()` |

### Returns

| Type | Description |
|------|-------------|
| `dict` | `{var1: ndarray, var2: ndarray, "Tc": ndarray}` where `Tc` has shape `(len(vals1), len(vals2))` |

### Example

```python
import numpy as np
import supermag

result = supermag.tc_phase_diagram(
    "gamma", np.linspace(0.05, 1.0, 30),
    "d_F", np.linspace(0.5, 15.0, 50),
    Tc0=9.2, d_S=50.0, E_ex=256.0, xi_S=38.0, xi_F=0.7
)
# result["Tc"].shape == (30, 50) — ready for pcolormesh
```
