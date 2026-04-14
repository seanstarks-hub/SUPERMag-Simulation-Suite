# supermag.proximity

S/F bilayer proximity effect solver. Computes pair amplitude $F(x)$ in the ferromagnet
and critical temperature $T_c(d_F)$ with oscillatory suppression from exchange splitting.

**C++ engine:** `critical_temp.cpp`, `pair_amplitude.cpp`

---

## `critical_temperature()`

Compute critical temperature $T_c$ as a function of ferromagnet thickness $d_F$.

```python
supermag.critical_temperature(
    Tc0, d_S, d_F_array, E_ex, xi_S, xi_F,
    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
    model="thin_s", phase="zero", depairing=None
)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `Tc0` | float | *required* | Bulk superconductor critical temperature (K) |
| `d_S` | float | *required* | Superconductor layer thickness (nm) |
| `d_F_array` | array_like | *required* | Ferromagnet thicknesses to evaluate (nm) |
| `E_ex` | float | *required* | Exchange energy in ferromagnet (meV) |
| `xi_S` | float | *required* | Superconductor coherence length (nm) |
| `xi_F` | float | *required* | Ferromagnet coherence length (nm) |
| `gamma` | float | `0.3` | Interface transparency parameter (dimensionless) |
| `gamma_B` | float | `0.0` | Interface barrier parameter (denominator of α, both models) |
| `D_F` | float | `2.5e-4` | Diffusion coefficient in ferromagnet (m²/s) |
| `model` | str | `"thin_s"` | Self-consistency model: `"thin_s"` (EQ-4) or `"fominov"` (EQ-5) |
| `phase` | str | `"zero"` | Junction phase: `"zero"` (coth kernel) or `"pi"` (tanh kernel) |
| `depairing` | dict or None | `None` | Pair-breaking channels (see below) |

#### Depairing channels

When `depairing` is a dict, the total depairing parameter $\lambda_\text{dep}$ (EQ-7) is the sum of:

| Key | Description |
|-----|-------------|
| `"ag"` | Abrikosov-Gor'kov magnetic impurity scattering |
| `"zeeman"` | Zeeman pair breaking |
| `"orbital"` | Orbital depairing |
| `"spin_orbit"` | Spin-orbit scattering |

All default to `0.0` if omitted.

### Returns

| Type | Description |
|------|-------------|
| `numpy.ndarray` | Critical temperature (K) for each $d_F$ value, shape matching `d_F_array` |

### Equations

- **`model="thin_s"`** uses EQ-4: $F(T) = \ln(T_{c0}/T) - \text{Re}[\psi(\tfrac{1}{2} + \alpha) - \psi(\tfrac{1}{2})]$ with $\alpha = \frac{\gamma}{\gamma_B + K} \cdot \frac{T_{c0}}{2\pi T} + \lambda_\text{dep}$
- **`model="fominov"`** uses EQ-5: same form with $\alpha = \frac{\gamma}{\gamma_B + K + \Omega_S(T)} \cdot \frac{T_{c0}}{2\pi T} + \lambda_\text{dep}$ where $\Omega_S(T) = \sqrt{T/T_{c0}} \coth(\sqrt{T/T_{c0}} \cdot d_S/\xi_S)$

### Example

```python
import numpy as np
import supermag

d_F = np.linspace(0.5, 20.0, 200)
Tc = supermag.critical_temperature(
    Tc0=9.2, d_S=50.0, d_F_array=d_F,
    E_ex=256.0, xi_S=38.0, xi_F=0.7, gamma=0.5
)
```

---

## `pair_amplitude()`

Compute the pair amplitude $F(x)$ in the ferromagnet layer of an S/F bilayer.

```python
supermag.pair_amplitude(d_F, xi_F, phase="zero", n_points=500)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `d_F` | float | *required* | Ferromagnet layer thickness (nm) |
| `xi_F` | float | *required* | Ferromagnet coherence length (nm). Typically 0.5–10 nm |
| `phase` | str | `"zero"` | `"zero"` (0-junction) or `"pi"` (π-junction) |
| `n_points` | int | `500` | Number of spatial grid points |

### Returns

| Type | Description |
|------|-------------|
| `x` : `numpy.ndarray` | Position array from 0 to $d_F$ (nm), shape `(n_points,)` |
| `F` : `numpy.ndarray` | Pair amplitude (dimensionless), shape `(n_points,)` |

### Equations

EQ-6:
- `phase="zero"`: $F(x) = e^{-x/\xi_F} \cos(x/\xi_F)$
- `phase="pi"`: $F(x) = e^{-x/\xi_F} \sin(x/\xi_F)$

### Example

```python
import supermag

x, F = supermag.pair_amplitude(d_F=10.0, xi_F=2.3, phase="zero")
print(f"F at interface: {F[0]:.3f}")  # 1.000
```
