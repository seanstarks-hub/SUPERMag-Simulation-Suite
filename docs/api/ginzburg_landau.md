# supermag.ginzburg_landau

Ginzburg-Landau free energy solver. Relaxes the TDGL equation on a 2D grid
using double-buffered snapshot Euler integration (no magnetic field).

**C++ engine:** `ginzburg_landau.cpp`

---

## `ginzburg_landau.minimize()`

Minimize the Ginzburg-Landau free energy on a 2D grid.

```python
supermag.ginzburg_landau.minimize(alpha, beta, kappa, nx, ny, dx)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `alpha` | float | *required* | GL $\alpha$ parameter ($< 0$ below $T_c$) |
| `beta` | float | *required* | GL $\beta$ parameter ($> 0$ for stability) |
| `kappa` | float | *required* | GL $\kappa = \lambda_L / \xi$. Type-II when $\kappa > 1/\sqrt{2}$ |
| `nx` | int | *required* | Grid width (number of columns) |
| `ny` | int | *required* | Grid height (number of rows) |
| `dx` | float | *required* | Grid spacing (nm) |

### Returns

| Type | Description |
|------|-------------|
| `numpy.ndarray` | Complex order parameter $\psi$ on 2D grid, shape `(ny, nx)` |

### Equation

EQ-11 — TDGL relaxation:

$$\frac{\partial \psi}{\partial t} = -\alpha\psi - \beta|\psi|^2\psi + \xi^2 \nabla^2\psi$$

Equilibrium: $|\psi|^2 = -\alpha/\beta$ (uniform, $\alpha < 0$). Simulation length scale $\xi = dx / \kappa$.

### Example

```python
import numpy as np
import supermag

psi = supermag.ginzburg_landau.minimize(
    alpha=-1.0, beta=1.0, kappa=5.0,
    nx=64, ny=64, dx=1.0
)
# |psi|^2 ≈ 1.0 everywhere for uniform ground state
print(f"Mean |ψ|² = {np.mean(np.abs(psi)**2):.3f}")
```
