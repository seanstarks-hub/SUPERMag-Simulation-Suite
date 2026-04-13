# supermag.depairing

Depairing channel computations and optimizer utilities. Computes
dimensionless pair-breaking parameters from physical laboratory inputs
and provides tools for matching experimental $T_c(d_F)$ data.

**C++ engine:** `depairing_models.cpp`, `optimizer.cpp`

---

## Individual Depairing Channels

### `depairing_ag()`

Abrikosov-Gor'kov pair-breaking from spin-flip scattering.

$$\lambda_\text{AG} = \Gamma_s / (2 k_B T)$$

```python
supermag.depairing_ag(gamma_s_meV, T_kelvin)
```

| Name | Type | Description |
|------|------|-------------|
| `gamma_s_meV` | float | Spin-flip scattering rate (meV) |
| `T_kelvin` | float | Temperature (K), must be > 0 |

**Returns:** `float` — dimensionless AG depairing parameter.

---

### `depairing_zeeman()`

Zeeman (Pauli paramagnetic) pair-breaking.

$$\lambda_Z = (\mu_B H)^2 / (2\pi k_B T)^2$$

```python
supermag.depairing_zeeman(H_tesla, T_kelvin)
```

| Name | Type | Description |
|------|------|-------------|
| `H_tesla` | float | Applied magnetic field (T) |
| `T_kelvin` | float | Temperature (K), must be > 0 |

**Returns:** `float` — dimensionless Zeeman depairing parameter.

---

### `depairing_orbital_perp()`

Orbital pair-breaking for perpendicular applied field.

$$\lambda_\text{orb\perp} = D (eH)^2 d^2 / (3\hbar^2 \cdot 2\pi k_B T)$$

```python
supermag.depairing_orbital_perp(D_nm2ps, H_tesla, thickness_nm, T_kelvin)
```

| Name | Type | Description |
|------|------|-------------|
| `D_nm2ps` | float | Diffusion coefficient (nm²/ps) |
| `H_tesla` | float | Applied magnetic field (T) |
| `thickness_nm` | float | Film thickness (nm) |
| `T_kelvin` | float | Temperature (K), must be > 0 |

**Returns:** `float` — dimensionless orbital depairing parameter (perp. field).

---

### `depairing_orbital_par()`

Orbital pair-breaking for parallel applied field.

$$\lambda_\text{orb\parallel} = D (eH)^2 d^2 / (12\hbar^2 \cdot 2\pi k_B T)$$

```python
supermag.depairing_orbital_par(D_nm2ps, H_tesla, thickness_nm, T_kelvin)
```

| Name | Type | Description |
|------|------|-------------|
| `D_nm2ps` | float | Diffusion coefficient (nm²/ps) |
| `H_tesla` | float | Applied magnetic field (T) |
| `thickness_nm` | float | Film thickness (nm) |
| `T_kelvin` | float | Temperature (K), must be > 0 |

**Returns:** `float` — dimensionless orbital depairing parameter (parallel field).

---

### `depairing_soc()`

Spin-orbit coupling pair-breaking.

$$\lambda_\text{SO} = \Gamma_\text{so} / (2 k_B T)$$

```python
supermag.depairing_soc(Gamma_so_meV, T_kelvin)
```

| Name | Type | Description |
|------|------|-------------|
| `Gamma_so_meV` | float | Spin-orbit scattering rate (meV) |
| `T_kelvin` | float | Temperature (K), must be > 0 |

**Returns:** `float` — dimensionless spin-orbit depairing parameter.

---

### `depairing_from_physical()`

Compute all depairing channels from physical laboratory inputs.

```python
supermag.depairing_from_physical(
    gamma_s_meV, H_tesla, D_nm2ps,
    thickness_nm, Gamma_so_meV, T_kelvin
)
```

**Returns:** `dict` — `{"ag": float, "zeeman": float, "orbital": float, "spin_orbit": float}`

---

## Optimizer Utilities

### `optimize_tc()`

Find the F-layer thickness that produces a target $T_c$ using golden-section search (EQ-20).

```python
supermag.optimize_tc(
    Tc0, d_S, xi_S, xi_F, E_ex,
    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
    model="thin_s", phase="zero",
    d_F_lo=0.5, d_F_hi=50.0, Tc_target=4.2,
    depairing=None
)
```

**Returns:** `float` — optimal $d_F$ (nm).

---

### `inverse_tc()`

Find $d_F$ that produces exactly $T_c = T_\text{target}$ via Brent's method (EQ-21).

```python
supermag.inverse_tc(
    Tc0, d_S, xi_S, xi_F, E_ex,
    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
    model="thin_s", phase="zero",
    Tc_target=4.2, d_F_lo=0.5, d_F_hi=50.0,
    depairing=None
)
```

**Returns:** `float` — $d_F$ (nm) giving $T_c = T_\text{target}$, or NaN if no root in bracket.

---

### `fit_tc()`

Fit proximity parameters to experimental $T_c(d_F)$ data using Nelder-Mead minimization (EQ-22).

```python
supermag.fit_tc(
    Tc0, d_S, xi_S, xi_F, E_ex,
    gamma=0.3, gamma_B=0.0, D_F=2.5e-4,
    model="thin_s", phase="zero",
    d_F_data=..., Tc_data=...,
    fit_gamma=True, fit_gamma_B=False,
    fit_E_ex=False, fit_xi_F=False,
    depairing=None
)
```

**Returns:** `dict` — `{"gamma": float, "gamma_B": float, "E_ex": float, "xi_F": float, "chi2": float}`

### Example

```python
import numpy as np
import supermag

# Experimental data
d_F_exp = np.array([1.0, 3.0, 5.0, 8.0, 12.0, 15.0])
Tc_exp = np.array([8.5, 6.2, 4.8, 5.1, 4.5, 4.3])

result = supermag.fit_tc(
    Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
    d_F_data=d_F_exp, Tc_data=Tc_exp,
    fit_gamma=True, fit_E_ex=True,
)
print(f"Fitted γ = {result['gamma']:.3f}, E_ex = {result['E_ex']:.1f} meV")
print(f"χ² = {result['chi2']:.4f}")
```
