# supermag.materials

Material parameter database. Ships with 3 superconductors and 6 ferromagnets
covering the most common S/F heterostructure combinations.

---

## `get_material()`

Retrieve material parameters by name.

```python
supermag.get_material(name)
```

### Parameters

| Name | Type | Description |
|------|------|-------------|
| `name` | str | Material name (case-sensitive) |

### Returns

| Type | Description |
|------|-------------|
| `dict` | Material parameters with units documented below |

### Raises

| Exception | Condition |
|-----------|-----------|
| `KeyError` | Material name not in database |

### Example

```python
import supermag

nb = supermag.get_material("Nb")
print(f"Nb Tc = {nb['Tc']} K, ξ_S = {nb['xi_S']} nm")
```

---

## `list_materials()`

List all available material names, grouped by type.

```python
supermag.list_materials()
```

### Returns

| Type | Description |
|------|-------------|
| `dict` | `{"superconductor": [...], "ferromagnet": [...]}` |

### Example

```python
import supermag

groups = supermag.list_materials()
print(groups["superconductor"])  # ['Al', 'Nb', 'Pb']
print(groups["ferromagnet"])     # ['Co', 'Cu0.43Ni0.57', 'CuNi', 'Fe', 'Ni', 'Py']
```

---

## Built-in Materials

### Superconductors

| Name | $T_c$ (K) | $\xi_S$ (nm) | $\lambda_L$ (nm) | $\Delta_0$ (meV) | $\rho$ (μΩ·cm) |
|------|-----------|--------------|------------------|------------------|----------------|
| Nb | 9.2 | 38.0 | 39.0 | 1.55 | 15.0 |
| Pb | 7.2 | 83.0 | 37.0 | 1.35 | 22.0 |
| Al | 1.2 | 1600.0 | 16.0 | 0.18 | 2.7 |

### Ferromagnets

| Name | $E_\text{ex}$ (meV) | $\xi_F$ (nm) | $D_F$ (m²/s) | $\rho$ (μΩ·cm) |
|------|---------------------|--------------|--------------|----------------|
| Fe | 256.0 | 0.7 | 2.5×10⁻⁴ | 10.0 |
| Co | 309.0 | 0.5 | 1.8×10⁻⁴ | 6.3 |
| Ni | 75.0 | 2.3 | 5.0×10⁻⁴ | 6.9 |
| Py (Permalloy) | 20.0 | 5.0 | 3.0×10⁻⁴ | 40.0 |
| CuNi | 5.0 | 10.0 | 4.0×10⁻⁴ | 35.0 |
| Cu₀.₄₃Ni₀.₅₇ | 11.2 | 4.2 | 4.0×10⁻⁴ | 50.0 |

### Dict keys

Superconductor entries contain: `type`, `Tc`, `xi_S`, `lambda_L`, `Delta_0`, `rho`.

Ferromagnet entries contain: `type`, `E_ex`, `xi_F`, `D_F`, `rho`.

### Custom materials

Register custom materials by adding to the `MATERIALS` dict:

```python
from supermag.materials import MATERIALS

MATERIALS["MgB2"] = {
    "type": "superconductor",
    "Tc": 39.0, "xi_S": 5.0, "lambda_L": 140.0, "Delta_0": 7.1,
    "rho": 5.0
}
```

---

## `get_interface()`

Retrieve interface coupling parameters for a superconductor/ferromagnet pair.

```python
supermag.get_interface(sc, fm, tier="sputtered")
```

### Parameters

| Name | Type | Description |
|------|------|-------------|
| `sc` | str | Superconductor name |
| `fm` | str | Ferromagnet name |
| `tier` | str | Fabrication quality tier: `"clean"`, `"sputtered"`, or `"oxidized"` |

### Returns

| Type | Description |
|------|-------------|
| `dict` | `{"gamma": float, "gamma_B": float, "source": str}` |

### Raises

| Exception | Condition |
|-----------|-----------|
| `KeyError` | (SC, FM) pair not in catalogue, or unknown tier |

### Example

```python
import supermag

iface = supermag.get_interface("Nb", "Fe", tier="sputtered")
print(f"γ = {iface['gamma']}, γ_B = {iface['gamma_B']}")
# γ = 0.3, γ_B = 0.05
```

---

## Interface Catalogue (`INTERFACES`)

Per-pair $(\gamma, \gamma_B)$ values. Three entries are validated against
published experimental data; the remaining 15 are estimated from the
validated pairs.

### Validated entries

| SC | FM | $\gamma$ | $\gamma_B$ | Reference |
|----|----|----------|-----------|-----------|
| Nb | Fe | 0.30 | 0.00 | Buzdin (1982) |
| Nb | Ni | 0.30 | 0.00 | Radovic (1991) |
| Nb | Cu₀.₄₃Ni₀.₅₇ | 0.15 | 0.30 | Fominov (2002) |

> $\gamma$ is a **phenomenological coupling parameter** — it cannot be
> reliably computed from the Kupriyanov–Lukichev mismatch ratio
> $\rho_S \xi_S / (\rho_F \xi_F)$. Validated values are obtained by
> fitting self-consistency models to experimental $T_c(d_F)$ curves.

---

## Fabrication Tiers (`FABRICATION_TIERS`)

Additive $\gamma_B$ corrections for interface quality.

| Tier | $\Delta\gamma_B$ | Description |
|------|------------------|-------------|
| `"clean"` | +0.00 | In-situ MBE, no vacuum break |
| `"sputtered"` | +0.05 | DC/RF sputtering, base < 10⁻⁷ Torr |
| `"oxidized"` | +0.30 | Ex-situ transfer with native oxide |

```python
from supermag.materials import FABRICATION_TIERS
FABRICATION_TIERS["sputtered"]["gamma_B_add"]  # 0.05
```
