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

| Name | $T_c$ (K) | $\xi_S$ (nm) | $\lambda_L$ (nm) | $\Delta_0$ (meV) |
|------|-----------|--------------|------------------|------------------|
| Nb | 9.2 | 38.0 | 39.0 | 1.55 |
| Pb | 7.2 | 83.0 | 37.0 | 1.35 |
| Al | 1.2 | 1600.0 | 16.0 | 0.18 |

### Ferromagnets

| Name | $E_\text{ex}$ (meV) | $\xi_F$ (nm) | $D_F$ (m²/s) |
|------|---------------------|--------------|--------------|
| Fe | 256.0 | 0.7 | 2.5×10⁻⁴ |
| Co | 309.0 | 0.5 | 1.8×10⁻⁴ |
| Ni | 75.0 | 2.3 | 5.0×10⁻⁴ |
| Py (Permalloy) | 20.0 | 5.0 | 3.0×10⁻⁴ |
| CuNi | 5.0 | 10.0 | 4.0×10⁻⁴ |
| Cu₀.₄₃Ni₀.₅₇ | 11.2 | 4.2 | 4.0×10⁻⁴ |

### Dict keys

Superconductor entries contain: `type`, `Tc`, `xi_S`, `lambda_L`, `Delta_0`.

Ferromagnet entries contain: `type`, `E_ex`, `xi_F`, `D_F`.

### Custom materials

Register custom materials by adding to the `MATERIALS` dict:

```python
from supermag.materials import MATERIALS

MATERIALS["MgB2"] = {
    "type": "superconductor",
    "Tc": 39.0, "xi_S": 5.0, "lambda_L": 140.0, "Delta_0": 7.1
}
```
