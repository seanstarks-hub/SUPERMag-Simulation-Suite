# supermag.bdg

BdG (Bogoliubov–de Gennes) tight-binding solver. Diagonalizes the Nambu
Hamiltonian on a 1D lattice to obtain the quasiparticle eigenvalue spectrum.

**C++ engine:** `bdg.cpp` (Jacobi eigensolver, $n \leq 2000$)

---

## `bdg.solve()`

Diagonalize the BdG Hamiltonian on a 1D tight-binding lattice.

```python
supermag.bdg.solve(n_sites, t_hop, Delta, E_ex, mu=0.0)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `n_sites` | int | *required* | Number of lattice sites |
| `t_hop` | float | *required* | Nearest-neighbor hopping energy (eV) |
| `Delta` | float | *required* | Superconducting pairing potential (meV) |
| `E_ex` | float | *required* | Exchange splitting in ferromagnet region (meV) |
| `mu` | float | `0.0` | Chemical potential (meV) |

### Returns

| Type | Description |
|------|-------------|
| `numpy.ndarray` | BdG eigenvalues (meV), sorted. Length $2 \times$ `n_sites` |

### Equation

EQ-10 — Nambu Hamiltonian:

$$H_\text{BdG} = \begin{pmatrix} -\mu I + E_\text{ex} I - t \cdot \text{tridiag} & \Delta I \\ \Delta^* I & +\mu I + E_\text{ex} I + t \cdot \text{tridiag} \end{pmatrix}$$

Eigenvalues are converted from eV to meV after diagonalization.

### Example

```python
import supermag

eigenvalues = supermag.bdg.solve(
    n_sites=100, t_hop=1.0, Delta=1.55, E_ex=256.0
)
# eigenvalues: sorted array of length 200
# Gap visible centered at E = 0, exchange-split bands
```
