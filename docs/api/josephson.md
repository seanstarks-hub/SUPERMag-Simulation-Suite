# supermag.josephson

Josephson junction solver for S/F/S heterostructures. Computes the
current-phase relation $I(\varphi)$ using a Matsubara frequency summation
that captures higher harmonics and the 0–π transition.

**C++ engine:** `josephson.cpp`

---

## `josephson.current_phase_relation()`

Compute the Josephson current-phase relation $I(\varphi)$ for an S/F/S junction.

```python
supermag.josephson.current_phase_relation(
    d_F, xi_F, E_ex, T, n_phases=100, Tc0=9.2, gamma_B=0.0
)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `d_F` | float | *required* | Ferromagnet thickness (nm) |
| `xi_F` | float | *required* | Ferromagnet coherence length (nm) |
| `E_ex` | float | *required* | Exchange energy (meV) |
| `T` | float | *required* | Temperature (K) |
| `n_phases` | int | `100` | Number of phase points from 0 to $2\pi$ |
| `Tc0` | float | `9.2` | Bulk superconductor critical temperature (K) |
| `gamma_B` | float | `0.0` | Interface barrier parameter (dimensionless). Attenuates the F-layer propagator as $1/(1 + \gamma_B |q_n| \xi_F)$ |

### Returns

| Type | Description |
|------|-------------|
| `phi` : `numpy.ndarray` | Phase difference array (rad), shape `(n_phases,)` |
| `I` : `numpy.ndarray` | Supercurrent (normalized to $\max\|I\| = 1$), shape `(n_phases,)` |

### Equation

EQ-9 — Matsubara frequency sum:

$$I(\varphi) = T \sum_n \text{Re}\!\left[ P_n \frac{\Delta^2 \sin\varphi}{\sqrt{(\omega_n^2 + \Delta^2\sin^2(\varphi/2))(\omega_n^2 + \Delta^2)}} \right]$$

where:
- $P_n = e^{-q_n d_F} e^{i\pi/4} / (1 + \gamma_B |q_n| \xi_F)$ — F-layer propagator with barrier damping
- $q_n = \sqrt{2(\omega_n/E_\text{ex} + i)} / \xi_F$ — complex wave vector
- $\omega_n = \pi k_B T (2n+1)$ — Matsubara frequency
- $\Delta(T) = 1.764\, k_B T_{c0} \sqrt{1 - T/T_{c0}}$
- Cutoff: $\omega_n > 20\Delta$ or $n > 500$

The $\sin^2(\varphi/2)$ denominator generates higher harmonics. The complex phase
of $P_n$ produces the 0–π oscillation as $d_F/\xi_F$ varies.

### Example

```python
import supermag

phi, I = supermag.josephson.current_phase_relation(
    d_F=3.0, xi_F=0.7, E_ex=256.0, T=4.0, gamma_B=0.3
)
# I > 0 at φ = π/2 → 0-junction
# I < 0 at φ = π/2 → π-junction
```
