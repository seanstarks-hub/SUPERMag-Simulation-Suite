# supermag.plotting

Visualization utilities for common SUPERMag outputs. All functions accept
an optional `ax` parameter for embedding in multi-panel figures.

---

## `plot_pair_amplitude()`

Plot pair amplitude $F(x)$ in the ferromagnet layer.

```python
supermag.plot_pair_amplitude(x, F, title=None, save_path=None, ax=None)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `x` | array_like | *required* | Position in ferromagnet (nm) |
| `F` | array_like | *required* | Pair amplitude (dimensionless) |
| `title` | str or None | `None` | Plot title |
| `save_path` | str or None | `None` | If given, save figure to this path |
| `ax` | matplotlib Axes or None | `None` | Axes to plot on. If `None`, creates new figure |

### Returns

| Type | Description |
|------|-------------|
| `(fig, ax)` | Matplotlib Figure and Axes |

### Example

```python
import supermag

x, F = supermag.pair_amplitude(d_F=10.0, xi_F=2.3)
supermag.plot_pair_amplitude(x, F, title="Nb/Ni F(x)")
```

---

## `plot_tc_vs_df()`

Plot critical temperature $T_c$ vs ferromagnet thickness $d_F$.

```python
supermag.plot_tc_vs_df(d_F, Tc, Tc0=None, title=None, save_path=None, ax=None)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `d_F` | array_like | *required* | Ferromagnet thickness (nm) |
| `Tc` | array_like | *required* | Critical temperature (K) |
| `Tc0` | float or None | `None` | Bulk $T_c$ for horizontal reference line |
| `title` | str or None | `None` | Plot title |
| `save_path` | str or None | `None` | If given, save figure to this path |
| `ax` | matplotlib Axes or None | `None` | Axes to plot on. If `None`, creates new figure |

### Returns

| Type | Description |
|------|-------------|
| `(fig, ax)` | Matplotlib Figure and Axes |

### Example

```python
import numpy as np
import supermag

d_F = np.linspace(0.5, 20.0, 200)
Tc = supermag.critical_temperature(
    Tc0=9.2, d_S=50.0, d_F_array=d_F,
    E_ex=256.0, xi_S=38.0, xi_F=0.7
)
supermag.plot_tc_vs_df(d_F, Tc, Tc0=9.2)
```
