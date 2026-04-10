# Quickstart

## Installation

```bash
pip install supermag
```

## Your First Calculation

```python
import numpy as np
import supermag

# Get material parameters for Niobium and Iron
nb = supermag.get_material("Nb")
fe = supermag.get_material("Fe")

# Compute Tc vs ferromagnet thickness
d_F = np.linspace(0.5, 20.0, 100)
Tc = supermag.critical_temperature(
    Tc0=nb["Tc"], d_S=50.0,
    d_F_array=d_F, E_ex=fe["E_ex"],
    xi_S=nb["xi_S"], xi_F=fe["xi_F"],
)

# Plot
from supermag.plotting import plot_tc_vs_df
plot_tc_vs_df(d_F, Tc, Tc0=nb["Tc"], save_path="tc_vs_df.png")
```
