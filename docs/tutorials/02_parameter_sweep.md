# Running a Parameter Sweep

## Tc(d_F) for Multiple Ferromagnets

```python
import numpy as np
import supermag
from supermag.plotting import plot_tc_vs_df

nb = supermag.get_material("Nb")
d_F = np.linspace(0.5, 30.0, 200)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 6))

for fm_name in ["Fe", "Co", "Ni", "Py"]:
    fm = supermag.get_material(fm_name)
    Tc = supermag.critical_temperature(
        Tc0=nb["Tc"], d_S=50.0, d_F_array=d_F,
        E_ex=fm["E_ex"], xi_S=nb["xi_S"], xi_F=fm["xi_F"],
    )
    ax.plot(d_F, Tc, label=fm_name)

ax.axhline(y=nb["Tc"], ls="--", color="gray", label="Tc0")
ax.set_xlabel("d_F (nm)")
ax.set_ylabel("Tc (K)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig("sweep.png", dpi=150)
```
