# Defining Custom Materials

## Using Custom Parameters

You don't need to use the built-in database. Pass parameters directly:

```python
import numpy as np
import supermag

# Custom weak ferromagnet parameters
d_F = np.linspace(1, 50, 200)
Tc = supermag.critical_temperature(
    Tc0=9.2,        # K (Nb)
    d_S=30.0,       # nm
    d_F_array=d_F,
    E_ex=10.0,      # meV (weak ferromagnet)
    xi_S=38.0,      # nm
    xi_F=8.0,       # nm (long coherence length)
)
```

## Adding to the Database

You can extend the materials database at runtime:

```python
from supermag.materials import MATERIALS

MATERIALS["MyAlloy"] = {
    "type": "ferromagnet",
    "E_ex": 15.0,    # meV
    "xi_F": 6.0,     # nm
    "D_F": 3.5e-4,   # m^2/s
}
```
