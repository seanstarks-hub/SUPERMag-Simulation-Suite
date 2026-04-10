"""
SUPERMag — Computational toolkit for superconductor/ferromagnet heterostructure research.

SUPERMag Lab, Texas State University.

Modules:
    proximity      — S/F bilayer proximity effect (pair amplitude, Tc vs d_F)
    usadel         — Usadel diffusive-limit solver (stub)
    eilenberger    — Eilenberger clean-limit solver (stub)
    bdg            — BdG tight-binding solver (stub)
    ginzburg_landau — Ginzburg-Landau free energy solver (stub)
    josephson      — Josephson junction CPR solver (stub)
    triplet        — Spin-triplet correlations solver (stub)
    materials      — Material parameter database
    plotting       — Visualization utilities
    themes         — Matplotlib figure styling presets
"""

__version__ = "0.1.0"

from supermag.proximity import pair_amplitude, critical_temperature
from supermag.materials import get_material, list_materials
from supermag.themes import apply_theme, get_theme, list_themes, theme_context
