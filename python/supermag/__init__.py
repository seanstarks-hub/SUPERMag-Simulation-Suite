"""
SUPERMag — Computational toolkit for superconductor/ferromagnet heterostructure research.

SUPERMag Lab, Texas State University.

Modules:
    proximity       — S/F bilayer proximity effect (pair amplitude, Tc vs d_F)
    usadel          — Usadel diffusive-limit solver
    eilenberger     — Eilenberger clean-limit solver
    bdg             — BdG tight-binding solver
    ginzburg_landau — Ginzburg-Landau free energy solver
    josephson       — Josephson junction CPR solver
    triplet         — Spin-triplet correlations solver
    materials       — Material parameter database
    plotting        — Visualization utilities
    themes          — Matplotlib figure styling presets
"""

__version__ = "0.1.0"

from supermag.proximity import pair_amplitude, critical_temperature
from supermag.materials import get_material, list_materials
from supermag.sweeps import tc_parameter_sweep, tc_phase_diagram
from supermag.themes import apply_theme, get_theme, list_themes, theme_context

# Solver modules — import submodules for convenient access
from supermag import bdg
from supermag import josephson
from supermag import usadel
from supermag import eilenberger
from supermag import ginzburg_landau
from supermag import triplet
