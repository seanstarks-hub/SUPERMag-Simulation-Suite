"""
Publication theme — APS/PRB journal-ready figure styling.

Single-column width: 3.375 inches (8.6 cm)
Double-column width: 6.75 inches (17.1 cm)
Serif font (Times-like), tight tick marks, thin lines.
"""

PARAMS = {
    # Figure size: single-column APS format
    "figure.figsize": (3.375, 2.53),
    "figure.dpi": 300,

    # Font: serif for journal articles
    "font.family": "serif",
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,

    # Lines and markers
    "lines.linewidth": 1.0,
    "lines.markersize": 3,

    # Axes
    "axes.linewidth": 0.6,
    "axes.grid": False,
    "axes.prop_cycle": __import__("cycler").cycler(
        color=["#000000", "#e41a1c", "#377eb8", "#4daf4a",
               "#984ea3", "#ff7f00", "#a65628", "#f781bf"]),

    # Ticks
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.top": True,
    "ytick.right": True,

    # Legend
    "legend.frameon": False,
    "legend.borderpad": 0.3,

    # Saving
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}
