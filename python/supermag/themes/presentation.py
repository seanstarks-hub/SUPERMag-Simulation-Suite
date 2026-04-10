"""
Presentation theme — Large fonts, bold colors for slides and talks.

16:9 aspect ratio, large labels, thick lines, high-contrast palette.
"""

PARAMS = {
    # Figure size: 16:9 slide-friendly
    "figure.figsize": (10, 5.625),
    "figure.dpi": 150,

    # Font: sans-serif for readability on screen
    "font.family": "sans-serif",
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 20,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,

    # Lines and markers
    "lines.linewidth": 2.5,
    "lines.markersize": 8,

    # Axes
    "axes.linewidth": 1.5,
    "axes.grid": True,
    "axes.prop_cycle": __import__("cycler").cycler(
        color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
               "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]),

    # Grid
    "grid.alpha": 0.3,
    "grid.linewidth": 0.8,

    # Ticks
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,

    # Legend
    "legend.frameon": True,
    "legend.framealpha": 0.8,
    "legend.edgecolor": "0.8",

    # Saving
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
}
