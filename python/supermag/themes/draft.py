"""
Draft theme — Fast rendering for iterative development.

Rasterized elements, lower DPI, simplified styling for quick previews.
"""

PARAMS = {
    # Figure size: moderate
    "figure.figsize": (8, 5),
    "figure.dpi": 72,

    # Font: default sans-serif, moderate size
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,

    # Lines — slightly thicker for visibility at low DPI
    "lines.linewidth": 1.5,
    "lines.markersize": 5,

    # Axes
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "axes.facecolor": "#fafafa",

    # Grid
    "grid.alpha": 0.4,
    "grid.linestyle": "--",

    # Ticks
    "xtick.direction": "out",
    "ytick.direction": "out",

    # Legend
    "legend.frameon": True,

    # Enable rasterization for speed
    "agg.path.chunksize": 10000,

    # Saving — low DPI for drafts
    "savefig.dpi": 100,
    "savefig.bbox": "tight",
}
