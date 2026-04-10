"""
Dark theme — Dark background for IDE notebooks and screen presentations.

Dark axes/figure background with light-colored text and bright plot lines.
"""

PARAMS = {
    # Figure
    "figure.figsize": (8, 5),
    "figure.dpi": 100,
    "figure.facecolor": "#1e1e1e",
    "figure.edgecolor": "#1e1e1e",

    # Font
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,

    # Colors — light text on dark background
    "text.color": "#d4d4d4",
    "axes.facecolor": "#252526",
    "axes.edgecolor": "#555555",
    "axes.labelcolor": "#d4d4d4",
    "xtick.color": "#d4d4d4",
    "ytick.color": "#d4d4d4",

    # Bright line colors for contrast on dark
    "axes.prop_cycle": __import__("cycler").cycler(
        color=["#569cd6", "#ce9178", "#6a9955", "#dcdcaa",
               "#c586c0", "#4ec9b0", "#d7ba7d", "#9cdcfe"]),

    # Lines
    "lines.linewidth": 1.8,
    "lines.markersize": 5,

    # Axes
    "axes.linewidth": 0.8,
    "axes.grid": True,

    # Grid — subtle on dark
    "grid.color": "#404040",
    "grid.alpha": 0.5,
    "grid.linewidth": 0.5,

    # Ticks
    "xtick.direction": "out",
    "ytick.direction": "out",

    # Legend
    "legend.frameon": True,
    "legend.facecolor": "#333333",
    "legend.edgecolor": "#555555",
    "legend.framealpha": 0.9,

    # Saving
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.facecolor": "#1e1e1e",
}
