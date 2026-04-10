"""
Matplotlib visualization utilities for SUPERMag.

Publication-quality plots with labeled axes, units, grids, and legends.
"""

import numpy as np


def plot_pair_amplitude(x, F, title=None, save_path=None, ax=None):
    """
    Plot pair amplitude F(x) in the ferromagnet layer.

    Parameters
    ----------
    x : array_like
        Position in ferromagnet (nm).
    F : array_like
        Pair amplitude (dimensionless).
    title : str, optional
        Plot title.
    save_path : str, optional
        If given, save figure to this path.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    ax.plot(x, F, "b-", linewidth=1.5, label=r"$F(x)$")
    ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Position in F layer, $x$ (nm)")
    ax.set_ylabel("Pair amplitude $F(x)$ (dimensionless)")
    ax.set_title(title or "Pair Amplitude in Ferromagnet Layer")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig, ax


def plot_tc_vs_df(d_F, Tc, Tc0=None, title=None, save_path=None, ax=None):
    """
    Plot critical temperature Tc vs ferromagnet thickness d_F.

    Parameters
    ----------
    d_F : array_like
        Ferromagnet thickness (nm).
    Tc : array_like
        Critical temperature (K).
    Tc0 : float, optional
        Bulk Tc for reference line.
    title : str, optional
        Plot title.
    save_path : str, optional
        If given, save figure to this path.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    ax.plot(d_F, Tc, "r-o", markersize=3, linewidth=1.5, label=r"$T_c(d_F)$")
    if Tc0 is not None:
        ax.axhline(y=Tc0, color="blue", linewidth=1, linestyle="--",
                    label=f"$T_{{c0}} = {Tc0}$ K")
    ax.set_xlabel("Ferromagnet thickness $d_F$ (nm)")
    ax.set_ylabel("Critical temperature $T_c$ (K)")
    ax.set_title(title or r"$T_c$ vs. Ferromagnet Thickness")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig, ax
