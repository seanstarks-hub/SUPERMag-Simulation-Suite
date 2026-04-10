"""
SUPERMag Plot Themes — Matplotlib rcParams presets for consistent styling.

Built-in themes:
    publication  — APS/PRB-ready: serif fonts, 3.375" single-column width
    presentation — Large fonts, bold colors for talks
    draft        — Speed-optimized with rasterized elements
    dark         — Dark background for IDE/notebook use

Usage:
    >>> from supermag.themes import apply_theme
    >>> apply_theme("publication")

    Or use as a context manager:
    >>> with theme_context("presentation"):
    ...     fig, ax = plt.subplots()
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from contextlib import contextmanager

from supermag.themes.publication import PARAMS as _pub
from supermag.themes.presentation import PARAMS as _pres
from supermag.themes.draft import PARAMS as _draft
from supermag.themes.dark import PARAMS as _dark

_REGISTRY = {
    "publication": _pub,
    "presentation": _pres,
    "draft": _draft,
    "dark": _dark,
}


def list_themes():
    """Return sorted list of available theme names."""
    return sorted(_REGISTRY.keys())


def get_theme(name):
    """Return a copy of the rcParams dict for the named theme.

    Parameters
    ----------
    name : str
        Theme name (one of ``list_themes()``).

    Returns
    -------
    dict
        Matplotlib rcParams overrides.

    Raises
    ------
    KeyError
        If theme name is not registered.
    """
    if name not in _REGISTRY:
        available = ", ".join(sorted(_REGISTRY.keys()))
        raise KeyError(f"Theme '{name}' not found. Available: {available}")
    return dict(_REGISTRY[name])


def apply_theme(name):
    """Apply a named theme globally (modifies ``matplotlib.rcParams`` in place).

    Parameters
    ----------
    name : str
        Theme name.
    """
    params = get_theme(name)
    mpl.rcParams.update(params)


@contextmanager
def theme_context(name):
    """Context manager that temporarily applies a theme, then restores defaults.

    Parameters
    ----------
    name : str
        Theme name.

    Yields
    ------
    dict
        The theme's rcParams dict.

    Examples
    --------
    >>> with theme_context("publication"):
    ...     fig, ax = plt.subplots()
    ...     ax.plot([0, 1], [0, 1])
    """
    params = get_theme(name)
    with mpl.rc_context(rc=params):
        yield params


def register_theme(name, params):
    """Register a custom theme.

    Parameters
    ----------
    name : str
        Theme name (will overwrite existing if same name).
    params : dict
        Matplotlib rcParams overrides.
    """
    _REGISTRY[name] = dict(params)
