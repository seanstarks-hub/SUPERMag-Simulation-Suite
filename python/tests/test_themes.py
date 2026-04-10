"""Tests for the SUPERMag themes system."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest
from supermag.themes import (
    apply_theme,
    get_theme,
    list_themes,
    register_theme,
    theme_context,
)


class TestListThemes:
    def test_returns_sorted_list(self):
        themes = list_themes()
        assert isinstance(themes, list)
        assert themes == sorted(themes)

    def test_all_builtin_present(self):
        themes = list_themes()
        for name in ("publication", "presentation", "draft", "dark"):
            assert name in themes


class TestGetTheme:
    def test_returns_dict(self):
        params = get_theme("publication")
        assert isinstance(params, dict)

    def test_returns_copy(self):
        p1 = get_theme("publication")
        p2 = get_theme("publication")
        p1["font.size"] = 999
        assert p2["font.size"] != 999

    def test_unknown_theme(self):
        with pytest.raises(KeyError, match="not found"):
            get_theme("nonexistent")

    def test_publication_serif(self):
        params = get_theme("publication")
        assert params["font.family"] == "serif"

    def test_publication_dpi(self):
        params = get_theme("publication")
        assert params["figure.dpi"] == 300

    def test_presentation_large_font(self):
        params = get_theme("presentation")
        assert params["font.size"] >= 14

    def test_draft_low_dpi(self):
        params = get_theme("draft")
        assert params["figure.dpi"] <= 100

    def test_dark_background(self):
        params = get_theme("dark")
        assert params["figure.facecolor"] == "#1e1e1e"


class TestApplyTheme:
    def setup_method(self):
        """Save original rcParams before each test."""
        self._original = mpl.rcParams.copy()

    def teardown_method(self):
        """Restore original rcParams after each test."""
        mpl.rcParams.update(self._original)

    def test_apply_modifies_rcparams(self):
        original_size = mpl.rcParams["font.size"]
        apply_theme("presentation")
        assert mpl.rcParams["font.size"] != original_size

    def test_apply_publication(self):
        apply_theme("publication")
        family = mpl.rcParams["font.family"]
        # font.family can be a string or a list
        if isinstance(family, list):
            assert "serif" in family
        else:
            assert family == "serif"


class TestThemeContext:
    def test_restores_on_exit(self):
        original_size = mpl.rcParams["font.size"]
        with theme_context("presentation"):
            assert mpl.rcParams["font.size"] >= 14
        assert mpl.rcParams["font.size"] == original_size

    def test_yields_params(self):
        with theme_context("draft") as params:
            assert isinstance(params, dict)
            assert "figure.dpi" in params


class TestRegisterTheme:
    def test_register_and_use(self):
        custom = {"font.size": 42, "lines.linewidth": 5.0}
        register_theme("my_custom", custom)
        assert "my_custom" in list_themes()
        params = get_theme("my_custom")
        assert params["font.size"] == 42

    def test_register_returns_copy(self):
        custom = {"font.size": 42}
        register_theme("temp_test", custom)
        custom["font.size"] = 99
        assert get_theme("temp_test")["font.size"] == 42
