"""Tests for plotting utilities."""

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")

from supermag.plotting import plot_pair_amplitude, plot_tc_vs_df


class TestPlotPairAmplitude:
    def test_returns_fig_ax(self):
        x = np.linspace(0, 10, 50)
        F = np.exp(-x) * np.cos(x)
        fig, ax = plot_pair_amplitude(x, F)
        assert fig is not None
        assert ax is not None

    def test_custom_title(self):
        x = np.linspace(0, 10, 50)
        F = np.exp(-x) * np.cos(x)
        fig, ax = plot_pair_amplitude(x, F, title="Test Title")
        assert ax.get_title() == "Test Title"

    def test_existing_axes(self):
        import matplotlib.pyplot as plt
        fig0, ax0 = plt.subplots()
        x = np.linspace(0, 10, 50)
        F = np.exp(-x) * np.cos(x)
        fig, ax = plot_pair_amplitude(x, F, ax=ax0)
        assert ax is ax0

    def test_save_to_file(self, tmp_path):
        x = np.linspace(0, 10, 50)
        F = np.exp(-x) * np.cos(x)
        path = str(tmp_path / "test.png")
        fig, ax = plot_pair_amplitude(x, F, save_path=path)
        import os
        assert os.path.exists(path)


class TestPlotTcVsDf:
    def test_returns_fig_ax(self):
        d_F = np.linspace(0.5, 20, 50)
        Tc = 9.2 * np.exp(-d_F / 5)
        fig, ax = plot_tc_vs_df(d_F, Tc)
        assert fig is not None
        assert ax is not None

    def test_tc0_reference_line(self):
        d_F = np.linspace(0.5, 20, 50)
        Tc = 9.2 * np.exp(-d_F / 5)
        fig, ax = plot_tc_vs_df(d_F, Tc, Tc0=9.2)
        # Should have at least 2 artists (data line + reference line)
        assert len(ax.lines) >= 2

    def test_existing_axes(self):
        import matplotlib.pyplot as plt
        fig0, ax0 = plt.subplots()
        d_F = np.linspace(0.5, 20, 50)
        Tc = 9.2 * np.exp(-d_F / 5)
        fig, ax = plot_tc_vs_df(d_F, Tc, ax=ax0)
        assert ax is ax0

    def test_save_to_file(self, tmp_path):
        d_F = np.linspace(0.5, 20, 50)
        Tc = 9.2 * np.exp(-d_F / 5)
        path = str(tmp_path / "tc_plot.png")
        fig, ax = plot_tc_vs_df(d_F, Tc, save_path=path)
        import os
        assert os.path.exists(path)

    def test_custom_title(self):
        d_F = np.linspace(0.5, 20, 50)
        Tc = 9.2 * np.exp(-d_F / 5)
        fig, ax = plot_tc_vs_df(d_F, Tc, title="Custom")
        assert ax.get_title() == "Custom"
