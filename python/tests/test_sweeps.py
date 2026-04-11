"""Tests for supermag.sweeps — generic Tc parameter sweep utilities."""

import numpy as np
import pytest

from supermag.sweeps import tc_parameter_sweep, tc_phase_diagram


# Common parameters for Nb/Fe bilayer
NB_FE = dict(Tc0=9.2, d_S=50.0, E_ex=256.0, xi_S=38.0, xi_F=0.7)


class TestTcParameterSweep:
    """Tests for tc_parameter_sweep."""

    def test_sweep_gamma_shape(self):
        gammas = np.linspace(0.05, 0.5, 5)
        d_F = np.linspace(0.5, 10.0, 20)
        result = tc_parameter_sweep("gamma", gammas, d_F_array=d_F, **NB_FE)
        assert result["Tc"].shape == (5, 20)
        assert result["sweep_var"] == "gamma"
        assert len(result["sweep_values"]) == 5

    def test_sweep_d_F_shape(self):
        d_F_vals = np.linspace(1.0, 15.0, 10)
        result = tc_parameter_sweep("d_F", d_F_vals, **NB_FE)
        assert result["Tc"].shape == (10,)

    def test_sweep_d_F_bounded(self):
        d_F_vals = np.linspace(0.5, 20.0, 15)
        result = tc_parameter_sweep("d_F", d_F_vals, **NB_FE)
        assert np.all(result["Tc"] >= 0)
        assert np.all(result["Tc"] <= 9.2)

    def test_sweep_gamma_monotonic(self):
        """Larger gamma → stronger suppression at a given d_F."""
        gammas = np.array([0.05, 0.2, 0.5])
        d_F = np.array([5.0])
        result = tc_parameter_sweep("gamma", gammas, d_F_array=d_F, **NB_FE)
        Tc_at_5nm = result["Tc"][:, 0]
        # Tc should decrease (or stay flat) as gamma increases
        assert Tc_at_5nm[0] >= Tc_at_5nm[-1]

    def test_sweep_E_ex_auto_xi_F(self):
        """When sweeping E_ex without specifying xi_F, it auto-recomputes."""
        E_ex_vals = np.array([50.0, 150.0, 256.0])
        d_F = np.linspace(0.5, 10.0, 10)
        params = {k: v for k, v in NB_FE.items() if k != "xi_F"}
        result = tc_parameter_sweep(
            "E_ex", E_ex_vals, d_F_array=d_F,
            D_F=2.5e-4, **params,
        )
        assert result["Tc"].shape == (3, 10)

    def test_missing_d_F_array_raises(self):
        with pytest.raises(ValueError, match="d_F_array is required"):
            tc_parameter_sweep("gamma", [0.1, 0.2], **NB_FE)


class TestTcPhaseDiagram:
    """Tests for tc_phase_diagram."""

    def test_2d_shape(self):
        d_F = np.linspace(1, 10, 5)
        gammas = np.linspace(0.05, 0.3, 4)
        result = tc_phase_diagram(
            "d_F", d_F, "gamma", gammas, **NB_FE)
        assert result["Tc"].shape == (5, 4)

    def test_non_df_axes_need_d_F_value(self):
        with pytest.raises(ValueError, match="d_F_value is required"):
            tc_phase_diagram(
                "gamma", [0.1, 0.2], "gamma_B", [0.0, 0.3], **NB_FE)

    def test_non_df_axes_with_d_F_value(self):
        gammas = np.array([0.1, 0.3])
        gamma_Bs = np.array([0.0, 0.5])
        result = tc_phase_diagram(
            "gamma", gammas, "gamma_B", gamma_Bs,
            d_F_value=5.0, **NB_FE)
        assert result["Tc"].shape == (2, 2)
        assert np.all(result["Tc"] > 0)
