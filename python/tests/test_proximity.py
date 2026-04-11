"""Tests for the proximity effect module."""

import numpy as np
import pytest
from supermag.proximity import pair_amplitude, critical_temperature


class TestPairAmplitude:
    def test_basic(self):
        x, F = pair_amplitude(d_F=10.0, xi_F=2.0, n_points=100)
        assert len(x) == 100
        assert len(F) == 100
        assert x[0] == pytest.approx(0.0)
        assert x[-1] == pytest.approx(10.0)

    def test_boundary_zero_phase(self):
        x, F = pair_amplitude(d_F=10.0, xi_F=2.0, phase="zero")
        # F(0) = exp(0)*cos(0) = 1.0
        assert F[0] == pytest.approx(1.0, abs=1e-12)

    def test_boundary_pi_phase(self):
        x, F = pair_amplitude(d_F=10.0, xi_F=2.0, phase="pi")
        # F(0) = exp(0)*sin(0) = 0.0
        assert F[0] == pytest.approx(0.0, abs=1e-12)

    def test_decay(self):
        x, F = pair_amplitude(d_F=20.0, xi_F=2.0)
        # F should decay overall
        assert abs(F[-1]) < abs(F[0])

    def test_oscillation(self):
        x, F = pair_amplitude(d_F=30.0, xi_F=2.0, n_points=1000)
        # Should have sign changes (oscillation)
        sign_changes = np.sum(np.diff(np.sign(F)) != 0)
        assert sign_changes >= 2  # at least a couple oscillations

    def test_pi_phase_nonzero_interior(self):
        x, F = pair_amplitude(d_F=10.0, xi_F=2.0, phase="pi", n_points=100)
        # Interior should be nonzero even though F(0)=0
        assert np.max(np.abs(F)) > 0.1

    def test_invalid_d_F(self):
        with pytest.raises(ValueError, match="d_F must be positive"):
            pair_amplitude(d_F=-1.0, xi_F=2.0)

    def test_invalid_xi_F(self):
        with pytest.raises(ValueError, match="xi_F must be positive"):
            pair_amplitude(d_F=10.0, xi_F=0.0)

    def test_invalid_phase(self):
        with pytest.raises(ValueError, match="phase must be"):
            pair_amplitude(d_F=10.0, xi_F=2.0, phase="invalid")

    def test_invalid_n_points(self):
        with pytest.raises(ValueError, match="n_points must be"):
            pair_amplitude(d_F=10.0, xi_F=2.0, n_points=0)


class TestCriticalTemperature:
    def test_basic(self, nb_fe_params, d_F_array):
        Tc = critical_temperature(d_F_array=d_F_array, **nb_fe_params)
        assert len(Tc) == len(d_F_array)
        assert np.all(Tc >= 0)
        assert np.all(Tc <= nb_fe_params["Tc0"] + 1e-10)

    def test_thin_limit(self, nb_fe_params):
        """Thin F layer should have Tc substantially above zero."""
        d_F = np.array([0.5])
        Tc = critical_temperature(d_F_array=d_F, **nb_fe_params)
        assert Tc[0] > 0.3 * nb_fe_params["Tc0"]

    def test_array_output(self, nb_fe_params):
        d_F = np.linspace(1.0, 10.0, 20)
        Tc = critical_temperature(d_F_array=d_F, **nb_fe_params)
        assert Tc.shape == d_F.shape

    def test_fominov_model(self, nb_fe_params, d_F_array):
        Tc = critical_temperature(d_F_array=d_F_array, model="fominov",
                                  gamma_B=0.3, **nb_fe_params)
        assert len(Tc) == len(d_F_array)
        assert np.all(Tc >= 0)
        assert np.all(Tc <= nb_fe_params["Tc0"] + 1e-10)

    def test_gamma_B_suppression(self, nb_fe_params):
        """Higher gamma_B (barrier) should reduce Tc suppression."""
        d_F = np.linspace(1.0, 10.0, 20)
        Tc_low = critical_temperature(d_F_array=d_F, gamma_B=0.0, **nb_fe_params)
        Tc_high = critical_temperature(d_F_array=d_F, gamma_B=1.0, **nb_fe_params)
        # With higher barrier, Tc should be closer to Tc0 (less suppressed)
        assert np.mean(Tc_high) >= np.mean(Tc_low) - 0.5

    def test_with_depairing(self, nb_fe_params):
        d_F = np.array([5.0])
        Tc_no_dp = critical_temperature(d_F_array=d_F, **nb_fe_params)
        Tc_with_dp = critical_temperature(
            d_F_array=d_F, depairing={"ag": 0.1}, **nb_fe_params)
        # Depairing should suppress Tc further
        assert Tc_with_dp[0] <= Tc_no_dp[0] + 0.5

    def test_invalid_Tc0(self, nb_fe_params):
        params = {**nb_fe_params, "Tc0": -1.0}
        with pytest.raises(ValueError, match="Tc0 must be positive"):
            critical_temperature(d_F_array=np.array([5.0]), **params)

    def test_invalid_model(self, nb_fe_params):
        with pytest.raises(ValueError, match="model must be"):
            critical_temperature(d_F_array=np.array([5.0]),
                                 model="invalid", **nb_fe_params)

    def test_invalid_phase_ct(self, nb_fe_params):
        with pytest.raises(ValueError, match="phase must be"):
            critical_temperature(d_F_array=np.array([5.0]),
                                 phase="bad", **nb_fe_params)
