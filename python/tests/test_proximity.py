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

    def test_boundary(self):
        x, F = pair_amplitude(d_F=10.0, xi_F=2.0, F0=1.0)
        assert F[0] == pytest.approx(1.0, abs=1e-12)

    def test_decay(self):
        x, F = pair_amplitude(d_F=20.0, xi_F=2.0)
        # F should decay overall
        assert abs(F[-1]) < abs(F[0])

    def test_oscillation(self):
        x, F = pair_amplitude(d_F=30.0, xi_F=2.0, n_points=1000)
        # Should have sign changes (oscillation)
        sign_changes = np.sum(np.diff(np.sign(F)) != 0)
        assert sign_changes >= 2  # at least a couple oscillations

    def test_f0_scaling(self):
        _, F1 = pair_amplitude(d_F=10.0, xi_F=2.0, F0=1.0)
        _, F2 = pair_amplitude(d_F=10.0, xi_F=2.0, F0=2.0)
        np.testing.assert_allclose(F2, 2.0 * F1, rtol=1e-12)


class TestCriticalTemperature:
    def test_basic(self, nb_fe_params, d_F_array):
        Tc = critical_temperature(d_F_array=d_F_array, **nb_fe_params)
        assert len(Tc) == len(d_F_array)
        assert np.all(Tc >= 0)
        assert np.all(Tc <= nb_fe_params["Tc0"] + 1e-10)

    def test_thin_limit(self, nb_fe_params):
        """Very thin F layer should have Tc close to Tc0."""
        d_F = np.array([0.01])
        Tc = critical_temperature(d_F_array=d_F, **nb_fe_params)
        assert Tc[0] > 0.8 * nb_fe_params["Tc0"]

    def test_array_output(self, nb_fe_params):
        d_F = np.linspace(1.0, 10.0, 20)
        Tc = critical_temperature(d_F_array=d_F, **nb_fe_params)
        assert Tc.shape == d_F.shape
