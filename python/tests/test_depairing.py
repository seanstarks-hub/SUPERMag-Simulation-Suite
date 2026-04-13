"""Tests for depairing channel computations and optimizer utilities."""

import numpy as np
import pytest
from supermag.depairing import (
    depairing_ag, depairing_zeeman,
    depairing_orbital_perp, depairing_orbital_par,
    depairing_soc, depairing_from_physical,
    optimize_tc, inverse_tc, fit_tc,
)


# ───────────────────── Individual Channels ─────────────────────

class TestDepairingAG:
    def test_positive(self):
        lam = depairing_ag(gamma_s_meV=0.1, T_kelvin=5.0)
        assert lam > 0

    def test_proportional_to_gamma_s(self):
        lam1 = depairing_ag(gamma_s_meV=0.1, T_kelvin=5.0)
        lam2 = depairing_ag(gamma_s_meV=0.2, T_kelvin=5.0)
        assert lam2 == pytest.approx(2.0 * lam1, rel=1e-10)

    def test_inversely_proportional_to_T(self):
        lam1 = depairing_ag(gamma_s_meV=0.1, T_kelvin=5.0)
        lam2 = depairing_ag(gamma_s_meV=0.1, T_kelvin=10.0)
        assert lam1 == pytest.approx(2.0 * lam2, rel=1e-10)

    def test_zero_rate(self):
        lam = depairing_ag(gamma_s_meV=0.0, T_kelvin=5.0)
        assert lam == pytest.approx(0.0)

    def test_invalid_T(self):
        with pytest.raises(ValueError, match="T_kelvin must be positive"):
            depairing_ag(gamma_s_meV=0.1, T_kelvin=0.0)


class TestDepairingZeeman:
    def test_positive(self):
        lam = depairing_zeeman(H_tesla=1.0, T_kelvin=5.0)
        assert lam > 0

    def test_quadratic_in_H(self):
        lam1 = depairing_zeeman(H_tesla=1.0, T_kelvin=5.0)
        lam2 = depairing_zeeman(H_tesla=2.0, T_kelvin=5.0)
        assert lam2 == pytest.approx(4.0 * lam1, rel=1e-10)

    def test_zero_field(self):
        lam = depairing_zeeman(H_tesla=0.0, T_kelvin=5.0)
        assert lam == pytest.approx(0.0)


class TestDepairingOrbital:
    def test_perp_positive(self):
        lam = depairing_orbital_perp(D_nm2ps=18.0, H_tesla=1.0,
                                      thickness_nm=50.0, T_kelvin=5.0)
        assert lam > 0

    def test_par_positive(self):
        lam = depairing_orbital_par(D_nm2ps=18.0, H_tesla=1.0,
                                     thickness_nm=50.0, T_kelvin=5.0)
        assert lam > 0

    def test_perp_larger_than_par(self):
        """Perpendicular field gives 4x stronger orbital depairing than parallel."""
        lam_perp = depairing_orbital_perp(D_nm2ps=18.0, H_tesla=1.0,
                                           thickness_nm=50.0, T_kelvin=5.0)
        lam_par = depairing_orbital_par(D_nm2ps=18.0, H_tesla=1.0,
                                         thickness_nm=50.0, T_kelvin=5.0)
        assert lam_perp == pytest.approx(4.0 * lam_par, rel=1e-10)

    def test_zero_field(self):
        lam = depairing_orbital_perp(D_nm2ps=18.0, H_tesla=0.0,
                                      thickness_nm=50.0, T_kelvin=5.0)
        assert lam == pytest.approx(0.0)


class TestDepairingSOC:
    def test_positive(self):
        lam = depairing_soc(Gamma_so_meV=0.05, T_kelvin=5.0)
        assert lam > 0

    def test_zero_rate(self):
        lam = depairing_soc(Gamma_so_meV=0.0, T_kelvin=5.0)
        assert lam == pytest.approx(0.0)


# ─────────────────── Composite Depairing ───────────────────────

class TestDepairingFromPhysical:
    def test_all_channels_nonneg(self):
        result = depairing_from_physical(
            gamma_s_meV=0.1, H_tesla=1.0, D_nm2ps=18.0,
            thickness_nm=50.0, Gamma_so_meV=0.05, T_kelvin=5.0)
        assert result["ag"] > 0
        assert result["zeeman"] > 0
        assert result["orbital"] > 0
        assert result["spin_orbit"] > 0

    def test_zero_field_channels(self):
        result = depairing_from_physical(
            gamma_s_meV=0.1, H_tesla=0.0, D_nm2ps=18.0,
            thickness_nm=50.0, Gamma_so_meV=0.05, T_kelvin=5.0)
        assert result["zeeman"] == pytest.approx(0.0)
        assert result["orbital"] == pytest.approx(0.0)
        assert result["ag"] > 0
        assert result["spin_orbit"] > 0

    def test_returns_dict(self):
        result = depairing_from_physical(
            gamma_s_meV=0.0, H_tesla=0.0, D_nm2ps=0.0,
            thickness_nm=50.0, Gamma_so_meV=0.0, T_kelvin=5.0)
        assert set(result.keys()) == {"ag", "zeeman", "orbital", "spin_orbit"}

    def test_invalid_T(self):
        with pytest.raises(ValueError, match="T_kelvin must be positive"):
            depairing_from_physical(
                gamma_s_meV=0.1, H_tesla=1.0, D_nm2ps=18.0,
                thickness_nm=50.0, Gamma_so_meV=0.05, T_kelvin=-1.0)


# ────────────────────── Optimizer ──────────────────────────────

class TestOptimizeTc:
    def test_returns_in_bounds(self):
        d_F = optimize_tc(
            Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            d_F_lo=0.5, d_F_hi=20.0, Tc_target=7.0)
        assert 0.5 <= d_F <= 20.0

    def test_invalid_target(self):
        with pytest.raises(ValueError, match="Tc_target must be positive"):
            optimize_tc(
                Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
                Tc_target=-1.0)

    def test_invalid_bounds(self):
        with pytest.raises(ValueError, match="d_F_lo must be < d_F_hi"):
            optimize_tc(
                Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
                d_F_lo=20.0, d_F_hi=1.0, Tc_target=7.0)


class TestInverseTc:
    def test_returns_in_bounds(self):
        d_F = inverse_tc(
            Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            Tc_target=7.0, d_F_lo=0.5, d_F_hi=20.0)
        assert 0.5 <= d_F <= 20.0

    def test_invalid_target(self):
        with pytest.raises(ValueError, match="Tc_target must be positive"):
            inverse_tc(
                Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
                Tc_target=-1.0, d_F_lo=0.5, d_F_hi=20.0)


class TestFitTc:
    def test_basic_fit(self):
        # Generate synthetic data with known gamma
        from supermag.proximity import critical_temperature
        d_F_data = np.linspace(1.0, 10.0, 8)
        Tc_data = critical_temperature(
            Tc0=9.2, d_S=50.0, d_F_array=d_F_data,
            E_ex=256.0, xi_S=38.0, xi_F=0.7, gamma=0.2)

        result = fit_tc(
            Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            gamma=0.1,  # wrong initial guess
            d_F_data=d_F_data, Tc_data=Tc_data,
            fit_gamma=True)

        assert "gamma" in result
        assert "chi2" in result
        assert result["chi2"] < 1.0  # should fit well

    def test_requires_data(self):
        with pytest.raises(ValueError, match="d_F_data and Tc_data are required"):
            fit_tc(Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0)

    def test_requires_min_points(self):
        with pytest.raises(ValueError, match="at least 2"):
            fit_tc(Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
                   d_F_data=np.array([1.0]), Tc_data=np.array([8.0]))

    def test_length_mismatch(self):
        with pytest.raises(ValueError, match="same length"):
            fit_tc(Tc0=9.2, d_S=50.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
                   d_F_data=np.array([1.0, 2.0]),
                   Tc_data=np.array([8.0]))
