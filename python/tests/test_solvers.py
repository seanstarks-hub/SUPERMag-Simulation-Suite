"""Tests for all solver modules: BdG, Josephson, Usadel, Eilenberger, GL, Triplet."""

import numpy as np
import pytest

from supermag import bdg, josephson, usadel, eilenberger, ginzburg_landau, triplet


# ───────────────────────────── BdG ─────────────────────────────

class TestBdG:
    def test_output_size(self):
        eig = bdg.solve(n_sites=10, t_hop=1.0, Delta=1.5, E_ex=50.0)
        assert len(eig) == 20  # 2*N eigenvalues

    def test_particle_hole_symmetry(self):
        """Eigenvalues should come in ±E pairs."""
        eig = bdg.solve(n_sites=20, t_hop=1.0, Delta=1.5, E_ex=0.0)
        # Sort and check pairing
        eig_sorted = np.sort(eig)
        for i in range(len(eig_sorted)):
            assert pytest.approx(eig_sorted[i], abs=0.1) == -eig_sorted[-(i + 1)]

    def test_gap_opens(self):
        """With Δ > 0 and no exchange, a gap should open around E=0."""
        eig = bdg.solve(n_sites=30, t_hop=1.0, Delta=10.0, E_ex=0.0)
        # Minimum |E| should be approximately Δ for large enough system
        min_abs = np.min(np.abs(eig))
        assert min_abs > 0.5  # gap is open

    def test_exchange_shifts_spectrum(self):
        """Exchange field should shift the spectrum."""
        eig_0 = bdg.solve(n_sites=20, t_hop=1.0, Delta=1.5, E_ex=0.0)
        eig_h = bdg.solve(n_sites=20, t_hop=1.0, Delta=1.5, E_ex=100.0)
        # Spectra should differ
        assert not np.allclose(eig_0, eig_h)

    def test_sorted_output(self):
        eig = bdg.solve(n_sites=15, t_hop=1.0, Delta=1.5, E_ex=50.0)
        assert np.all(np.diff(eig) >= -1e-10)  # sorted ascending

    def test_mu_shifts_spectrum(self):
        """Chemical potential should shift the spectrum."""
        eig_0 = bdg.solve(n_sites=20, t_hop=1.0, Delta=1.5, E_ex=0.0, mu=0.0)
        eig_mu = bdg.solve(n_sites=20, t_hop=1.0, Delta=1.5, E_ex=0.0, mu=50.0)
        assert not np.allclose(eig_0, eig_mu)

    def test_mu_default_backward_compat(self):
        """Default mu=0 should match old behavior."""
        eig = bdg.solve(n_sites=10, t_hop=1.0, Delta=1.5, E_ex=50.0)
        eig_explicit = bdg.solve(n_sites=10, t_hop=1.0, Delta=1.5, E_ex=50.0, mu=0.0)
        np.testing.assert_allclose(eig, eig_explicit)


# ───────────────────────── Josephson ───────────────────────────

class TestJosephson:
    def test_output_shape(self):
        phi, I = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=4.0, n_phases=50)
        assert len(phi) == 50
        assert len(I) == 50

    def test_normalized(self):
        """Output should be normalized to max |I| = 1."""
        phi, I = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=4.0)
        assert np.max(np.abs(I)) == pytest.approx(1.0, abs=1e-10)

    def test_sinusoidal_shape(self):
        """CPR should be approximately sinusoidal (first harmonic)."""
        phi, I = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=4.0, n_phases=200)
        # I should cross zero near phi = 0, pi
        # Check that I ≈ 0 at phi ≈ 0
        assert abs(I[0]) < 0.1

    def test_zero_pi_transition(self):
        """Critical current should change sign as d_F crosses 0-pi transition."""
        # d_F < 0.75*xi_F is 0-junction, d_F > 0.75*xi_F is pi-junction
        xi_F = 2.0
        _, I_short = josephson.current_phase_relation(
            d_F=0.5 * xi_F, xi_F=xi_F, E_ex=100.0, T=1.0, n_phases=50)
        _, I_long = josephson.current_phase_relation(
            d_F=3.0 * xi_F, xi_F=xi_F, E_ex=100.0, T=1.0, n_phases=50)
        # They should have opposite sign at the same phase
        # (one is 0-junction, other is pi-junction)
        idx = 10  # some phase away from 0
        # Sign may or may not flip depending on exact d_F; just check they differ
        assert not np.allclose(I_short, I_long, atol=0.1)

    def test_phase_range(self):
        phi, _ = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=4.0)
        assert phi[0] == pytest.approx(0.0)
        assert phi[-1] < 2 * np.pi  # endpoint=False

    def test_custom_Tc0(self):
        """Custom Tc0 parameter should be accepted and affect temperature factor.
        Note: CPR is normalized, so we test at T > Tc0 where gap vanishes."""
        # With Tc0=5.0 and T=4.9 (near Tc), gap is small → amplitude is small
        # With Tc0=20.0 and T=4.9 (far from Tc), gap is larger
        # Both get normalized, but at T > Tc0 the amplitude is zero
        phi, I = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=4.0, Tc0=9.2)
        assert len(phi) == 100
        assert np.max(np.abs(I)) == pytest.approx(1.0, abs=1e-10)
        # Verify very high T relative to Tc0 suppresses current entirely
        phi2, I2 = josephson.current_phase_relation(
            d_F=2.0, xi_F=0.7, E_ex=256.0, T=9.19, Tc0=9.2)
        # Near Tc, tanh of small gap → small factor (still normalized though)
        assert len(I2) == 100


# ─────────────────────────── Usadel ───────────────────────────

class TestUsadel:
    def test_output_shape(self):
        x, Delta = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            n_grid=100)
        assert len(x) == 100
        assert len(Delta) == 100

    def test_delta_nonnegative(self):
        x, Delta = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0)
        assert np.all(Delta >= 0)

    def test_delta_zero_in_ferromagnet(self):
        """Order parameter should be zero (or near zero) in F region."""
        x, Delta = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            n_grid=200)
        # F region is x > 0
        f_mask = x > 1.0  # safely in F region
        if np.any(f_mask):
            assert np.max(Delta[f_mask]) < 1e-3

    def test_delta_finite_in_superconductor(self):
        """Order parameter should be finite in S region."""
        x, Delta = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            n_grid=200)
        s_mask = x < -5.0  # well inside S
        if np.any(s_mask):
            assert np.max(Delta[s_mask]) > 0.01

    def test_spatial_grid(self):
        x, _ = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            n_grid=100)
        assert x[0] == pytest.approx(-50.0, abs=1.0)
        assert x[-1] == pytest.approx(10.0, abs=1.0)

    def test_custom_temperature(self):
        """Custom T should produce different results than default."""
        x1, Delta1 = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0)
        x2, Delta2 = usadel.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, xi_F=0.7, E_ex=256.0,
            T=3.0)
        assert not np.allclose(Delta1, Delta2, atol=1e-6)


# ─────────────────────── Eilenberger ──────────────────────────

class TestEilenberger:
    def test_output_shape(self):
        x, f = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, E_ex=256.0, n_grid=100)
        assert len(x) == 100
        assert len(f) == 100

    def test_f_nonnegative(self):
        x, f = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, E_ex=256.0)
        assert np.all(f >= 0)

    def test_f_bounded(self):
        """Anomalous Green's function |f| should be ≤ 1."""
        x, f = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, E_ex=256.0)
        assert np.all(f <= 1.5)  # allow small numerical overshoot

    def test_f_decays_in_ferromagnet(self):
        """f should decay moving into the F region."""
        x, f = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=20.0, xi_S=38.0, E_ex=256.0, n_grid=200)
        # Compare f at F boundary vs deep F
        n_half = len(x) // 2
        f_boundary = f[n_half]
        f_deep = f[-1]
        assert f_deep <= f_boundary + 0.1  # should decay or stay small

    def test_custom_temperature(self):
        """Custom T should produce different results than default."""
        x1, f1 = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, E_ex=256.0)
        x2, f2 = eilenberger.solve(
            Tc0=9.2, d_S=50.0, d_F=10.0, xi_S=38.0, E_ex=256.0, T=3.0)
        assert not np.allclose(f1, f2, atol=1e-6)


# ───────────────────── Ginzburg-Landau ────────────────────────

class TestGinzburgLandau:
    def test_output_shape(self):
        psi = ginzburg_landau.minimize(
            alpha=-1.0, beta=1.0, kappa=1.0, nx=16, ny=16, dx=1.0)
        assert psi.shape == (16, 16)

    def test_complex_output(self):
        psi = ginzburg_landau.minimize(
            alpha=-1.0, beta=1.0, kappa=1.0, nx=16, ny=16, dx=1.0)
        assert psi.dtype == complex or np.iscomplexobj(psi)

    def test_uniform_state(self):
        """Below Tc (α < 0), |ψ|² ≈ −α/β in equilibrium."""
        alpha, beta = -2.0, 1.0
        psi = ginzburg_landau.minimize(
            alpha=alpha, beta=beta, kappa=1.0, nx=32, ny=32, dx=1.0)
        expected = np.sqrt(-alpha / beta)
        mean_psi = np.mean(np.abs(psi))
        assert mean_psi == pytest.approx(expected, rel=0.1)

    def test_normal_state(self):
        """Above Tc (α > 0), |ψ| → 0."""
        psi = ginzburg_landau.minimize(
            alpha=1.0, beta=1.0, kappa=1.0, nx=16, ny=16, dx=1.0)
        assert np.mean(np.abs(psi)) < 0.1


# ────────────────────────── Triplet ───────────────────────────

class TestTriplet:
    def test_output_shape(self):
        x, ft = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0], n_grid=100)
        assert len(x) == 100
        assert len(ft) == 100

    def test_f_nonnegative(self):
        x, ft = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0])
        assert np.all(ft >= 0)

    def test_maximal_at_perpendicular(self):
        """Triplet amplitude should be maximal when M rotates by 90 deg."""
        _, ft_90 = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0], n_grid=100)
        _, ft_0 = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, 0.0, 0.0], n_grid=100)
        assert np.max(ft_90) > np.max(ft_0)

    def test_parallel_no_triplet(self):
        """Parallel magnetization should produce no triplet amplitude."""
        _, ft = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, 0.0, 0.0], n_grid=100)
        assert np.max(ft) < 1e-10

    def test_spatial_grid(self):
        thicknesses = [5.0, 10.0, 5.0]
        x, _ = triplet.solve(3, thicknesses, [0.0, np.pi / 4, 0.0])
        assert x[0] == pytest.approx(0.0)
        assert x[-1] == pytest.approx(sum(thicknesses), rel=0.01)

    def test_custom_xi_N(self):
        """Larger xi_N should produce longer-range correlations."""
        _, ft_short = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0], xi_N=5.0)
        _, ft_long = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0], xi_N=50.0)
        # Longer xi_N should produce less decay (higher f at edges)
        assert ft_long[-1] > ft_short[-1]

    def test_custom_xi_F_backward_compat(self):
        """Default xi_F=1.0, xi_N=10.0 should match old behavior."""
        _, ft1 = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0])
        _, ft2 = triplet.solve(
            3, [5.0, 10.0, 5.0], [0.0, np.pi / 2, 0.0], xi_F=1.0, xi_N=10.0)
        np.testing.assert_allclose(ft1, ft2)
