// Tests for individual depairing channel functions.
#define _USE_MATH_DEFINES
#include "supermag/depairing.h"
#include "supermag/proximity.h"
#include "supermag/constants.h"
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstring>

static bool approx(double a, double b, double tol = 1e-8) {
    return std::fabs(a - b) < tol * (1.0 + std::fabs(b));
}

int main() {
    std::printf("Running depairing model tests...\n");

    double T = 9.2;  // Temperature in K
    double kB = supermag_const_kB();
    double mu_B = supermag_const_mu_B();
    double hbar = supermag_const_hbar();
    double e = supermag_const_e();

    // ── AG depairing ────────────────────────────────────────────
    {
        double gamma_s = 0.1;  // meV
        double lambda = supermag_depairing_ag(gamma_s, T);
        double expected = gamma_s / (2.0 * kB * T);
        assert(approx(lambda, expected));
        assert(lambda > 0.0);
        std::printf("  AG depairing: lambda = %.6f  PASS\n", lambda);
    }

    // AG with T <= 0 returns 0
    assert(supermag_depairing_ag(0.1, 0.0) == 0.0);
    assert(supermag_depairing_ag(0.1, -1.0) == 0.0);

    // ── Zeeman depairing ────────────────────────────────────────
    {
        double H = 1.0;  // Tesla
        double lambda = supermag_depairing_zeeman(H, T);
        double denom = 2.0 * M_PI * kB * T;
        double expected = (mu_B * H) * (mu_B * H) / (denom * denom);
        assert(approx(lambda, expected));
        assert(lambda > 0.0);
        std::printf("  Zeeman depairing: lambda = %.6e  PASS\n", lambda);
    }

    assert(supermag_depairing_zeeman(1.0, 0.0) == 0.0);

    // ── Orbital depairing (perpendicular) ───────────────────────
    {
        double D = 1.0;         // diffusion coefficient (nm^2/ps)
        double H = 1.0;         // Tesla
        double thickness = 10.0; // nm
        double lambda = supermag_depairing_orbital_perp(D, H, thickness, T);
        double denom = 3.0 * hbar * hbar * 2.0 * M_PI * kB * T;
        double expected = D * (e * H) * (e * H) * thickness * thickness / denom;
        assert(approx(lambda, expected));
        assert(lambda > 0.0);
        std::printf("  Orbital perp depairing: lambda = %.6e  PASS\n", lambda);
    }

    assert(supermag_depairing_orbital_perp(1.0, 1.0, 10.0, 0.0) == 0.0);

    // ── Orbital depairing (parallel) ────────────────────────────
    {
        double D = 1.0;
        double H = 1.0;
        double thickness = 10.0;
        double lambda_perp = supermag_depairing_orbital_perp(D, H, thickness, T);
        double lambda_par = supermag_depairing_orbital_par(D, H, thickness, T);
        // Parallel should be 1/4 of perpendicular (factor 12 vs 3)
        assert(approx(lambda_par, lambda_perp / 4.0));
        assert(lambda_par > 0.0);
        std::printf("  Orbital par depairing: lambda = %.6e  PASS\n", lambda_par);
    }

    assert(supermag_depairing_orbital_par(1.0, 1.0, 10.0, 0.0) == 0.0);

    // ── Spin-orbit coupling depairing ───────────────────────────
    {
        double Gamma_so = 0.05;
        double lambda = supermag_depairing_soc(Gamma_so, T);
        double expected = Gamma_so / (2.0 * kB * T);
        assert(approx(lambda, expected));
        assert(lambda > 0.0);
        std::printf("  SOC depairing: lambda = %.6f  PASS\n", lambda);
    }

    assert(supermag_depairing_soc(0.05, 0.0) == 0.0);

    // ── Consistency: AG and SOC use same formula ────────────────
    {
        double rate = 0.2;
        double ag = supermag_depairing_ag(rate, T);
        double soc = supermag_depairing_soc(rate, T);
        assert(approx(ag, soc));
        std::printf("  AG == SOC for same rate: PASS\n");
    }

    // ── Zero fields give zero depairing ─────────────────────────
    assert(supermag_depairing_ag(0.0, T) == 0.0);
    assert(supermag_depairing_zeeman(0.0, T) == 0.0);
    assert(supermag_depairing_orbital_perp(1.0, 0.0, 10.0, T) == 0.0);
    assert(supermag_depairing_orbital_par(1.0, 0.0, 10.0, T) == 0.0);
    assert(supermag_depairing_soc(0.0, T) == 0.0);

    // ── from_physical fills all channels ────────────────────────
    {
        supermag_depairing_t out;
        std::memset(&out, 0, sizeof(out));
        int rc = supermag_depairing_from_physical(
            0.1, 1.0, 2.5, 50.0, 0.05, T, &out);
        assert(rc == SUPERMAG_OK);
        assert(out.ag > 0.0);
        assert(out.zeeman > 0.0);
        assert(out.orbital > 0.0);
        assert(out.spin_orbit > 0.0);
        std::printf("  from_physical: AG=%.4e Z=%.4e O=%.4e SO=%.4e  PASS\n",
                    out.ag, out.zeeman, out.orbital, out.spin_orbit);
    }

    // from_physical null output
    {
        int rc = supermag_depairing_from_physical(0.1, 1.0, 2.5, 50.0, 0.05, T, nullptr);
        assert(rc == SUPERMAG_ERR_NULL_PTR);
        std::printf("  from_physical null: PASS\n");
    }

    std::printf("All depairing model tests passed.\n");
    return 0;
}
