// Tests for proximity optimizer, inverse solver, and fit functions.
#include "supermag/optimizer.h"
#include "supermag/proximity.h"
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstring>

static bool approx(double a, double b, double tol = 0.05) {
    return std::fabs(a - b) < tol * (1.0 + std::fabs(b));
}

static supermag_proximity_params_t make_params() {
    supermag_proximity_params_t p;
    std::memset(&p, 0, sizeof(p));
    p.Tc0 = 9.2;
    p.d_F = 5.0;
    p.xi_F = 3.0;
    p.gamma = 0.5;
    p.gamma_B = 0.3;
    p.E_ex = 100.0;
    p.model = SUPERMAG_MODEL_THIN_S;
    p.phase = SUPERMAG_PHASE_ZERO;
    p.geometry = SUPERMAG_GEOM_BILAYER;
    p.D_S = 1.0;
    p.geom_params = nullptr;
    p.spin_active = nullptr;
    return p;
}

int main() {
    std::printf("Running optimizer tests...\n");

    // ── Optimize null pointer checks ────────────────────────────
    {
        double d_F_out;
        int rc = supermag_optimize_tc(nullptr, nullptr,
                                       1.0, 20.0, 5.0, &d_F_out);
        assert(rc == SUPERMAG_ERR_NULL_PTR);

        auto p = make_params();
        rc = supermag_optimize_tc(&p, nullptr,
                                   1.0, 20.0, 5.0, nullptr);
        assert(rc == SUPERMAG_ERR_NULL_PTR);
        std::printf("  Optimize null checks: PASS\n");
    }

    // ── Optimize invalid range ──────────────────────────────────
    {
        auto p = make_params();
        double d_F_out;
        int rc = supermag_optimize_tc(&p, nullptr,
                                       20.0, 1.0, 5.0, &d_F_out);
        assert(rc == SUPERMAG_ERR_INVALID_DIM);
        std::printf("  Optimize invalid range: PASS\n");
    }

    // ── Optimize returns d_F within bounds ──────────────────────
    {
        auto p = make_params();
        double d_F_out;
        int rc = supermag_optimize_tc(&p, nullptr,
                                       1.0, 30.0, 5.0, &d_F_out);
        assert(rc == SUPERMAG_OK);
        assert(d_F_out >= 1.0 && d_F_out <= 30.0);
        std::printf("  Optimize d_F: result = %.4f  PASS\n", d_F_out);
    }

    // ── Inverse solver null checks ──────────────────────────────
    {
        double d_F_out;
        int rc = supermag_inverse_tc(nullptr, nullptr,
                                      5.0, 1.0, 20.0, &d_F_out);
        assert(rc == SUPERMAG_ERR_NULL_PTR);

        auto p = make_params();
        rc = supermag_inverse_tc(&p, nullptr,
                                  5.0, 1.0, 20.0, nullptr);
        assert(rc == SUPERMAG_ERR_NULL_PTR);
        std::printf("  Inverse null checks: PASS\n");
    }

    // ── Inverse invalid inputs ──────────────────────────────────
    {
        auto p = make_params();
        double d_F_out;
        int rc = supermag_inverse_tc(&p, nullptr,
                                      5.0, 20.0, 1.0, &d_F_out);
        assert(rc == SUPERMAG_ERR_INVALID_DIM);

        rc = supermag_inverse_tc(&p, nullptr,
                                  -1.0, 1.0, 20.0, &d_F_out);
        assert(rc == SUPERMAG_ERR_INVALID_DIM);
        std::printf("  Inverse invalid inputs: PASS\n");
    }

    // ── Inverse solver returns d_F within bounds ────────────────
    {
        auto p = make_params();
        double d_F_out;
        int rc = supermag_inverse_tc(&p, nullptr,
                                      5.0, 1.0, 30.0, &d_F_out);
        assert(rc == SUPERMAG_OK);
        assert(d_F_out >= 1.0 && d_F_out <= 30.0);
        std::printf("  Inverse d_F: result = %.4f  PASS\n", d_F_out);
    }

    // ── Fit null checks ─────────────────────────────────────────
    {
        double chi2;
        double d_F_data[] = {1.0, 5.0, 10.0};
        double Tc_data[] = {8.0, 6.0, 4.0};
        int rc = supermag_fit_tc(nullptr, nullptr,
                                  d_F_data, Tc_data, 3, 1, 1, 0, 0, &chi2);
        assert(rc == SUPERMAG_ERR_NULL_PTR);

        auto p = make_params();
        rc = supermag_fit_tc(&p, nullptr,
                              nullptr, Tc_data, 3, 1, 1, 0, 0, &chi2);
        assert(rc == SUPERMAG_ERR_NULL_PTR);

        rc = supermag_fit_tc(&p, nullptr,
                              d_F_data, nullptr, 3, 1, 1, 0, 0, &chi2);
        assert(rc == SUPERMAG_ERR_NULL_PTR);
        std::printf("  Fit null checks: PASS\n");
    }

    // ── Fit too few data points ─────────────────────────────────
    {
        auto p = make_params();
        double chi2;
        double d_F_data[] = {5.0};
        double Tc_data[] = {6.0};
        int rc = supermag_fit_tc(&p, nullptr,
                                  d_F_data, Tc_data, 1, 1, 1, 0, 0, &chi2);
        assert(rc == SUPERMAG_ERR_INVALID_DIM);
        std::printf("  Fit too few points: PASS\n");
    }

    // ── Fit returns finite chi2 ─────────────────────────────────
    {
        auto p = make_params();
        double chi2;
        double d_F_data[] = {2.0, 5.0, 10.0, 15.0, 20.0};
        double Tc_data[] = {8.5, 7.0, 5.0, 3.5, 2.5};
        int rc = supermag_fit_tc(&p, nullptr,
                                  d_F_data, Tc_data, 5, 1, 1, 0, 0, &chi2);
        assert(rc == SUPERMAG_OK);
        assert(std::isfinite(chi2));
        assert(chi2 >= 0.0);
        assert(p.gamma > 0.0);
        std::printf("  Fit: gamma = %.4f, chi2 = %.4f  PASS\n", p.gamma, chi2);
    }

    std::printf("All optimizer tests passed.\n");
    return 0;
}
