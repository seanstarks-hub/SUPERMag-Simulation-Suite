// Unit tests for proximity solver (pair_amplitude, solve_tc, solve_tc_batch)
#include "supermag/proximity.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstring>

// ---------- pair_amplitude tests (new signature) ----------

void test_pair_amplitude_zero_phase() {
    const int n = 100;
    std::vector<double> x(n), F(n);
    int rc = supermag_proximity_pair_amplitude(20.0, 5.0, SUPERMAG_PHASE_ZERO,
                                               n, x.data(), F.data());
    assert(rc == SUPERMAG_OK);
    // F(0) = exp(0)*cos(0) = 1.0
    assert(std::abs(F[0] - 1.0) < 1e-12);
    // F should decay overall
    assert(std::abs(F[n-1]) < std::abs(F[0]));
    // x spans [0, d_F]
    assert(std::abs(x[0]) < 1e-12);
    assert(std::abs(x[n-1] - 20.0) < 1e-12);
    std::printf("  PASS: test_pair_amplitude_zero_phase\n");
}

void test_pair_amplitude_pi_phase() {
    const int n = 100;
    std::vector<double> x(n), F(n);
    int rc = supermag_proximity_pair_amplitude(20.0, 5.0, SUPERMAG_PHASE_PI,
                                               n, x.data(), F.data());
    assert(rc == SUPERMAG_OK);
    // F(0) = exp(0)*sin(0) = 0.0
    assert(std::abs(F[0]) < 1e-12);
    // Should be nonzero at interior points
    assert(std::abs(F[n/4]) > 1e-6);
    std::printf("  PASS: test_pair_amplitude_pi_phase\n");
}

void test_pair_amplitude_null() {
    int rc = supermag_proximity_pair_amplitude(20.0, 5.0, SUPERMAG_PHASE_ZERO,
                                               10, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_pair_amplitude_null\n");
}

void test_pair_amplitude_invalid() {
    double x[2], F[2];
    // Negative d_F
    int rc = supermag_proximity_pair_amplitude(-1.0, 5.0, SUPERMAG_PHASE_ZERO, 2, x, F);
    assert(rc == SUPERMAG_ERR_INVALID_DIM);
    std::printf("  PASS: test_pair_amplitude_invalid\n");
}

// ---------- solve_tc tests ----------

static supermag_proximity_params_t make_default_params() {
    supermag_proximity_params_t p;
    std::memset(&p, 0, sizeof(p));
    p.Tc0 = 9.2;       // Nb bulk Tc (K)
    p.d_S = 50.0;      // nm
    p.d_F = 10.0;      // nm
    p.xi_S = 38.0;     // nm
    p.xi_F = 5.0;      // nm
    p.gamma = 0.3;
    p.gamma_B = 0.3;
    p.E_ex = 10.0;     // meV
    p.D_F = 2.5;       // nm^2/ps
    p.model = SUPERMAG_MODEL_THIN_S;
    p.phase = SUPERMAG_PHASE_ZERO;
    return p;
}

void test_solve_tc_thin_s() {
    supermag_proximity_params_t p = make_default_params();
    p.model = SUPERMAG_MODEL_THIN_S;
    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_OK);
    // Tc should be suppressed below Tc0 but remain positive
    assert(tc > 0.0);
    assert(tc <= p.Tc0);
    std::printf("  PASS: test_solve_tc_thin_s (Tc = %.4f K)\n", tc);
}

void test_solve_tc_fominov() {
    supermag_proximity_params_t p = make_default_params();
    p.model = SUPERMAG_MODEL_FOMINOV;
    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_OK);
    assert(tc >= 0.0);
    assert(tc <= p.Tc0);
    std::printf("  PASS: test_solve_tc_fominov (Tc = %.4f K)\n", tc);
}

void test_solve_tc_null() {
    double tc;
    int rc = supermag_proximity_solve_tc(nullptr, nullptr, &tc);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_solve_tc_null\n");
}

void test_solve_tc_invalid_model() {
    supermag_proximity_params_t p = make_default_params();
    p.model = static_cast<supermag_model_t>(99);
    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_ERR_INVALID_MODEL);
    std::printf("  PASS: test_solve_tc_invalid_model\n");
}

void test_solve_tc_with_depairing() {
    supermag_proximity_params_t p = make_default_params();
    supermag_depairing_t dp = {0.01, 0.0, 0.0, 0.0};
    double tc_no_dp, tc_with_dp;
    supermag_proximity_solve_tc(&p, nullptr, &tc_no_dp);
    supermag_proximity_solve_tc(&p, &dp, &tc_with_dp);
    // Depairing should further suppress Tc
    assert(tc_with_dp <= tc_no_dp + 1e-6);
    std::printf("  PASS: test_solve_tc_with_depairing (%.4f -> %.4f K)\n", tc_no_dp, tc_with_dp);
}

// ---------- solve_tc_batch tests ----------

void test_solve_tc_batch() {
    supermag_proximity_params_t p = make_default_params();
    p.model = SUPERMAG_MODEL_THIN_S;

    const int n = 20;
    std::vector<double> d_F(n), Tc(n);
    for (int i = 0; i < n; ++i) d_F[i] = 1.0 + i * 2.0;

    int rc = supermag_proximity_solve_tc_batch(&p, d_F.data(), n, nullptr, Tc.data());
    assert(rc == SUPERMAG_OK);
    for (int i = 0; i < n; ++i) {
        assert(Tc[i] >= 0.0);
        assert(Tc[i] <= p.Tc0 + 1e-10);
    }
    std::printf("  PASS: test_solve_tc_batch\n");
}

void test_solve_tc_batch_null() {
    int rc = supermag_proximity_solve_tc_batch(nullptr, nullptr, 1, nullptr, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_solve_tc_batch_null\n");
}

// ---------- depairing tests ----------

void test_depairing_total() {
    supermag_depairing_t dp = {0.1, 0.2, 0.3, 0.05};
    double total = supermag_depairing_total(&dp);
    assert(std::abs(total - 0.65) < 1e-12);
    std::printf("  PASS: test_depairing_total\n");
}

void test_depairing_null() {
    double total = supermag_depairing_total(nullptr);
    assert(total == 0.0);
    std::printf("  PASS: test_depairing_null\n");
}

int main() {
    std::printf("Running proximity tests...\n");

    // pair amplitude
    test_pair_amplitude_zero_phase();
    test_pair_amplitude_pi_phase();
    test_pair_amplitude_null();
    test_pair_amplitude_invalid();

    // solve_tc
    test_solve_tc_thin_s();
    test_solve_tc_fominov();
    test_solve_tc_null();
    test_solve_tc_invalid_model();
    test_solve_tc_with_depairing();

    // batch
    test_solve_tc_batch();
    test_solve_tc_batch_null();

    // depairing
    test_depairing_total();
    test_depairing_null();

    std::printf("All proximity tests passed!\n");
    return 0;
}
