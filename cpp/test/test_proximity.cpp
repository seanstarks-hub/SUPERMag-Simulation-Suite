// Unit tests for proximity solver (pair_amplitude, solve_tc, solve_tc_batch)
#include "supermag/proximity.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstring>
#include <complex>

// Forward-declare internal kernel functions for overflow tests
namespace supermag {
std::complex<double> kernel_coth(double d_F, double xi_F);
std::complex<double> kernel_tanh(double d_F, double xi_F);
}

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
    p.geometry = SUPERMAG_GEOM_BILAYER;
    p.geom_params = nullptr;
    p.spin_active = nullptr;
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

// ---------- overflow-safe kernel tests ----------

void test_kernel_overflow_coth() {
    // d_F/xi_F = 2000 would overflow raw sinh/cosh
    auto K = supermag::kernel_coth(2000.0, 1.0);
    assert(std::isfinite(K.real()));
    assert(std::isfinite(K.imag()));
    // In the asymptotic regime coth(z) → 1, so K → q = (1+i)/xi_F
    auto q = std::complex<double>(1.0, 1.0) / 1.0;
    assert(std::abs(K.real() - q.real()) < 1e-10);
    assert(std::abs(K.imag() - q.imag()) < 1e-10);
    std::printf("  PASS: test_kernel_overflow_coth\n");
}

void test_kernel_overflow_tanh() {
    auto K = supermag::kernel_tanh(2000.0, 1.0);
    assert(std::isfinite(K.real()));
    assert(std::isfinite(K.imag()));
    auto q = std::complex<double>(1.0, 1.0) / 1.0;
    assert(std::abs(K.real() - q.real()) < 1e-10);
    assert(std::abs(K.imag() - q.imag()) < 1e-10);
    std::printf("  PASS: test_kernel_overflow_tanh\n");
}

void test_solve_tc_large_df() {
    supermag_proximity_params_t p = make_default_params();
    p.d_F = 5000.0;  // extreme d_F >> xi_F
    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_OK);
    assert(std::isfinite(tc));
    assert(tc >= 0.0);
    assert(tc <= p.Tc0);
    std::printf("  PASS: test_solve_tc_large_df (Tc = %.4f K)\n", tc);
}

// ---------- physics-based depairing tests ----------

void test_depairing_compute_zero_field() {
    supermag_depairing_input_t input;
    std::memset(&input, 0, sizeof(input));
    input.Tc0 = 9.2;
    input.Gamma_s = 0.0;
    input.H = 0.0;
    input.D = 2.5;
    input.thickness = 50.0;
    input.Gamma_so = 0.0;

    supermag_depairing_t out;
    int rc = supermag_depairing_compute(&input, &out);
    assert(rc == SUPERMAG_OK);
    // Zero field → zero Zeeman and orbital
    assert(std::abs(out.zeeman) < 1e-15);
    assert(std::abs(out.orbital) < 1e-15);
    // Zero spin-flip/spin-orbit rates → zero AG and SO
    assert(std::abs(out.ag) < 1e-15);
    assert(std::abs(out.spin_orbit) < 1e-15);
    std::printf("  PASS: test_depairing_compute_zero_field\n");
}

void test_depairing_compute_known_values() {
    supermag_depairing_input_t input;
    std::memset(&input, 0, sizeof(input));
    input.Tc0 = 9.2;        // K (Nb)
    input.Gamma_s = 0.1;    // meV
    input.H = 1.0;          // Tesla
    input.D = 2.5;          // nm^2/ps
    input.thickness = 50.0; // nm
    input.Gamma_so = 0.05;  // meV

    supermag_depairing_t out;
    int rc = supermag_depairing_compute(&input, &out);
    assert(rc == SUPERMAG_OK);
    // All channels should be positive
    assert(out.ag > 0.0);
    assert(out.zeeman > 0.0);
    assert(out.orbital > 0.0);
    assert(out.spin_orbit > 0.0);
    // AG should be proportional to Gamma_s
    assert(out.ag > out.spin_orbit * 1.5);  // Gamma_s > Gamma_so
    std::printf("  PASS: test_depairing_compute_known_values (AG=%.4e, Z=%.4e, O=%.4e, SO=%.4e)\n",
                out.ag, out.zeeman, out.orbital, out.spin_orbit);
}

void test_depairing_compute_null() {
    supermag_depairing_t out;
    int rc = supermag_depairing_compute(nullptr, &out);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    supermag_depairing_input_t input;
    std::memset(&input, 0, sizeof(input));
    input.Tc0 = 9.2;
    rc = supermag_depairing_compute(&input, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_depairing_compute_null\n");
}

void test_depairing_roundtrip() {
    supermag_depairing_input_t input;
    std::memset(&input, 0, sizeof(input));
    input.Tc0 = 9.2;
    input.Gamma_s = 0.1;
    input.H = 0.5;
    input.D = 3.0;
    input.thickness = 30.0;
    input.Gamma_so = 0.02;

    supermag_depairing_t out;
    int rc = supermag_depairing_compute(&input, &out);
    assert(rc == SUPERMAG_OK);
    double total = supermag_depairing_total(&out);
    double sum = out.ag + out.zeeman + out.orbital + out.spin_orbit;
    assert(std::abs(total - sum) < 1e-15);
    std::printf("  PASS: test_depairing_roundtrip (total=%.4e)\n", total);
}

// ---------- kernel integration tests (new model dispatch) ----------

void test_solve_tc_trilayer() {
    supermag_proximity_params_t p = make_default_params();
    p.geometry = SUPERMAG_GEOM_TRILAYER;

    supermag_trilayer_params_t tri;
    std::memset(&tri, 0, sizeof(tri));
    tri.d_N = 5.0;       // 5 nm normal metal interlayer
    tri.xi_N = 50.0;     // long coherence length
    tri.R_B = 0.0;
    p.geom_params = &tri;

    double tc_tri;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc_tri);
    assert(rc == SUPERMAG_OK);
    assert(tc_tri >= 0.0);
    assert(tc_tri <= p.Tc0);

    // Compare with bilayer — trilayer with N interlayer should differ
    supermag_proximity_params_t p2 = make_default_params();
    double tc_bi;
    supermag_proximity_solve_tc(&p2, nullptr, &tc_bi);
    // Tc values should be distinct (N layer modifies kernel)
    assert(std::abs(tc_tri - tc_bi) > 1e-6 || tc_tri == 0.0 || tc_bi == 0.0);

    std::printf("  PASS: test_solve_tc_trilayer (Tc=%.4f vs bilayer %.4f K)\n", tc_tri, tc_bi);
}

void test_solve_tc_graded() {
    supermag_proximity_params_t p = make_default_params();
    p.geometry = SUPERMAG_GEOM_GRADED;

    supermag_graded_params_t grade;
    std::memset(&grade, 0, sizeof(grade));
    grade.E_ex_surface = 10.0;    // meV at S/F interface
    grade.E_ex_bulk = 5.0;    // meV at vacuum surface
    grade.n_slabs = 10;
    grade.profile = SUPERMAG_GRADE_LINEAR;
    p.geom_params = &grade;

    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_OK);
    assert(tc >= 0.0);
    assert(tc <= p.Tc0);
    std::printf("  PASS: test_solve_tc_graded (Tc=%.4f K)\n", tc);
}

void test_solve_tc_domains() {
    supermag_proximity_params_t p = make_default_params();
    p.geometry = SUPERMAG_GEOM_DOMAINS;

    supermag_domain_params_t dom;
    std::memset(&dom, 0, sizeof(dom));
    dom.domain_width = 5.0;
    p.d_F = 10.0;  // total thickness = 2 domains
    p.geom_params = &dom;

    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_OK);
    assert(tc >= 0.0);
    assert(tc <= p.Tc0);
    std::printf("  PASS: test_solve_tc_domains (Tc=%.4f K)\n", tc);
}

void test_solve_tc_trilayer_null_params() {
    supermag_proximity_params_t p = make_default_params();
    p.geometry = SUPERMAG_GEOM_TRILAYER;
    p.geom_params = nullptr;  // should error
    double tc;
    int rc = supermag_proximity_solve_tc(&p, nullptr, &tc);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_solve_tc_trilayer_null_params\n");
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

    // overflow-safe kernels
    test_kernel_overflow_coth();
    test_kernel_overflow_tanh();
    test_solve_tc_large_df();

    // physics-based depairing
    test_depairing_compute_zero_field();
    test_depairing_compute_known_values();
    test_depairing_compute_null();
    test_depairing_roundtrip();

    // kernel integration (model dispatch)
    test_solve_tc_trilayer();
    test_solve_tc_graded();
    test_solve_tc_domains();
    test_solve_tc_trilayer_null_params();

    std::printf("All proximity tests passed!\n");
    return 0;
}
