/**
 * pybind11 module wrapping SUPERMag C++ solvers for Python.
 *
 * Module name: _native
 * Wraps C-linkage functions from cpp/include/supermag/ into Python-callable
 * functions that accept and return NumPy arrays.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
#include <cstring>

// Include C-linkage headers (each has its own extern "C" guard)
#include "supermag/proximity.h"
#include "supermag/constants.h"
#include "supermag/error.h"
#include "supermag/depairing.h"
#include "supermag/optimizer.h"
#include "supermag/bdg.h"
#include "supermag/usadel.h"
#include "supermag/eilenberger.h"
#include "supermag/ginzburg_landau.h"
#include "supermag/josephson.h"
#include "supermag/triplet.h"
#include "supermag/solver_options.h"

namespace py = pybind11;

static std::pair<py::array_t<double>, py::array_t<double>>
py_pair_amplitude(double d_F, double xi_F, int phase, int n_points) {
    auto x = py::array_t<double>(n_points);
    auto F = py::array_t<double>(n_points);
    auto x_buf = x.mutable_unchecked<1>();
    auto F_buf = F.mutable_unchecked<1>();

    supermag_phase_t ph = (phase == 1) ? SUPERMAG_PHASE_PI : SUPERMAG_PHASE_ZERO;
    int rc = supermag_proximity_pair_amplitude(
        d_F, xi_F, ph, n_points,
        x_buf.mutable_data(0), F_buf.mutable_data(0));

    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x, F};
}

static py::array_t<double>
py_solve_tc_batch(double Tc0, double d_S, double xi_S, double xi_F,
                  double gamma, double gamma_B, double E_ex, double D_F,
                  int model, int phase,
                  double dp_ag, double dp_zeeman, double dp_orbital, double dp_spin_orbit,
                  py::array_t<double> d_F_arr) {
    auto d_F_buf = d_F_arr.unchecked<1>();
    int n = static_cast<int>(d_F_buf.shape(0));

    auto Tc = py::array_t<double>(n);
    auto Tc_buf = Tc.mutable_unchecked<1>();

    supermag_proximity_params_t params;
    std::memset(&params, 0, sizeof(params));
    params.Tc0 = Tc0;
    params.d_S = d_S;
    params.xi_S = xi_S;
    params.xi_F = xi_F;
    params.gamma = gamma;
    params.gamma_B = gamma_B;
    params.E_ex = E_ex;
    params.D_F = D_F;
    params.model = (model == 1) ? SUPERMAG_MODEL_FOMINOV : SUPERMAG_MODEL_THIN_S;
    params.phase = (phase == 1) ? SUPERMAG_PHASE_PI : SUPERMAG_PHASE_ZERO;

    supermag_depairing_t dp;
    dp.ag = dp_ag;
    dp.zeeman = dp_zeeman;
    dp.orbital = dp_orbital;
    dp.spin_orbit = dp_spin_orbit;

    bool has_dp = (dp_ag != 0.0 || dp_zeeman != 0.0 ||
                   dp_orbital != 0.0 || dp_spin_orbit != 0.0);

    int rc = supermag_proximity_solve_tc_batch(
        &params, d_F_buf.data(0), n,
        has_dp ? &dp : nullptr,
        Tc_buf.mutable_data(0));

    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return Tc;
}

// --- BdG solver wrapper ---
static py::array_t<double>
py_bdg_solve(int n_sites, double t_hop, double Delta, double E_ex, double mu) {
    int dim = 2 * n_sites;
    auto eigenvalues = py::array_t<double>(dim);
    auto ev_buf = eigenvalues.mutable_unchecked<1>();
    int n_eigenvalues = 0;

    int rc = supermag_bdg_solve(n_sites, t_hop, Delta, E_ex, mu,
                                ev_buf.mutable_data(0), &n_eigenvalues,
                                nullptr);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return eigenvalues;
}

// --- Usadel solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_usadel_solve(double Tc0, double d_S, double d_F,
                double xi_S, double xi_F, double E_ex,
                double T,
                int n_grid) {
    auto Delta_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto d_buf = Delta_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    int rc = supermag_usadel_solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex,
                                   T, SUPERMAG_USADEL_LINEARIZED,
                                   nullptr,
                                   n_grid,
                                   d_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, Delta_out};
}

// --- Eilenberger solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_eilenberger_solve(double Tc0, double d_S, double d_F,
                     double xi_S, double E_ex,
                     double T,
                     int n_grid) {
    auto f_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto f_buf = f_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    int rc = supermag_eilenberger_solve(Tc0, d_S, d_F, xi_S, E_ex,
                                        T, nullptr, n_grid,
                                        f_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, f_out};
}

// --- Ginzburg-Landau solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_gl_minimize(double alpha, double beta, double kappa,
               int nx, int ny, double dx,
               int mode = 0, double H_applied = 0.0) {
    int N = nx * ny;
    auto psi_real = py::array_t<double>(N);
    auto psi_imag = py::array_t<double>(N);
    auto pr_buf = psi_real.mutable_unchecked<1>();
    auto pi_buf = psi_imag.mutable_unchecked<1>();

    // Zero-initialize so C++ solver uses its own initial condition
    // (solver checks psi_real[0] != 0 for user-supplied IC)
    std::memset(pr_buf.mutable_data(0), 0, N * sizeof(double));
    std::memset(pi_buf.mutable_data(0), 0, N * sizeof(double));

    supermag_gl_mode_t gl_mode = (mode == 1) ? SUPERMAG_GL_GAUGE : SUPERMAG_GL_SCALAR;

    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                  gl_mode, H_applied, nullptr,
                                  pr_buf.mutable_data(0), pi_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {psi_real, psi_imag};
}

// --- Josephson CPR solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_josephson_cpr(double d_F, double xi_F, double E_ex, double T,
                 double Tc0, int n_phases, double gamma_B = 0.0) {
    const double pi = 3.14159265358979323846;
    auto phase_arr = py::array_t<double>(n_phases);
    auto current_out = py::array_t<double>(n_phases);
    auto ph_buf = phase_arr.mutable_unchecked<1>();
    auto cur_buf = current_out.mutable_unchecked<1>();

    // Generate uniform phase grid [0, 2π)
    for (int i = 0; i < n_phases; ++i)
        ph_buf(i) = 2.0 * pi * i / n_phases;

    int rc = supermag_josephson_cpr(d_F, xi_F, E_ex, T, Tc0, gamma_B,
                                    n_phases, ph_buf.data(0), nullptr,
                                    cur_buf.mutable_data(0), nullptr);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {phase_arr, current_out};
}

// --- Triplet solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_triplet_solve(int n_layers, py::array_t<double> thicknesses,
                 py::array_t<double> magnetization_angles,
                 double xi_F, double xi_N,
                 int n_grid, double T = 4.2, double Tc0 = 9.2,
                 py::object E_ex_per_layer_obj = py::none(),
                 py::object D_per_layer_obj = py::none(),
                 int mode = 0) {
    auto thick_buf = thicknesses.unchecked<1>();
    auto mag_buf = magnetization_angles.unchecked<1>();

    auto f_triplet_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto f_buf = f_triplet_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    const double* E_ex_ptr = nullptr;
    const double* D_ptr = nullptr;
    py::array_t<double> E_ex_arr, D_arr;

    if (!E_ex_per_layer_obj.is_none()) {
        E_ex_arr = E_ex_per_layer_obj.cast<py::array_t<double>>();
        E_ex_ptr = E_ex_arr.unchecked<1>().data(0);
    }
    if (!D_per_layer_obj.is_none()) {
        D_arr = D_per_layer_obj.cast<py::array_t<double>>();
        D_ptr = D_arr.unchecked<1>().data(0);
    }

    supermag_triplet_mode_t trip_mode =
        (mode == 1) ? SUPERMAG_TRIPLET_DIFFUSIVE : SUPERMAG_TRIPLET_PHENOMENOLOGICAL;

    int rc = supermag_triplet_solve(n_layers, thick_buf.data(0), mag_buf.data(0),
                                    E_ex_ptr, D_ptr,
                                    xi_F, xi_N, T, Tc0, trip_mode,
                                    n_grid, f_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, f_triplet_out};
}

// --- Depairing individual channel wrappers ---
static double py_depairing_ag(double gamma_s_meV, double T_kelvin) {
    return supermag_depairing_ag(gamma_s_meV, T_kelvin);
}

static double py_depairing_zeeman(double H_tesla, double T_kelvin) {
    return supermag_depairing_zeeman(H_tesla, T_kelvin);
}

static double py_depairing_orbital_perp(double D_nm2ps, double H_tesla,
                                         double thickness_nm, double T_kelvin) {
    return supermag_depairing_orbital_perp(D_nm2ps, H_tesla, thickness_nm, T_kelvin);
}

static double py_depairing_orbital_par(double D_nm2ps, double H_tesla,
                                        double thickness_nm, double T_kelvin) {
    return supermag_depairing_orbital_par(D_nm2ps, H_tesla, thickness_nm, T_kelvin);
}

static double py_depairing_soc(double Gamma_so_meV, double T_kelvin) {
    return supermag_depairing_soc(Gamma_so_meV, T_kelvin);
}

static py::tuple py_depairing_from_physical(double gamma_s_meV, double H_tesla,
                                             double D_nm2ps, double thickness_nm,
                                             double Gamma_so_meV, double T_kelvin) {
    supermag_depairing_t dp;
    int rc = supermag_depairing_from_physical(
        gamma_s_meV, H_tesla, D_nm2ps, thickness_nm,
        Gamma_so_meV, T_kelvin, &dp);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));
    return py::make_tuple(dp.ag, dp.zeeman, dp.orbital, dp.spin_orbit);
}

// --- Optimizer wrappers ---
static double py_optimize_tc(double Tc0, double d_S, double xi_S, double xi_F,
                              double gamma, double gamma_B, double E_ex, double D_F,
                              int model, int phase,
                              double dp_ag, double dp_zeeman,
                              double dp_orbital, double dp_spin_orbit,
                              double d_F_lo, double d_F_hi, double Tc_target) {
    supermag_proximity_params_t params;
    std::memset(&params, 0, sizeof(params));
    params.Tc0 = Tc0; params.d_S = d_S; params.xi_S = xi_S; params.xi_F = xi_F;
    params.gamma = gamma; params.gamma_B = gamma_B; params.E_ex = E_ex; params.D_F = D_F;
    params.model = (model == 1) ? SUPERMAG_MODEL_FOMINOV : SUPERMAG_MODEL_THIN_S;
    params.phase = (phase == 1) ? SUPERMAG_PHASE_PI : SUPERMAG_PHASE_ZERO;

    supermag_depairing_t dp = {dp_ag, dp_zeeman, dp_orbital, dp_spin_orbit};
    bool has_dp = (dp_ag != 0.0 || dp_zeeman != 0.0 ||
                   dp_orbital != 0.0 || dp_spin_orbit != 0.0);

    double d_F_out;
    int rc = supermag_optimize_tc(&params, has_dp ? &dp : nullptr,
                                  d_F_lo, d_F_hi, Tc_target, &d_F_out);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));
    return d_F_out;
}

static double py_inverse_tc(double Tc0, double d_S, double xi_S, double xi_F,
                             double gamma, double gamma_B, double E_ex, double D_F,
                             int model, int phase,
                             double dp_ag, double dp_zeeman,
                             double dp_orbital, double dp_spin_orbit,
                             double Tc_target, double d_F_lo, double d_F_hi) {
    supermag_proximity_params_t params;
    std::memset(&params, 0, sizeof(params));
    params.Tc0 = Tc0; params.d_S = d_S; params.xi_S = xi_S; params.xi_F = xi_F;
    params.gamma = gamma; params.gamma_B = gamma_B; params.E_ex = E_ex; params.D_F = D_F;
    params.model = (model == 1) ? SUPERMAG_MODEL_FOMINOV : SUPERMAG_MODEL_THIN_S;
    params.phase = (phase == 1) ? SUPERMAG_PHASE_PI : SUPERMAG_PHASE_ZERO;

    supermag_depairing_t dp = {dp_ag, dp_zeeman, dp_orbital, dp_spin_orbit};
    bool has_dp = (dp_ag != 0.0 || dp_zeeman != 0.0 ||
                   dp_orbital != 0.0 || dp_spin_orbit != 0.0);

    double d_F_out;
    int rc = supermag_inverse_tc(&params, has_dp ? &dp : nullptr,
                                 Tc_target, d_F_lo, d_F_hi, &d_F_out);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));
    return d_F_out;
}

static py::dict py_fit_tc(double Tc0, double d_S, double xi_S, double xi_F,
                           double gamma, double gamma_B, double E_ex, double D_F,
                           int model, int phase,
                           double dp_ag, double dp_zeeman,
                           double dp_orbital, double dp_spin_orbit,
                           py::array_t<double> d_F_data, py::array_t<double> Tc_data,
                           int fit_gamma, int fit_gamma_B,
                           int fit_E_ex, int fit_xi_F) {
    auto df_buf = d_F_data.unchecked<1>();
    auto tc_buf = Tc_data.unchecked<1>();
    int n = static_cast<int>(df_buf.shape(0));

    supermag_proximity_params_t params;
    std::memset(&params, 0, sizeof(params));
    params.Tc0 = Tc0; params.d_S = d_S; params.xi_S = xi_S; params.xi_F = xi_F;
    params.gamma = gamma; params.gamma_B = gamma_B; params.E_ex = E_ex; params.D_F = D_F;
    params.model = (model == 1) ? SUPERMAG_MODEL_FOMINOV : SUPERMAG_MODEL_THIN_S;
    params.phase = (phase == 1) ? SUPERMAG_PHASE_PI : SUPERMAG_PHASE_ZERO;

    supermag_depairing_t dp = {dp_ag, dp_zeeman, dp_orbital, dp_spin_orbit};
    bool has_dp = (dp_ag != 0.0 || dp_zeeman != 0.0 ||
                   dp_orbital != 0.0 || dp_spin_orbit != 0.0);

    double chi2;
    int rc = supermag_fit_tc(&params, has_dp ? &dp : nullptr,
                              df_buf.data(0), tc_buf.data(0), n,
                              fit_gamma, fit_gamma_B, fit_E_ex, fit_xi_F,
                              &chi2);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    py::dict result;
    result["gamma"] = params.gamma;
    result["gamma_B"] = params.gamma_B;
    result["E_ex"] = params.E_ex;
    result["xi_F"] = params.xi_F;
    result["chi2"] = chi2;
    return result;
}

PYBIND11_MODULE(_native, m) {
    m.doc() = "SUPERMag native C++ bindings";

    m.def("_pair_amplitude", &py_pair_amplitude,
          "Compute pair amplitude F(x) in the F layer",
          py::arg("d_F"), py::arg("xi_F"), py::arg("phase"), py::arg("n_points"));

    m.def("_solve_tc_batch", &py_solve_tc_batch,
          "Solve Tc(d_F) for S/F bilayer using digamma self-consistency",
          py::arg("Tc0"), py::arg("d_S"), py::arg("xi_S"), py::arg("xi_F"),
          py::arg("gamma"), py::arg("gamma_B"), py::arg("E_ex"), py::arg("D_F"),
          py::arg("model"), py::arg("phase"),
          py::arg("dp_ag"), py::arg("dp_zeeman"),
          py::arg("dp_orbital"), py::arg("dp_spin_orbit"),
          py::arg("d_F_arr"));

    m.def("_bdg_solve", &py_bdg_solve,
          "Diagonalize BdG Hamiltonian on tight-binding lattice",
          py::arg("n_sites"), py::arg("t_hop"), py::arg("Delta"),
          py::arg("E_ex"), py::arg("mu"));

    m.def("_usadel_solve", &py_usadel_solve,
          "Solve Usadel equation for S/F bilayer",
          py::arg("Tc0"), py::arg("d_S"), py::arg("d_F"),
          py::arg("xi_S"), py::arg("xi_F"), py::arg("E_ex"),
          py::arg("T"),
          py::arg("n_grid"));

    m.def("_eilenberger_solve", &py_eilenberger_solve,
          "Solve Eilenberger equation for S/F bilayer",
          py::arg("Tc0"), py::arg("d_S"), py::arg("d_F"),
          py::arg("xi_S"), py::arg("E_ex"),
          py::arg("T"),
          py::arg("n_grid"));

    m.def("_gl_minimize", &py_gl_minimize,
          "Minimize Ginzburg-Landau free energy on 2D grid",
          py::arg("alpha"), py::arg("beta"), py::arg("kappa"),
          py::arg("nx"), py::arg("ny"), py::arg("dx"),
          py::arg("mode") = 0,
          py::arg("H_applied") = 0.0);

    m.def("_josephson_cpr", &py_josephson_cpr,
          "Compute Josephson current-phase relation",
          py::arg("d_F"), py::arg("xi_F"), py::arg("E_ex"),
          py::arg("T"), py::arg("Tc0"), py::arg("n_phases"),
          py::arg("gamma_B") = 0.0);

    m.def("_triplet_solve", &py_triplet_solve,
          "Compute spin-triplet pair correlations",
          py::arg("n_layers"), py::arg("thicknesses"),
          py::arg("magnetization_angles"),
          py::arg("xi_F"), py::arg("xi_N"),
          py::arg("n_grid"), py::arg("T") = 4.2,
          py::arg("Tc0") = 9.2,
          py::arg("E_ex_per_layer") = py::none(),
          py::arg("D_per_layer") = py::none(),
          py::arg("mode") = 0);

    // Depairing individual channels
    m.def("_depairing_ag", &py_depairing_ag,
          "AG pair-breaking parameter",
          py::arg("gamma_s_meV"), py::arg("T_kelvin"));

    m.def("_depairing_zeeman", &py_depairing_zeeman,
          "Zeeman pair-breaking parameter",
          py::arg("H_tesla"), py::arg("T_kelvin"));

    m.def("_depairing_orbital_perp", &py_depairing_orbital_perp,
          "Orbital pair-breaking (perpendicular field)",
          py::arg("D_nm2ps"), py::arg("H_tesla"),
          py::arg("thickness_nm"), py::arg("T_kelvin"));

    m.def("_depairing_orbital_par", &py_depairing_orbital_par,
          "Orbital pair-breaking (parallel field)",
          py::arg("D_nm2ps"), py::arg("H_tesla"),
          py::arg("thickness_nm"), py::arg("T_kelvin"));

    m.def("_depairing_soc", &py_depairing_soc,
          "Spin-orbit coupling pair-breaking parameter",
          py::arg("Gamma_so_meV"), py::arg("T_kelvin"));

    m.def("_depairing_from_physical", &py_depairing_from_physical,
          "Compute all depairing channels from physical inputs",
          py::arg("gamma_s_meV"), py::arg("H_tesla"),
          py::arg("D_nm2ps"), py::arg("thickness_nm"),
          py::arg("Gamma_so_meV"), py::arg("T_kelvin"));

    // Optimizer utilities
    m.def("_optimize_tc", &py_optimize_tc,
          "Golden-section search for d_F matching target Tc",
          py::arg("Tc0"), py::arg("d_S"), py::arg("xi_S"), py::arg("xi_F"),
          py::arg("gamma"), py::arg("gamma_B"), py::arg("E_ex"), py::arg("D_F"),
          py::arg("model"), py::arg("phase"),
          py::arg("dp_ag"), py::arg("dp_zeeman"),
          py::arg("dp_orbital"), py::arg("dp_spin_orbit"),
          py::arg("d_F_lo"), py::arg("d_F_hi"), py::arg("Tc_target"));

    m.def("_inverse_tc", &py_inverse_tc,
          "Brent root-finding for d_F matching target Tc",
          py::arg("Tc0"), py::arg("d_S"), py::arg("xi_S"), py::arg("xi_F"),
          py::arg("gamma"), py::arg("gamma_B"), py::arg("E_ex"), py::arg("D_F"),
          py::arg("model"), py::arg("phase"),
          py::arg("dp_ag"), py::arg("dp_zeeman"),
          py::arg("dp_orbital"), py::arg("dp_spin_orbit"),
          py::arg("Tc_target"), py::arg("d_F_lo"), py::arg("d_F_hi"));

    m.def("_fit_tc", &py_fit_tc,
          "Nelder-Mead fit of proximity params to Tc(d_F) data",
          py::arg("Tc0"), py::arg("d_S"), py::arg("xi_S"), py::arg("xi_F"),
          py::arg("gamma"), py::arg("gamma_B"), py::arg("E_ex"), py::arg("D_F"),
          py::arg("model"), py::arg("phase"),
          py::arg("dp_ag"), py::arg("dp_zeeman"),
          py::arg("dp_orbital"), py::arg("dp_spin_orbit"),
          py::arg("d_F_data"), py::arg("Tc_data"),
          py::arg("fit_gamma"), py::arg("fit_gamma_B"),
          py::arg("fit_E_ex"), py::arg("fit_xi_F"));
}
