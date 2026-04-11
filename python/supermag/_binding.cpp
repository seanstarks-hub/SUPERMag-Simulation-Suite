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
#include "supermag/bdg.h"
#include "supermag/usadel.h"
#include "supermag/eilenberger.h"
#include "supermag/ginzburg_landau.h"
#include "supermag/josephson.h"
#include "supermag/triplet.h"

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
py_bdg_solve(int n_sites, double t_hop, double Delta, double E_ex) {
    int dim = 2 * n_sites;
    auto eigenvalues = py::array_t<double>(dim);
    auto ev_buf = eigenvalues.mutable_unchecked<1>();
    int n_eigenvalues = 0;

    int rc = supermag_bdg_solve(n_sites, t_hop, Delta, E_ex,
                                ev_buf.mutable_data(0), &n_eigenvalues);
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return eigenvalues;
}

// --- Usadel solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_usadel_solve(double Tc0, double d_S, double d_F,
                double xi_S, double xi_F, double E_ex,
                int n_grid) {
    auto Delta_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto d_buf = Delta_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    int rc = supermag_usadel_solve(Tc0, d_S, d_F, xi_S, xi_F, E_ex,
                                   n_grid, d_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, Delta_out};
}

// --- Eilenberger solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_eilenberger_solve(double Tc0, double d_S, double d_F,
                     double xi_S, double E_ex,
                     int n_grid) {
    auto f_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto f_buf = f_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    int rc = supermag_eilenberger_solve(Tc0, d_S, d_F, xi_S, E_ex,
                                        n_grid, f_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, f_out};
}

// --- Ginzburg-Landau solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_gl_minimize(double alpha, double beta, double kappa,
               int nx, int ny, double dx) {
    int N = nx * ny;
    auto psi_real = py::array_t<double>(N);
    auto psi_imag = py::array_t<double>(N);
    auto pr_buf = psi_real.mutable_unchecked<1>();
    auto pi_buf = psi_imag.mutable_unchecked<1>();

    int rc = supermag_gl_minimize(alpha, beta, kappa, nx, ny, dx,
                                  pr_buf.mutable_data(0), pi_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {psi_real, psi_imag};
}

// --- Josephson CPR solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_josephson_cpr(double d_F, double xi_F, double E_ex, double T,
                 int n_phases) {
    auto phase_arr = py::array_t<double>(n_phases);
    auto current_out = py::array_t<double>(n_phases);
    auto ph_buf = phase_arr.mutable_unchecked<1>();
    auto cur_buf = current_out.mutable_unchecked<1>();

    int rc = supermag_josephson_cpr(d_F, xi_F, E_ex, T,
                                    n_phases, ph_buf.mutable_data(0), cur_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {phase_arr, current_out};
}

// --- Triplet solver wrapper ---
static std::pair<py::array_t<double>, py::array_t<double>>
py_triplet_solve(int n_layers, py::array_t<double> thicknesses,
                 py::array_t<double> magnetization_angles,
                 int n_grid) {
    auto thick_buf = thicknesses.unchecked<1>();
    auto mag_buf = magnetization_angles.unchecked<1>();

    auto f_triplet_out = py::array_t<double>(n_grid);
    auto x_out = py::array_t<double>(n_grid);
    auto f_buf = f_triplet_out.mutable_unchecked<1>();
    auto x_buf = x_out.mutable_unchecked<1>();

    int rc = supermag_triplet_solve(n_layers, thick_buf.data(0), mag_buf.data(0),
                                    n_grid, f_buf.mutable_data(0), x_buf.mutable_data(0));
    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x_out, f_triplet_out};
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
          py::arg("E_ex"));

    m.def("_usadel_solve", &py_usadel_solve,
          "Solve Usadel equation for S/F bilayer",
          py::arg("Tc0"), py::arg("d_S"), py::arg("d_F"),
          py::arg("xi_S"), py::arg("xi_F"), py::arg("E_ex"),
          py::arg("n_grid"));

    m.def("_eilenberger_solve", &py_eilenberger_solve,
          "Solve Eilenberger equation for S/F bilayer",
          py::arg("Tc0"), py::arg("d_S"), py::arg("d_F"),
          py::arg("xi_S"), py::arg("E_ex"),
          py::arg("n_grid"));

    m.def("_gl_minimize", &py_gl_minimize,
          "Minimize Ginzburg-Landau free energy on 2D grid",
          py::arg("alpha"), py::arg("beta"), py::arg("kappa"),
          py::arg("nx"), py::arg("ny"), py::arg("dx"));

    m.def("_josephson_cpr", &py_josephson_cpr,
          "Compute Josephson current-phase relation",
          py::arg("d_F"), py::arg("xi_F"), py::arg("E_ex"),
          py::arg("T"), py::arg("n_phases"));

    m.def("_triplet_solve", &py_triplet_solve,
          "Compute spin-triplet pair correlations",
          py::arg("n_layers"), py::arg("thicknesses"),
          py::arg("magnetization_angles"),
          py::arg("n_grid"));
}
