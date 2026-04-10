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

// Include C-linkage headers
extern "C" {
#include "supermag/proximity.h"
#include "supermag/constants.h"
#include "supermag/error.h"
}

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
}
