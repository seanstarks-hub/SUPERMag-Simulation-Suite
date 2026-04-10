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

// Include C-linkage headers
extern "C" {
#include "supermag/proximity.h"
#include "supermag/constants.h"
#include "supermag/error.h"
}

namespace py = pybind11;

static std::pair<py::array_t<double>, py::array_t<double>>
py_pair_amplitude(double F0, double xi_F, double d_F, int n_points) {
    auto x = py::array_t<double>(n_points);
    auto F = py::array_t<double>(n_points);
    auto x_buf = x.mutable_unchecked<1>();
    auto F_buf = F.mutable_unchecked<1>();

    int rc = supermag_proximity_pair_amplitude(
        F0, xi_F, d_F, n_points, x_buf.mutable_data(0), F_buf.mutable_data(0));

    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return {x, F};
}

static py::array_t<double>
py_critical_temp(double Tc0, double d_S, double xi_S, double xi_F, double E_ex,
                 py::array_t<double> d_F_arr) {
    auto d_F_buf = d_F_arr.unchecked<1>();
    int n = static_cast<int>(d_F_buf.shape(0));

    auto Tc = py::array_t<double>(n);
    auto Tc_buf = Tc.mutable_unchecked<1>();

    int rc = supermag_proximity_critical_temp(
        Tc0, d_S, xi_S, xi_F, E_ex,
        d_F_buf.data(0), n, Tc_buf.mutable_data(0));

    if (rc != SUPERMAG_OK)
        throw std::runtime_error(supermag_error_string(rc));

    return Tc;
}

PYBIND11_MODULE(_native, m) {
    m.doc() = "SUPERMag native C++ bindings";

    m.def("_pair_amplitude", &py_pair_amplitude,
          "Compute pair amplitude F(x) in the F layer",
          py::arg("F0"), py::arg("xi_F"), py::arg("d_F"), py::arg("n_points"));

    m.def("_critical_temp", &py_critical_temp,
          "Compute Tc(d_F) for S/F bilayer",
          py::arg("Tc0"), py::arg("d_S"), py::arg("xi_S"),
          py::arg("xi_F"), py::arg("E_ex"), py::arg("d_F_arr"));
}
