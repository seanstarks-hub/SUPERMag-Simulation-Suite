// Optimizer, inverse solver, and Nelder-Mead fit for proximity Tc.
//
// supermag_optimize_tc:  Golden-section search on d_F to match Tc_target
// supermag_inverse_tc:   Brent inversion for d_F given target Tc
// supermag_fit_tc:        Nelder-Mead for multi-parameter least-squares fit

#include "supermag/optimizer.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>

extern "C" {

// ── Golden-section optimizer for d_F ────────────────────────────

int supermag_optimize_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double d_F_lo, double d_F_hi,
    double Tc_target,
    double *d_F_out)
{
    if (!params || !d_F_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (d_F_lo >= d_F_hi)
        return SUPERMAG_ERR_INVALID_DIM;

    const double gr = (std::sqrt(5.0) + 1.0) / 2.0;
    double a = d_F_lo, b = d_F_hi;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    auto objective = [&](double val) -> double {
        params->d_F = val;
        double tc;
        int rc = supermag_proximity_solve_tc(params, depairing, &tc);
        if (rc != SUPERMAG_OK) return 1e30;
        return std::fabs(tc - Tc_target);
    };

    for (int iter = 0; iter < 100; ++iter) {
        if (std::fabs(b - a) < 1e-10) break;
        if (objective(c) < objective(d))
            b = d;
        else
            a = c;
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }

    *d_F_out = (a + b) / 2.0;
    params->d_F = *d_F_out;
    return SUPERMAG_OK;
}

// ── Brent inversion for d_F ─────────────────────────────────────

int supermag_inverse_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    double Tc_target,
    double d_F_lo, double d_F_hi,
    double *d_F_out)
{
    if (!params || !d_F_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (d_F_lo >= d_F_hi || Tc_target <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Brent's method on f(d_F) = Tc(d_F) - Tc_target
    double a = d_F_lo, b = d_F_hi;
    double tc_a, tc_b;

    params->d_F = a;
    supermag_proximity_solve_tc(params, depairing, &tc_a);
    double fa = tc_a - Tc_target;

    params->d_F = b;
    supermag_proximity_solve_tc(params, depairing, &tc_b);
    double fb = tc_b - Tc_target;

    if (fa * fb > 0.0) {
        // No sign change: return best endpoint
        *d_F_out = (std::fabs(fa) < std::fabs(fb)) ? a : b;
        params->d_F = *d_F_out;
        return SUPERMAG_OK;
    }

    double c_val = a, fc = fa;
    double d_val = 0.0, e_val = 0.0;
    bool mflag = true;

    if (std::fabs(fc) < std::fabs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    c_val = a; fc = fa;
    d_val = b - a;
    e_val = d_val;

    for (int iter = 0; iter < 100; ++iter) {
        if (std::fabs(fb) < 1e-12 || std::fabs(b - a) < 1e-12) break;

        double s;
        if (std::fabs(fa - fc) > 1e-30 && std::fabs(fb - fc) > 1e-30) {
            // Inverse quadratic interpolation
            s = a*fb*fc / ((fa-fb)*(fa-fc))
              + b*fa*fc / ((fb-fa)*(fb-fc))
              + c_val*fa*fb / ((fc-fa)*(fc-fb));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        // Conditions for bisection fallback
        bool cond1 = (s < std::min((3*a+b)/4.0, b) || s > std::max((3*a+b)/4.0, b));
        bool cond2 = (mflag && std::fabs(s - b) >= std::fabs(b - c_val) / 2.0);
        bool cond3 = (!mflag && std::fabs(s - b) >= std::fabs(c_val - d_val) / 2.0);

        if (cond1 || cond2 || cond3) {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }

        params->d_F = s;
        double tc_s;
        supermag_proximity_solve_tc(params, depairing, &tc_s);
        double fs = tc_s - Tc_target;

        d_val = c_val;
        c_val = b; fc = fb;

        if (fa * fs < 0.0) {
            b = s; fb = fs;
        } else {
            a = s; fa = fs;
        }

        if (std::fabs(fa) < std::fabs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    *d_F_out = b;
    params->d_F = *d_F_out;
    return SUPERMAG_OK;
}

// ── Nelder-Mead fit ─────────────────────────────────────────────

// Helper: get/set a fitted parameter by index in the active set
static double* get_param_ptr(supermag_proximity_params_t *p, int which) {
    switch (which) {
        case 0: return &p->gamma;
        case 1: return &p->gamma_B;
        case 2: return &p->E_ex;
        case 3: return &p->xi_F;
        default: return nullptr;
    }
}

int supermag_fit_tc(
    supermag_proximity_params_t *params,
    const supermag_depairing_t *depairing,
    const double *d_F_data, const double *Tc_data, int n_data,
    int fit_gamma, int fit_gamma_B, int fit_E_ex, int fit_xi_F,
    double *chi2_out)
{
    if (!params || !d_F_data || !Tc_data)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_data < 2)
        return SUPERMAG_ERR_INVALID_DIM;

    // Build list of active parameter indices
    std::vector<int> active;
    if (fit_gamma)   active.push_back(0);
    if (fit_gamma_B) active.push_back(1);
    if (fit_E_ex)    active.push_back(2);
    if (fit_xi_F)    active.push_back(3);

    int ndim = static_cast<int>(active.size());
    if (ndim == 0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Chi^2 objective function
    auto chi2 = [&](const std::vector<double>& x) -> double {
        for (int i = 0; i < ndim; ++i)
            *get_param_ptr(params, active[i]) = x[i];

        double sum = 0.0;
        for (int k = 0; k < n_data; ++k) {
            params->d_F = d_F_data[k];
            double tc;
            int rc = supermag_proximity_solve_tc(params, depairing, &tc);
            if (rc != SUPERMAG_OK) return 1e30;
            double diff = tc - Tc_data[k];
            sum += diff * diff;
        }
        return sum;
    };

    // Initialize simplex from current parameter values
    int n_vertices = ndim + 1;
    std::vector<std::vector<double>> simplex(n_vertices, std::vector<double>(ndim));
    std::vector<double> f_vals(n_vertices);

    // Vertex 0: current values
    for (int i = 0; i < ndim; ++i)
        simplex[0][i] = *get_param_ptr(params, active[i]);
    f_vals[0] = chi2(simplex[0]);

    // Remaining vertices: perturb each dimension by 10%
    for (int j = 0; j < ndim; ++j) {
        simplex[j+1] = simplex[0];
        double step = simplex[0][j] * 0.1;
        if (std::fabs(step) < 1e-6) step = 0.1;
        simplex[j+1][j] += step;
        f_vals[j+1] = chi2(simplex[j+1]);
    }

    // Nelder-Mead iteration
    const double alpha = 1.0, gamma_nm = 2.0, rho = 0.5, sigma = 0.5;
    const int max_iter = 500;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Sort vertices by function value
        for (int i = 0; i < n_vertices - 1; ++i)
            for (int j = i + 1; j < n_vertices; ++j)
                if (f_vals[j] < f_vals[i]) {
                    std::swap(simplex[i], simplex[j]);
                    std::swap(f_vals[i], f_vals[j]);
                }

        // Check convergence
        double range = f_vals[n_vertices-1] - f_vals[0];
        if (range < tol) break;

        // Centroid (excluding worst)
        std::vector<double> centroid(ndim, 0.0);
        for (int i = 0; i < n_vertices - 1; ++i)
            for (int d = 0; d < ndim; ++d)
                centroid[d] += simplex[i][d];
        for (int d = 0; d < ndim; ++d)
            centroid[d] /= (n_vertices - 1);

        // Reflection
        std::vector<double> xr(ndim);
        for (int d = 0; d < ndim; ++d)
            xr[d] = centroid[d] + alpha * (centroid[d] - simplex[n_vertices-1][d]);
        double fr = chi2(xr);

        if (fr >= f_vals[0] && fr < f_vals[n_vertices-2]) {
            simplex[n_vertices-1] = xr;
            f_vals[n_vertices-1] = fr;
        } else if (fr < f_vals[0]) {
            // Expansion
            std::vector<double> xe(ndim);
            for (int d = 0; d < ndim; ++d)
                xe[d] = centroid[d] + gamma_nm * (xr[d] - centroid[d]);
            double fe = chi2(xe);
            if (fe < fr) {
                simplex[n_vertices-1] = xe;
                f_vals[n_vertices-1] = fe;
            } else {
                simplex[n_vertices-1] = xr;
                f_vals[n_vertices-1] = fr;
            }
        } else {
            // Contraction
            std::vector<double> xc(ndim);
            for (int d = 0; d < ndim; ++d)
                xc[d] = centroid[d] + rho * (simplex[n_vertices-1][d] - centroid[d]);
            double fc = chi2(xc);
            if (fc < f_vals[n_vertices-1]) {
                simplex[n_vertices-1] = xc;
                f_vals[n_vertices-1] = fc;
            } else {
                // Shrink
                for (int i = 1; i < n_vertices; ++i) {
                    for (int d = 0; d < ndim; ++d)
                        simplex[i][d] = simplex[0][d] + sigma * (simplex[i][d] - simplex[0][d]);
                    f_vals[i] = chi2(simplex[i]);
                }
            }
        }
    }

    // Best vertex
    int best = 0;
    for (int i = 1; i < n_vertices; ++i)
        if (f_vals[i] < f_vals[best]) best = i;

    for (int i = 0; i < ndim; ++i)
        *get_param_ptr(params, active[i]) = simplex[best][i];

    if (chi2_out)
        *chi2_out = f_vals[best];

    return SUPERMAG_OK;
}

} // extern "C"
