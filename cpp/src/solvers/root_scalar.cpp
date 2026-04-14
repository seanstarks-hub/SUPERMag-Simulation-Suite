// Root-finding for scalar equations F(T) = 0
// Strategy: coarse grid scan → bracket sign changes → Brent's method refinement
// Returns the highest root (the physical Tc).

#include "root_scalar.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace supermag {

// Brent's method for finding a root of f in [a, b] where f(a) and f(b) have opposite signs.
static double brent(double (*f)(double, void*), void *ctx,
                    double a, double b, double tol) {
    double fa = f(a, ctx);
    double fb = f(b, ctx);

    double c = a, fc = fa;
    double d = b - a, e = d;

    for (int iter = 0; iter < 200; ++iter) {
        if (fb * fc > 0.0) {
            c = a; fc = fa;
            d = b - a; e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol1 = 2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) + 0.5 * tol;
        double m = 0.5 * (c - b);

        if (std::abs(m) <= tol1 || fb == 0.0) {
            return b;
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            double s = fb / fa;
            double p, q;
            if (a == c) {
                // Secant method
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation
                double r = fb / fc;
                q = fa / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;

            if (2.0 * p < 3.0 * m * q - std::abs(tol1 * q) &&
                2.0 * p < std::abs(e * q)) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = m;
            }
        } else {
            d = m;
            e = m;
        }

        a = b; fa = fb;
        if (std::abs(d) > tol1) {
            b += d;
        } else {
            b += (m > 0.0 ? tol1 : -tol1);
        }
        fb = f(b, ctx);
    }

    return b;
}

double root_scalar_solve(
    double (*f)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol)
{
    const int N = 1000;
    double dT = (T_max - T_min) / N;

    // Evaluate on coarse grid
    std::vector<double> T_grid(N + 1);
    std::vector<double> f_grid(N + 1);
    for (int i = 0; i <= N; ++i) {
        T_grid[i] = T_min + i * dT;
        f_grid[i] = f(T_grid[i], context);
    }

    // Find all sign-change brackets and exact grid-point roots
    struct Bracket { double a, b; };
    std::vector<Bracket> brackets;
    double highest_root = -std::numeric_limits<double>::infinity();
    bool found_exact = false;

    for (int i = 0; i <= N; ++i) {
        // Check for exact zero at grid point
        if (f_grid[i] == 0.0) {
            if (T_grid[i] > highest_root) {
                highest_root = T_grid[i];
                found_exact = true;
            }
        }
        // Check for sign change between consecutive grid points
        if (i < N && f_grid[i] * f_grid[i + 1] < 0.0) {
            brackets.push_back({T_grid[i], T_grid[i + 1]});
        }
    }

    // Refine each bracket and track highest root
    for (auto &br : brackets) {
        double root = brent(f, context, br.a, br.b, tol);
        if (root > highest_root) {
            highest_root = root;
            found_exact = true;
        }
    }

    if (!found_exact) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return highest_root;
}

double root_scalar_solve_log(
    double (*f)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol)
{
    // Log-spaced grid gives much better resolution at low T.
    // At T=0.1 K with Tc0=9 K, spacing is ~0.001 K vs ~0.009 K for linear.
    const int N = 1000;
    double log_min = std::log(T_min);
    double log_max = std::log(T_max);
    double d_log = (log_max - log_min) / N;

    std::vector<double> T_grid(N + 1);
    std::vector<double> f_grid(N + 1);
    for (int i = 0; i <= N; ++i) {
        T_grid[i] = std::exp(log_min + i * d_log);
        f_grid[i] = f(T_grid[i], context);
    }

    struct Bracket { double a, b; };
    std::vector<Bracket> brackets;
    double highest_root = -std::numeric_limits<double>::infinity();
    bool found_exact = false;

    for (int i = 0; i <= N; ++i) {
        if (f_grid[i] == 0.0) {
            if (T_grid[i] > highest_root) {
                highest_root = T_grid[i];
                found_exact = true;
            }
        }
        if (i < N && f_grid[i] * f_grid[i + 1] < 0.0) {
            brackets.push_back({T_grid[i], T_grid[i + 1]});
        }
    }

    for (auto &br : brackets) {
        double root = brent(f, context, br.a, br.b, tol);
        if (root > highest_root) {
            highest_root = root;
            found_exact = true;
        }
    }

    if (!found_exact) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return highest_root;
}

} // namespace supermag
