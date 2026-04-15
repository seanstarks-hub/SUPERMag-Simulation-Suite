#ifndef SUPERMAG_ROOT_SCALAR_H
#define SUPERMAG_ROOT_SCALAR_H

namespace supermag {

// Find the highest root of f(T, context) = 0 in [T_min, T_max].
// Uses coarse grid scan to bracket sign changes,
// then refines each bracket with Brent's method to tolerance tol.
// grid_points: number of scan points (default 1000).
// Returns the highest root, or NaN if no root found.
double root_scalar_solve(
    double (*f)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol,
    int grid_points = 1000
);

// Log-spaced variant for problems where roots can span orders of magnitude.
// Requires T_min > 0. Uses log-spaced grid for better resolution at low T.
double root_scalar_solve_log(
    double (*f)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol,
    int grid_points = 1000
);

} // namespace supermag

#endif
