#ifndef SUPERMAG_ROOT_SCALAR_H
#define SUPERMAG_ROOT_SCALAR_H

namespace supermag {

// Find the highest root of f(T, context) = 0 in [T_min, T_max].
// Uses coarse grid scan (N=1000) to bracket sign changes,
// then refines each bracket with Brent's method to tolerance tol.
// Returns the highest root, or NaN if no root found.
double root_scalar_solve(
    double (*f)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol
);

} // namespace supermag

#endif
