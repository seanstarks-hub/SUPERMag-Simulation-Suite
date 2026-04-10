// Root-finding for determinant equations det M(T) = 0
// Same strategy as root_scalar: coarse grid scan → bracket → Brent's method.
// Used for the Fominov model where det M(T) = 0 determines Tc.

#include "root_scalar.h"

// The determinant solver uses the exact same algorithm as root_scalar_solve.
// The function pointer f(T, context) returns the determinant value.
// The caller is responsible for constructing the determinant function.
//
// We re-export root_scalar_solve under a determinant-specific name
// for API clarity, but the implementation is identical.

namespace supermag {

double determinant_solve(
    double (*det_func)(double, void*),
    void *context,
    double T_min,
    double T_max,
    double tol)
{
    return root_scalar_solve(det_func, context, T_min, T_max, tol);
}

} // namespace supermag
