#include "supermag/proximity.h"
#include <cmath>

extern "C" {

int supermag_proximity_pair_amplitude(
    double F0, double xi_F, double d_F,
    int n_points, double* x_out, double* F_out)
{
    if (!x_out || !F_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_points < 2 || d_F <= 0.0 || xi_F <= 0.0)
        return SUPERMAG_ERR_INVALID_DIM;

    double dx = d_F / (n_points - 1);
    for (int i = 0; i < n_points; ++i) {
        double x = i * dx;
        x_out[i] = x;
        // FFLO-like decaying oscillation in the ferromagnet
        // F(x) = F0 * exp(-x/xi_F) * cos(x/xi_F)
        F_out[i] = F0 * std::exp(-x / xi_F) * std::cos(x / xi_F);
    }
    return SUPERMAG_OK;
}

}
