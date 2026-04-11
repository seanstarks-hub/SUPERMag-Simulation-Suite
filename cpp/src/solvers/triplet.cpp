// Spin-triplet superconductivity solver
// Long-range triplet correlations from inhomogeneous magnetization.
// Equal-spin triplets generated at non-collinear interfaces.

#include "supermag/triplet.h"
#include <cmath>
#include <algorithm>

extern "C" {

int supermag_triplet_solve(
    int n_layers, const double* thicknesses, const double* magnetization_angles,
    int n_grid, double* f_triplet_out, double* x_out)
{
    if (!thicknesses || !magnetization_angles || !f_triplet_out || !x_out)
        return SUPERMAG_ERR_NULL_PTR;
    if (n_layers < 2 || n_grid <= 0)
        return SUPERMAG_ERR_INVALID_DIM;

    // Hardcoded coherence lengths  [KNOWN-LIMIT-4]
    const double xi_F = 1.0;   // nm
    const double xi_N = 10.0;  // nm

    // Total thickness
    double total = 0.0;
    for (int i = 0; i < n_layers; ++i)
        total += thicknesses[i];

    // Build spatial grid
    for (int i = 0; i < n_grid; ++i)
        x_out[i] = total * i / std::max(n_grid - 1, 1);

    // Initialize
    for (int i = 0; i < n_grid; ++i)
        f_triplet_out[i] = 0.0;

    // Interface positions
    double cumulative = 0.0;
    for (int lay = 0; lay < n_layers - 1; ++lay) {
        cumulative += thicknesses[lay];
        double x_int = cumulative;

        // Misalignment angle
        double alpha = magnetization_angles[lay + 1] - magnetization_angles[lay];
        double triplet_amp = std::fabs(std::sin(alpha));

        // Long-range triplet contribution from this interface
        for (int j = 0; j < n_grid; ++j) {
            double dist = std::fabs(x_out[j] - x_int);
            f_triplet_out[j] += triplet_amp * std::exp(-dist / xi_N);
        }
    }

    return SUPERMAG_OK;
}

}
