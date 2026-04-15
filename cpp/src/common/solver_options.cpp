// Default solver options implementation.

#include "supermag/solver_options.h"

extern "C" {

supermag_solver_options_t supermag_default_solver_options(void) {
    supermag_solver_options_t opts;
    opts.matsubara_max     = 500;
    opts.omega_cut_factor  = 20.0;
    opts.max_steps         = 5000;
    opts.max_iter          = 100;
    opts.conv_tol          = 1e-8;
    opts.root_grid_points  = 1000;
    return opts;
}

}
