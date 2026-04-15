#ifndef SUPERMAG_SOLVER_OPTIONS_H
#define SUPERMAG_SOLVER_OPTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Composable solver options for controlling numerical parameters.
 *
 * Pass NULL to any solver accepting this struct to use built-in defaults.
 * Call supermag_default_solver_options() to get a struct pre-filled with
 * defaults, then modify individual fields as needed.
 *
 * Which fields each solver reads:
 *   josephson:    matsubara_max, omega_cut_factor
 *   usadel:       matsubara_max, omega_cut_factor, max_iter, conv_tol
 *   eilenberger:  matsubara_max, omega_cut_factor
 *   ginzburg_landau: max_steps, conv_tol
 *   root_scalar:  root_grid_points
 */
typedef struct {
    int    matsubara_max;       /* Max Matsubara frequencies.          Default: 500   */
    double omega_cut_factor;   /* omega_cut = factor · Δ.             Default: 20.0  */
    int    max_steps;          /* Max relaxation steps (GL).          Default: 5000  */
    int    max_iter;           /* Max self-consistency iterations.    Default: 100   */
    double conv_tol;           /* Convergence tolerance.              Default: 1e-8  */
    int    root_grid_points;   /* Root-finder grid density.           Default: 1000  */
} supermag_solver_options_t;

/* Return a struct with all fields set to their defaults. */
supermag_solver_options_t supermag_default_solver_options(void);

#ifdef __cplusplus
}
#endif

#endif
