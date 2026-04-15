(** Signature for FFI stubs — ctypes bindings to C++ solvers. *)

(** Solver options for controlling numerical parameters. *)
type solver_options = {
  matsubara_max : int;
  omega_cut_factor : float;
  max_steps : int;
  max_iter : int;
  conv_tol : float;
  root_grid_points : int;
}

(** Return default solver options matching C defaults. *)
val default_solver_options : unit -> solver_options

val supermag_const_hbar : unit -> float
val supermag_const_kB : unit -> float

(** Sum depairing channels. *)
val depairing_total :
  ag:float -> zeeman:float -> orbital:float -> spin_orbit:float -> float

(** Compute pair amplitude F(x) in the ferromagnet layer.
    [phase]: 0 = zero-phase (coth kernel), 1 = pi-phase (tanh kernel). *)
val pair_amplitude :
  d_f:float -> xi_f:float -> phase:int -> n_points:int ->
  float array * float array

(** Solve for Tc at a single d_F value. *)
val solve_tc :
  tc0:float -> d_s:float -> d_f:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  d_s_coeff:float -> model:int -> phase:int -> geometry:int ->
  depairing:(float * float * float * float) -> float

(** Solve for Tc across a batch of d_F values.
    [model]: 0 = thin-S, 1 = Fominov, 2 = Fominov multimode.
    [phase]: 0 = zero, 1 = pi.
    [geometry]: 0 = bilayer, 1 = trilayer, 2 = graded, 3 = domains.
    [depairing]: (ag, zeeman, orbital, spin_orbit) tuple. *)
val solve_tc_batch :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  d_s_coeff:float -> model:int -> phase:int -> geometry:int ->
  depairing:(float * float * float * float) ->
  d_f_arr:float array -> float array

(** Usadel diffusive-limit solver.
    [t]: temperature (K), must be > 0.
    [mode]: 0 = linearized, 1 = nonlinear.
    Returns (Delta_out, x_out) arrays of length [n_grid]. *)
val usadel_solve :
  tc0:float -> d_s:float -> d_f:float ->
  xi_s:float -> xi_f:float -> e_ex:float ->
  t:float -> mode:int -> ?opts:solver_options -> n_grid:int -> float array * float array

(** Eilenberger clean-limit solver.
    [t]: temperature (K), must be > 0.
    Returns (f_out, x_out) arrays of length [n_grid]. *)
val eilenberger_solve :
  tc0:float -> d_s:float -> d_f:float ->
  xi_s:float -> e_ex:float ->
  t:float -> ?opts:solver_options -> n_grid:int -> float array * float array

(** BdG tight-binding Hamiltonian diagonalization.
    [mu]: chemical potential (meV); defaults to 0.0.
    Returns eigenvalues array (sorted, length = 2*n_sites). *)
val bdg_solve :
  n_sites:int -> t_hop:float -> delta:float -> e_ex:float ->
  ?mu:float -> unit -> float array

(** Ginzburg-Landau free energy functional minimization.
    [mode]: 0 = scalar, 1 = gauge.
    [h_applied]: applied field (GL units). Ignored in scalar mode.
    Returns (psi_real, psi_imag) arrays of length [nx*ny]. *)
val gl_minimize :
  alpha:float -> beta:float -> kappa:float ->
  nx:int -> ny:int -> dx:float ->
  mode:int -> ?opts:solver_options -> h_applied:float ->
  float array * float array

(** Josephson current-phase relation for S/F/S junctions.
    [tc0]: bulk Tc (K); defaults to 9.2 (Nb).
    [gamma_b]: interface barrier parameter (dimensionless).
    Returns (phase_out, current_out) arrays of length [n_phases]. *)
val josephson_cpr :
  d_f:float -> xi_f:float -> e_ex:float -> t:float ->
  ?tc0:float -> gamma_b:float -> ?opts:solver_options -> n_phases:int -> float array * float array

(** Spin-triplet superconductivity solver.
    [xi_f]: ferromagnetic coherence length (nm); defaults to 1.0.
    [xi_n]: triplet coherence length (nm); defaults to 10.0.
    [t]: temperature (K), must be > 0.
    [mode]: 0 = phenomenological, 1 = diffusive.
    Returns (f_triplet_out, x_out) arrays of length [n_grid]. *)
val triplet_solve :
  n_layers:int -> thicknesses:float array ->
  magnetization_angles:float array ->
  ?xi_f:float -> ?xi_n:float ->
  t:float -> ?tc0:float -> mode:int ->
  n_grid:int -> float array * float array

(** Individual depairing channel: Abrikosov-Gor'kov *)
val depairing_ag : gamma_s_mev:float -> t_kelvin:float -> float

(** Individual depairing channel: Zeeman *)
val depairing_zeeman : h_tesla:float -> t_kelvin:float -> float

(** Individual depairing channel: orbital (perpendicular field) *)
val depairing_orbital_perp :
  d_nm2ps:float -> h_tesla:float -> thickness_nm:float -> t_kelvin:float -> float

(** Individual depairing channel: orbital (parallel field) *)
val depairing_orbital_par :
  d_nm2ps:float -> h_tesla:float -> thickness_nm:float -> t_kelvin:float -> float

(** Individual depairing channel: spin-orbit coupling *)
val depairing_soc : gamma_so_mev:float -> t_kelvin:float -> float

(** Compute all depairing channels from physical inputs.
    Returns (ag, zeeman, orbital, spin_orbit) tuple. *)
val depairing_from_physical :
  gamma_s_mev:float -> h_tesla:float -> d_nm2ps:float ->
  thickness_nm:float -> gamma_so_mev:float -> t_kelvin:float ->
  float * float * float * float

(** Golden-section optimize d_F to match target Tc.
    Returns optimal d_F (nm). *)
val optimize_tc :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  d_s_coeff:float -> model:int -> phase:int -> geometry:int ->
  depairing:(float * float * float * float) ->
  d_f_lo:float -> d_f_hi:float -> tc_target:float -> float

(** Inverse solve: find d_F that gives target Tc via Brent's method.
    Returns d_F (nm). *)
val inverse_tc :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  d_s_coeff:float -> model:int -> phase:int -> geometry:int ->
  depairing:(float * float * float * float) ->
  tc_target:float -> d_f_lo:float -> d_f_hi:float -> float

(** Nelder-Mead least-squares fit.
    Returns (chi2, fitted_gamma, fitted_gamma_B, fitted_E_ex, fitted_xi_F). *)
val fit_tc :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  d_s_coeff:float -> model:int -> phase:int -> geometry:int ->
  depairing:(float * float * float * float) ->
  d_f_data:float array -> tc_data:float array ->
  fit_gamma:bool -> fit_gamma_b:bool -> fit_e_ex:bool -> fit_xi_f:bool ->
  float * float * float * float * float
