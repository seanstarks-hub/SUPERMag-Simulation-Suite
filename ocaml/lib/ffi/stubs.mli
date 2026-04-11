(** Signature for FFI stubs — ctypes bindings to C++ solvers. *)

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
  model:int -> phase:int ->
  depairing:(float * float * float * float) -> float

(** Solve for Tc across a batch of d_F values.
    [model]: 0 = thin-S, 1 = Fominov.
    [phase]: 0 = zero, 1 = pi.
    [depairing]: (ag, zeeman, orbital, spin_orbit) tuple. *)
val solve_tc_batch :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f_coeff:float ->
  model:int -> phase:int ->
  depairing:(float * float * float * float) ->
  d_f_arr:float array -> float array

(** Usadel diffusive-limit solver.
    Returns (Delta_out, x_out) arrays of length [n_grid]. *)
val usadel_solve :
  tc0:float -> d_s:float -> d_f:float ->
  xi_s:float -> xi_f:float -> e_ex:float ->
  n_grid:int -> float array * float array

(** Eilenberger clean-limit solver.
    Returns (f_out, x_out) arrays of length [n_grid]. *)
val eilenberger_solve :
  tc0:float -> d_s:float -> d_f:float ->
  xi_s:float -> e_ex:float ->
  n_grid:int -> float array * float array

(** BdG tight-binding Hamiltonian diagonalization.
    Returns eigenvalues array (sorted, length = 2*n_sites). *)
val bdg_solve :
  n_sites:int -> t_hop:float -> delta:float -> e_ex:float ->
  float array

(** Ginzburg-Landau free energy functional minimization.
    Returns (psi_real, psi_imag) arrays of length [nx*ny]. *)
val gl_minimize :
  alpha:float -> beta:float -> kappa:float ->
  nx:int -> ny:int -> dx:float ->
  float array * float array

(** Josephson current-phase relation for S/F/S junctions.
    Returns (phase_arr, current_out) arrays of length [n_phases]. *)
val josephson_cpr :
  d_f:float -> xi_f:float -> e_ex:float -> t:float ->
  n_phases:int -> float array * float array

(** Spin-triplet superconductivity solver.
    Returns (f_triplet_out, x_out) arrays of length [n_grid]. *)
val triplet_solve :
  n_layers:int -> thicknesses:float array ->
  magnetization_angles:float array ->
  n_grid:int -> float array * float array
