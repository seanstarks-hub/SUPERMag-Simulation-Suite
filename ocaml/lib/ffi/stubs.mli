(** Signature for FFI stubs — ctypes bindings to C++ solvers. *)

val supermag_const_hbar : unit -> float
val supermag_const_kB : unit -> float

(** Compute pair amplitude F(x) in the ferromagnet layer.
    [phase]: 0 = zero-phase (coth kernel), 1 = pi-phase (tanh kernel). *)
val pair_amplitude :
  xi_f:float -> d_f:float -> phase:int -> n_points:int ->
  float array * float array

(** Compute Tc for a batch of d_F values using the digamma self-consistency equation.
    [model]: 0 = thin-S, 1 = Fominov.
    [phase]: 0 = zero, 1 = pi.
    [depairing]: (ag, zeeman, orbital, spin_orbit) tuple. *)
val solve_tc_batch :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float ->
  gamma:float -> gamma_b:float -> e_ex:float -> d_f:float ->
  model:int -> phase:int ->
  depairing:(float * float * float * float) ->
  d_f_arr:float array -> float array

(** Sum depairing channels. *)
val depairing_total :
  ag:float -> zeeman:float -> orbital:float -> spin_orbit:float -> float
