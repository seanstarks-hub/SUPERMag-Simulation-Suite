(** Signature for FFI stubs — ctypes bindings to C++ solvers. *)

val supermag_const_hbar : unit -> float
val supermag_const_kB : unit -> float

val pair_amplitude :
  f0:float -> xi_f:float -> d_f:float -> n_points:int ->
  float array * float array

val critical_temp :
  tc0:float -> d_s:float -> xi_s:float -> xi_f:float -> e_ex:float ->
  d_f_arr:float array -> float array
