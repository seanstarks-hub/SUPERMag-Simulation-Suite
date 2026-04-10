(** Module signature for the pipeline interface. *)

type sweep_type = Grid | Random | Adaptive

type sweep_config = {
  param_name : string;
  min_val : float;
  max_val : float;
  n_points : int;
  sweep_type : sweep_type;
}

val grid_sweep : sweep_config -> float array
val random_sweep : sweep_config -> float array

type 'a solver_step = {
  name : string;
  run : 'a -> 'a;
}

val chain : 'a solver_step list -> 'a -> 'a
