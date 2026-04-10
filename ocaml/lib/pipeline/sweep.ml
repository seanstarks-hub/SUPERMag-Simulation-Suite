(** Parameter sweep logic: grid, random, and adaptive strategies.

    Generates parameter sets and dispatches solver evaluations.
    Supports parallelism via Domainslib (OCaml 5).

    TODO: Implement grid_sweep, random_sweep, adaptive_sweep. *)

type sweep_type = Grid | Random | Adaptive

type sweep_config = {
  param_name : string;
  min_val : float;
  max_val : float;
  n_points : int;
  sweep_type : sweep_type;
}

let grid_sweep (_config : sweep_config) : float array =
  (* TODO: Generate uniform grid of parameter values *)
  failwith "not implemented"

let random_sweep (_config : sweep_config) : float array =
  (* TODO: Generate random parameter samples *)
  failwith "not implemented"

let adaptive_sweep (_config : sweep_config) (_eval_fn : float -> float) : float array =
  (* TODO: Adaptive refinement based on curvature *)
  failwith "not implemented"
