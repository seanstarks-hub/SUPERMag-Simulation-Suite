(** Parameter sweep logic: grid, random, and adaptive strategies.

    Generates parameter sets and dispatches solver evaluations.
    Supports parallelism via Domain.spawn (OCaml 5 stdlib). *)

type sweep_type = Grid | Random | Adaptive

type sweep_config = {
  param_name : string;
  min_val : float;
  max_val : float;
  n_points : int;
  sweep_type : sweep_type;
}

(* ── Sweep strategies ────────────────────────────────── *)

let grid_sweep (config : sweep_config) : float array =
  let n = max 1 config.n_points in
  Array.init n (fun i ->
    if n = 1 then config.min_val
    else
      config.min_val +. (config.max_val -. config.min_val)
      *. Float.of_int i /. Float.of_int (n - 1))

let random_sweep (config : sweep_config) : float array =
  Array.init config.n_points (fun _ ->
    config.min_val +. Random.float (config.max_val -. config.min_val))

let adaptive_sweep (config : sweep_config) (eval_fn : float -> float) : float array =
  (* Start with a coarse grid, refine intervals with high curvature *)
  let coarse_n = max 5 (config.n_points / 3) in
  let coarse = grid_sweep { config with n_points = coarse_n } in
  let coarse_vals = Array.map eval_fn coarse in
  (* Compute second finite difference as curvature estimate *)
  let curvatures = Array.init (max 0 (coarse_n - 2)) (fun i ->
    Float.abs (coarse_vals.(i+2) -. 2.0 *. coarse_vals.(i+1) +. coarse_vals.(i))
  ) in
  let total_curv = Array.fold_left (+.) 0.0 curvatures in
  if total_curv < 1e-15 then
    (* Uniform — just return regular grid *)
    grid_sweep config
  else begin
    let remaining = config.n_points - coarse_n in
    (* Distribute extra points proportional to curvature *)
    let extra_per_interval = Array.map (fun c ->
      int_of_float (Float.of_int remaining *. c /. total_curv)
    ) curvatures in
    (* Build refined array *)
    let points = ref (Array.to_list coarse) in
    for i = 0 to Array.length curvatures - 1 do
      let n_extra = extra_per_interval.(i) in
      let lo = coarse.(i) in
      let hi = coarse.(i + 2) in
      for j = 1 to n_extra do
        let frac = Float.of_int j /. Float.of_int (n_extra + 1) in
        points := (lo +. frac *. (hi -. lo)) :: !points
      done
    done;
    let arr = Array.of_list !points in
    Array.sort Float.compare arr;
    arr
  end

(* ── Sweep parameter dispatch ────────────────────────── *)

type sweep_param =
  | Sweep_d_F
  | Sweep_d_S
  | Sweep_gamma
  | Sweep_gamma_B
  | Sweep_E_ex
  | Sweep_xi_F
  | Sweep_Tc0

let sweep_param_of_string = function
  | "d_F" -> Ok Sweep_d_F
  | "d_S" -> Ok Sweep_d_S
  | "gamma" -> Ok Sweep_gamma
  | "gamma_B" -> Ok Sweep_gamma_B
  | "E_ex" -> Ok Sweep_E_ex
  | "xi_F" -> Ok Sweep_xi_F
  | "Tc0" -> Ok Sweep_Tc0
  | s -> Error (Printf.sprintf "unknown sweep parameter: %s" s)

let sweep_param_to_string = function
  | Sweep_d_F -> "d_F"
  | Sweep_d_S -> "d_S"
  | Sweep_gamma -> "gamma"
  | Sweep_gamma_B -> "gamma_B"
  | Sweep_E_ex -> "E_ex"
  | Sweep_xi_F -> "xi_F"
  | Sweep_Tc0 -> "Tc0"

open Supermag_types

(** Apply a swept parameter value to a proximity_params record. *)
let apply_param (p : Params.proximity_params) param value =
  match param with
  | Sweep_d_F -> p  (* d_F is in the batch array, not in params *)
  | Sweep_d_S -> { p with d_s = value }
  | Sweep_gamma -> { p with gamma = value }
  | Sweep_gamma_B -> { p with gamma_b = value }
  | Sweep_E_ex ->
    (* Auto-recompute xi_F = sqrt(hbar * D_F / E_ex) when sweeping E_ex *)
    let hbar = Supermag_ffi.Stubs.supermag_const_hbar () in
    let xi_f_new = Float.sqrt (hbar *. p.d_f_coeff /. value) *. 1e9 in
    { p with e_ex = value; xi_f = xi_f_new }
  | Sweep_xi_F -> { p with xi_f = value }
  | Sweep_Tc0 -> { p with tc0 = value }

(** 1D parameter sweep for Tc vs. one parameter.

    When sweeping d_F, the sweep values ARE the d_F_array.
    For other parameters, [d_f_array] must be provided. *)
let tc_parameter_sweep param sweep_values ~(params : Params.proximity_params)
    ?(d_f_array = [|1.0|]) ?(depairing = Params.no_depairing) () =
  let config_name = sweep_param_to_string param in
  match param with
  | Sweep_d_F ->
    (* sweep_values IS the d_F array *)
    Supermag_ffi.Solvers.solve_tc ~params ~d_f_array:sweep_values ~depairing ()
  | _ ->
    (* For each swept value, do a full batch Tc for the d_f_array *)
    let n_sweep = Array.length sweep_values in
    let n_df = Array.length d_f_array in
    let all_tc = Array.init n_sweep (fun i ->
      let p = apply_param params param sweep_values.(i) in
      match Supermag_ffi.Solvers.solve_tc ~params:p ~d_f_array ~depairing () with
      | Ok r -> r.tc_values
      | Error msg ->
        Printf.eprintf "Warning: solver error at %s=%.6g: %s\n%!" config_name sweep_values.(i) msg;
        Array.make n_df Float.nan
    ) in
    (* Return the first d_F's Tc curve as function of swept param *)
    let tc_vals = Array.init n_sweep (fun i -> all_tc.(i).(0)) in
    Ok Result.{ d_f_values = sweep_values; tc_values = tc_vals; tc0 = params.tc0 }

(** 2D phase diagram: Tc as function of two parameters.

    Returns a 2D array [n1 x n2] of Tc values. *)
let tc_phase_diagram param1 values1 param2 values2
    ~(params : Params.proximity_params) ?(d_f_value = 1.0)
    ?(depairing = Params.no_depairing) () =
  let n1 = Array.length values1 in
  let n2 = Array.length values2 in
  Array.init n1 (fun i ->
    let p1 = apply_param params param1 values1.(i) in
    Array.init n2 (fun j ->
      let p2 = apply_param p1 param2 values2.(j) in
      let d_f_arr =
        if param1 = Sweep_d_F then [|values1.(i)|]
        else if param2 = Sweep_d_F then [|values2.(j)|]
        else [|d_f_value|]
      in
      match Supermag_ffi.Solvers.solve_tc ~params:p2 ~d_f_array:d_f_arr ~depairing () with
      | Ok r -> r.tc_values.(0)
      | Error msg ->
        Printf.eprintf "Warning: phase diagram solver error at (%s=%.6g, %s=%.6g): %s\n%!"
          (sweep_param_to_string param1) values1.(i)
          (sweep_param_to_string param2) values2.(j) msg;
        Float.nan
    )
  )
