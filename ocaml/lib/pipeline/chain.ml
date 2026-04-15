(** Solver chaining: compose sequential solver runs.

    Allows piping the output of one solver into the input of another.
    E.g., compute Delta(x) from Usadel, then feed into Josephson CPR. *)

type 'a solver_step = {
  name : string;
  run : 'a -> 'a;
}

let chain (steps : 'a solver_step list) (input : 'a) : 'a =
  List.fold_left (fun acc step -> step.run acc) input steps

let parallel_chain (steps : 'a solver_step list) (inputs : 'a list) : 'a list =
  let pairs = List.combine steps inputs in
  let domains = List.map (fun (step, input) ->
    Domain.spawn (fun () -> step.run input)
  ) pairs in
  List.map Domain.join domains

(* ── Result-monad pipeline ───────────────────────────── *)

open Supermag_types

type solver_step_r = {
  name : string;
  run_r : Result.solver_result -> (Result.solver_result, string) result;
}

(** Chain solver steps with error short-circuiting.
    On first error, the remaining steps are skipped. *)
let solver_chain (steps : solver_step_r list)
    (input : Result.solver_result) : (Result.solver_result, string) result =
  List.fold_left (fun acc step ->
    match acc with
    | Error _ -> acc
    | Ok v -> step.run_r v
  ) (Ok input) steps

(** Run multiple independent solver steps in parallel,
    collecting results or errors. *)
let parallel_solver_chain (steps : solver_step_r list)
    (inputs : Result.solver_result list)
    : (Result.solver_result, string) result list =
  let pairs = List.combine steps inputs in
  let domains = List.map (fun (step, input) ->
    Domain.spawn (fun () -> step.run_r input)
  ) pairs in
  List.map Domain.join domains

(* ── Typed adapter functions ─────────────────────────── *)

open Supermag_ffi

(** Find the index of the minimum value in a float array.
    Raises [Invalid_argument] for empty arrays. *)
let argmin_float arr =
  let n = Array.length arr in
  if n = 0 then invalid_arg "argmin_float: empty array"
  else
    let idx = ref 0 in
    for i = 1 to n - 1 do
      if arr.(i) < arr.(!idx) then idx := i
    done;
    !idx

(** Typed adapter: from a Tc curve, find the d_F at minimum Tc
    and return an optimization result. *)
let tc_to_minimum (tc : Result.tc_result)
    : (Result.optimization_result, string) result =
  let n = Array.length tc.tc_values in
  if n = 0 then Error "tc_result has no data points"
  else
    let i = argmin_float tc.tc_values in
    Ok Result.{ d_f_optimal = tc.d_f_values.(i) }

(** Typed adapter: compute Tc curve, then run Usadel at the d_F where
    Tc is minimized to extract the spatial order parameter profile.
    Returns (tc_result, profile) pair. *)
let tc_then_usadel ~(params : Params.proximity_params) ~d_f_array ~t
    ?(depairing = Params.no_depairing) ?(n_grid = 100) ()
    : (Result.tc_result * Result.profile, string) result =
  if Array.length d_f_array = 0 then Error "d_f_array must not be empty"
  else
    match Solvers.solve_tc ~params ~d_f_array ~depairing () with
    | Error e -> Error e
    | Ok tc ->
      let i = argmin_float tc.tc_values in
      let d_f_min = tc.d_f_values.(i) in
      match Solvers.usadel ~tc0:params.tc0 ~d_s:params.d_s ~d_f:d_f_min
          ~xi_s:params.xi_s ~xi_f:params.xi_f ~e_ex:params.e_ex
          ~t ~mode:Params.Linearized ~n_grid with
      | Error e -> Error e
      | Ok profile -> Ok (tc, profile)

(** Typed adapter: compute Delta(x) profile from Usadel, then evaluate
    Josephson CPR using the interface Delta value.
    Returns (profile, cpr_result) pair. *)
let usadel_then_josephson ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~t ~gamma_b
    ?(n_grid = 100) ?(n_phases = 200) ()
    : (Result.profile * Result.cpr_result, string) result =
  if d_f <= 0.0 then Error "d_f must be > 0"
  else
    match Solvers.usadel ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~t
        ~mode:Params.Nonlinear ~n_grid with
    | Error e -> Error e
    | Ok profile ->
      match Solvers.josephson ~d_f ~xi_f ~e_ex ~t ~gamma_b ~n_phases with
      | Error e -> Error e
      | Ok cpr -> Ok (profile, cpr)
