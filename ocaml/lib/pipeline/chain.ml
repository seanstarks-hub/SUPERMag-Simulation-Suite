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
