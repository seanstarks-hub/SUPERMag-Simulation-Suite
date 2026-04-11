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
