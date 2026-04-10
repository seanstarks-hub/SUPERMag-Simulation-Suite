(** Solver chaining: compose sequential solver runs.

    Allows piping the output of one solver into the input of another.
    E.g., compute Delta(x) from Usadel, then feed into Josephson CPR.

    TODO: Implement chain combinator and result threading. *)

type 'a solver_step = {
  name : string;
  run : 'a -> 'a;
}

let chain (steps : 'a solver_step list) (input : 'a) : 'a =
  (* TODO: Execute steps sequentially, threading results *)
  ignore steps;
  ignore input;
  failwith "not implemented"

let parallel_chain (_steps : 'a solver_step list) (_inputs : 'a list) : 'a list =
  (* TODO: Run independent steps in parallel using Domainslib *)
  failwith "not implemented"
