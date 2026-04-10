(** CLI entry point for parameter sweep orchestration.
    Reads sweep parameters, dispatches to C++ solvers via FFI,
    and aggregates results.

    Usage: sweep_driver --solver proximity --param d_F --range 0.5,20.0,50

    TODO: Implement CLI argument parsing and sweep dispatch. *)

let () =
  print_endline "SUPERMag sweep driver — not yet implemented";
  failwith "not implemented"
