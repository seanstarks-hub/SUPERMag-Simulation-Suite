(** Ctypes bindings to C headers in cpp/include/supermag/.

    Uses OCaml ctypes to call into the compiled C++ library.
    The C headers define the FFI boundary — only C types cross it.

    TODO: Bind all solver functions from proximity.h, usadel.h, etc. *)

(* open Ctypes *)

(* Physical constants *)
let supermag_const_hbar () : float =
  (* TODO: Call supermag_const_hbar via ctypes *)
  failwith "not implemented"

let supermag_const_kB () : float =
  failwith "not implemented"

(* Proximity solver *)
let pair_amplitude ~f0:_ ~xi_f:_ ~d_f:_ ~n_points:_ : float array * float array =
  (* TODO: Call supermag_proximity_pair_amplitude via ctypes *)
  failwith "not implemented"

let critical_temp ~tc0:_ ~d_s:_ ~xi_s:_ ~xi_f:_ ~e_ex:_ ~d_f_arr:_ : float array =
  (* TODO: Call supermag_proximity_critical_temp via ctypes *)
  failwith "not implemented"
