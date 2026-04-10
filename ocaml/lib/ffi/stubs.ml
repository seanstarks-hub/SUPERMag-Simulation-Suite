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

(* Proximity solver — pair amplitude *)
let pair_amplitude ~xi_f:_ ~d_f:_ ~phase:_ ~n_points:_ : float array * float array =
  (* TODO: Call supermag_proximity_pair_amplitude(d_F, xi_F, phase, n_points, x_out, F_out) via ctypes *)
  failwith "not implemented"

(* Proximity solver — batch Tc calculation *)
let solve_tc_batch ~tc0:_ ~d_s:_ ~xi_s:_ ~xi_f:_
    ~gamma:_ ~gamma_b:_ ~e_ex:_ ~d_f:_
    ~model:_ ~phase:_
    ~depairing:_ ~d_f_arr:_ : float array =
  (* TODO: Call supermag_proximity_solve_tc_batch via ctypes *)
  failwith "not implemented"

(* Depairing — total pair-breaking parameter *)
let depairing_total ~ag:_ ~zeeman:_ ~orbital:_ ~spin_orbit:_ : float =
  (* TODO: Call supermag_depairing_total via ctypes *)
  failwith "not implemented"
