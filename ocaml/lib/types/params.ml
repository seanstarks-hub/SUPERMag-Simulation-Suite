(** Typed parameter records and enums for SUPERMag solvers.

    Provides compile-time phase/model checking instead of raw int codes. *)

type phase = Zero | Pi
type model = Thin_S | Fominov

let phase_to_int = function Zero -> 0 | Pi -> 1
let model_to_int = function Thin_S -> 0 | Fominov -> 1

let phase_of_string = function
  | "zero" -> Ok Zero
  | "pi" -> Ok Pi
  | s -> Error (Printf.sprintf "unknown phase: %s" s)

let model_of_string = function
  | "thin_s" -> Ok Thin_S
  | "fominov" -> Ok Fominov
  | s -> Error (Printf.sprintf "unknown model: %s" s)

type depairing = {
  ag : float;
  zeeman : float;
  orbital : float;
  spin_orbit : float;
}

let no_depairing = { ag = 0.0; zeeman = 0.0; orbital = 0.0; spin_orbit = 0.0 }

let depairing_to_tuple dp = (dp.ag, dp.zeeman, dp.orbital, dp.spin_orbit)

type proximity_params = {
  tc0 : float;
  d_s : float;
  xi_s : float;
  xi_f : float;
  gamma : float;
  gamma_b : float;
  e_ex : float;
  d_f_coeff : float;  (** Diffusion coefficient m^2/s *)
  model : model;
  phase : phase;
}
