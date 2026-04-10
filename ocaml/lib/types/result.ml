(** Result type: order parameter profiles, Tc values, free energies. *)

type profile = {
  x : float array;     (** Position (nm) *)
  values : float array; (** Order parameter or pair amplitude *)
}

type tc_result = {
  d_f_values : float array;  (** F-layer thicknesses (nm) *)
  tc_values : float array;   (** Corresponding Tc (K) *)
  tc0 : float;               (** Bulk Tc (K) *)
}

type solver_result =
  | Profile of profile
  | TcCurve of tc_result
  | Error of string
