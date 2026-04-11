(** Result type: order parameter profiles, Tc values, eigenvalues, CPR. *)

type profile = {
  x : float array;     (** Position (nm) *)
  values : float array; (** Order parameter or pair amplitude *)
}

type tc_result = {
  d_f_values : float array;  (** F-layer thicknesses (nm) *)
  tc_values : float array;   (** Corresponding Tc (K) *)
  tc0 : float;               (** Bulk Tc (K) *)
}

type eigenvalue_result = {
  eigenvalues : float array;  (** meV, sorted *)
}

type order_parameter_2d = {
  nx : int;
  ny : int;
  psi_real : float array;
  psi_imag : float array;
}

type cpr_result = {
  phi : float array;      (** Phase (rad) *)
  current : float array;  (** Normalized current *)
}

type solver_result =
  | Profile of profile
  | TcCurve of tc_result
  | Eigenvalues of eigenvalue_result
  | OrderParameter2D of order_parameter_2d
  | CurrentPhase of cpr_result
  | Error of string
