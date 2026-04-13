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

type optimization_result = {
  d_f_optimal : float;  (** Optimal d_F thickness (nm) *)
}

type fit_result = {
  chi2 : float;  (** Final chi-squared residual *)
}

type solver_result =
  | Profile of profile
  | TcCurve of tc_result
  | Eigenvalues of eigenvalue_result
  | OrderParameter2D of order_parameter_2d
  | CurrentPhase of cpr_result
  | Optimization of optimization_result
  | Fit of fit_result
  | Error of string
