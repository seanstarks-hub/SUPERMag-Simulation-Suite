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

(** Result of evaluating a single device stack. *)
type device_result = {
  description : string;   (** Stack description, e.g. "Nb:30/Fe:8" *)
  tc : tc_result;         (** Tc vs d_F curve *)
  tc_min : float;         (** Minimum Tc achieved (K) *)
  d_f_at_min : float;     (** d_F at minimum Tc (nm) *)
}

(** Result of a combinatorial exploration. *)
type exploration_result = {
  ranked : device_result list;  (** Stacks ranked by figure of merit *)
  total_evaluated : int;        (** Total stacks evaluated *)
}

type solver_result =
  | Profile of profile
  | TcCurve of tc_result
  | Eigenvalues of eigenvalue_result
  | OrderParameter2D of order_parameter_2d
  | CurrentPhase of cpr_result
  | Optimization of optimization_result
  | Fit of fit_result
  | DeviceResult of device_result
  | ExplorationResult of exploration_result
  | Error of string
