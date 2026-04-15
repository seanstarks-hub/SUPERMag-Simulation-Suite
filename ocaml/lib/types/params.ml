(** Typed parameter records and enums for SUPERMag solvers.

    Provides compile-time phase/model checking instead of raw int codes. *)

type phase = Zero | Pi
type model = Thin_S | Fominov | Fominov_Multi

(** Geometry selection (orthogonal to equation model) *)
type geometry = Bilayer | Trilayer | Graded | Domains

(** Usadel solver mode *)
type usadel_mode = Linearized | Nonlinear

(** GL solver mode *)
type gl_mode = Scalar | Gauge

(** Triplet solver mode *)
type triplet_mode = Phenomenological | Diffusive

(** Graded ferromagnet exchange profile *)
type grade_profile = Linear | Exponential | Step

(** Geometry-specific parameters *)
type trilayer_params = {
  d_n : float;      (** Normal-metal interlayer thickness (nm) *)
  xi_n : float;     (** Normal-metal coherence length (nm) *)
  r_b : float;      (** Interface barrier resistance (Ohm·nm^2) *)
}

type graded_params = {
  e_ex_surface : float;    (** Exchange energy at S/F interface (meV) *)
  e_ex_bulk : float;       (** Exchange energy in F bulk (meV) *)
  profile : grade_profile;  (** Spatial profile shape *)
  n_slabs : int;           (** Number of discretization slabs *)
}

type domain_params = {
  domain_width : float;  (** Magnetic domain width (nm) *)
  domain_wall : float;   (** Domain wall thickness (nm), 0 = sharp *)
}

(** Spin-active interface parameters *)
type spin_active_params = {
  mixing_angle : float;   (** Spin-mixing angle (rad) *)
  polarization : float;   (** Interface spin polarization [0,1] *)
}

let phase_to_int = function Zero -> 0 | Pi -> 1
let model_to_int = function Thin_S -> 0 | Fominov -> 1 | Fominov_Multi -> 2

let geometry_to_int = function
  | Bilayer -> 0 | Trilayer -> 1 | Graded -> 2 | Domains -> 3

let usadel_mode_to_int = function Linearized -> 0 | Nonlinear -> 1
let gl_mode_to_int = function Scalar -> 0 | Gauge -> 1
let triplet_mode_to_int = function Phenomenological -> 0 | Diffusive -> 1
let grade_profile_to_int = function Linear -> 0 | Exponential -> 1 | Step -> 2

let phase_of_string = function
  | "zero" -> Ok Zero
  | "pi" -> Ok Pi
  | s -> Error (Printf.sprintf "unknown phase: %s" s)

let model_of_string = function
  | "thin_s" -> Ok Thin_S
  | "fominov" -> Ok Fominov
  | "fominov_multi" -> Ok Fominov_Multi
  | s -> Error (Printf.sprintf "unknown model: %s" s)

let geometry_of_string = function
  | "bilayer" -> Ok Bilayer
  | "trilayer" -> Ok Trilayer
  | "graded" -> Ok Graded
  | "domains" -> Ok Domains
  | s -> Error (Printf.sprintf "unknown geometry: %s" s)

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
  d_s_coeff : float;  (** Diffusion coefficient in S layer (nm^2/ps) *)
  xi_s : float;
  xi_f : float;
  gamma : float;
  gamma_b : float;
  e_ex : float;
  d_f_coeff : float;  (** Diffusion coefficient in F layer (nm^2/ps) *)
  model : model;
  phase : phase;
  geometry : geometry;
  geom_config : geom_config option;
  spin_active : spin_active_params option;
}

and geom_config =
  | Trilayer_config of trilayer_params
  | Graded_config of graded_params
  | Domain_config of domain_params

(** Fit flags: which parameters to include in least-squares fitting *)
type fit_flags = {
  fit_gamma : bool;
  fit_gamma_b : bool;
  fit_e_ex : bool;
  fit_xi_f : bool;
}

let default_fit_flags = {
  fit_gamma = true; fit_gamma_b = false;
  fit_e_ex = false; fit_xi_f = false;
}

(** Physical inputs for depairing computation *)
type depairing_input = {
  gamma_s_mev : float;   (** Spin-flip scattering rate (meV) *)
  h_tesla : float;       (** Applied magnetic field (T) *)
  d_nm2ps : float;       (** Diffusion coefficient (nm^2/ps) *)
  thickness_nm : float;  (** Film thickness (nm) *)
  gamma_so_mev : float;  (** Spin-orbit scattering rate (meV) *)
  t_kelvin : float;      (** Temperature (K) *)
}
