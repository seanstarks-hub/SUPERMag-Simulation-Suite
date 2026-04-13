(** Ctypes bindings to C headers in cpp/include/supermag/.

    Uses OCaml ctypes to call into the compiled C++ library.
    The C headers define the FFI boundary — only C types cross it.

    Binds all 10 C API functions:
      constants (2), proximity (4), solvers (6). *)

open Ctypes
open Foreign

(* ── Error codes ─────────────────────────────────────── *)

let supermag_ok = 0

let check_rc rc =
  if rc <> supermag_ok then
    failwith (Printf.sprintf "SUPERMag C error code: %d" rc)

(* ── Constants ───────────────────────────────────────── *)

let c_const_hbar =
  foreign "supermag_const_hbar" (void @-> returning double)

let c_const_kB =
  foreign "supermag_const_kB" (void @-> returning double)

let supermag_const_hbar () = c_const_hbar ()
let supermag_const_kB () = c_const_kB ()

(* ── Depairing struct ────────────────────────────────── *)

type depairing_t
let depairing_struct : depairing_t structure typ =
  structure "supermag_depairing_t"
let dp_ag         = field depairing_struct "ag" double
let dp_zeeman     = field depairing_struct "zeeman" double
let dp_orbital    = field depairing_struct "orbital" double
let dp_spin_orbit = field depairing_struct "spin_orbit" double
let () = seal depairing_struct

(* ── Proximity params struct ─────────────────────────── *)

type params_t
let params_struct : params_t structure typ =
  structure "supermag_proximity_params_t"
let p_Tc0       = field params_struct "Tc0" double
let p_d_S       = field params_struct "d_S" double
let p_d_F       = field params_struct "d_F" double
let p_xi_S      = field params_struct "xi_S" double
let p_xi_F      = field params_struct "xi_F" double
let p_gamma     = field params_struct "gamma" double
let p_gamma_B   = field params_struct "gamma_B" double
let p_E_ex      = field params_struct "E_ex" double
let p_D_F       = field params_struct "D_F" double
let p_D_S       = field params_struct "D_S" double
let p_model     = field params_struct "model" int
let p_phase     = field params_struct "phase" int
let p_geometry  = field params_struct "geometry" int
let p_geom_params  = field params_struct "geom_params" (ptr void)
let p_spin_active  = field params_struct "spin_active" (ptr void)
let () = seal params_struct

(* ── C function bindings ─────────────────────────────── *)

let c_depairing_total =
  foreign "supermag_depairing_total"
    (ptr depairing_struct @-> returning double)

let c_solve_tc =
  foreign "supermag_proximity_solve_tc"
    (ptr params_struct @-> ptr depairing_struct @-> ptr double @-> returning int)

let c_solve_tc_batch =
  foreign "supermag_proximity_solve_tc_batch"
    (ptr params_struct @-> ptr double @-> int @-> ptr depairing_struct
     @-> ptr double @-> returning int)

let c_pair_amplitude =
  foreign "supermag_proximity_pair_amplitude"
    (double @-> double @-> int @-> int @-> ptr double @-> ptr double
     @-> returning int)

let c_usadel_solve =
  foreign "supermag_usadel_solve"
    (double @-> double @-> double @-> double @-> double @-> double
     @-> double @-> int @-> int @-> ptr double @-> ptr double @-> returning int)

let c_eilenberger_solve =
  foreign "supermag_eilenberger_solve"
    (double @-> double @-> double @-> double @-> double
     @-> double @-> int @-> ptr double @-> ptr double @-> returning int)

let c_bdg_solve =
  foreign "supermag_bdg_solve"
    (int @-> double @-> double @-> double @-> double
     @-> ptr double @-> ptr int @-> ptr double @-> returning int)

let c_gl_minimize =
  foreign "supermag_gl_minimize"
    (double @-> double @-> double @-> int @-> int @-> double
     @-> int @-> double @-> ptr double @-> ptr double @-> returning int)

let c_josephson_cpr =
  foreign "supermag_josephson_cpr"
    (double @-> double @-> double @-> double @-> double
     @-> double @-> int @-> ptr double @-> ptr double
     @-> ptr double @-> returning int)

let c_triplet_solve =
  foreign "supermag_triplet_solve"
    (int @-> ptr double @-> ptr double @-> ptr double @-> ptr double
     @-> double @-> double @-> double @-> int
     @-> int @-> ptr double @-> ptr double @-> returning int)

(* ── Helper: CArray to OCaml array ──────────────────── *)

let carray_to_array n ca =
  Array.init n (fun i -> CArray.get ca i)

(* ── OCaml wrappers ──────────────────────────────────── *)

let depairing_total ~ag ~zeeman ~orbital ~spin_orbit =
  let dp = make depairing_struct in
  setf dp dp_ag ag;
  setf dp dp_zeeman zeeman;
  setf dp dp_orbital orbital;
  setf dp dp_spin_orbit spin_orbit;
  c_depairing_total (addr dp)

let pair_amplitude ~d_f ~xi_f ~phase ~n_points =
  let x_out = CArray.make double n_points in
  let f_out = CArray.make double n_points in
  let rc = c_pair_amplitude d_f xi_f phase n_points
      (CArray.start x_out) (CArray.start f_out) in
  check_rc rc;
  (carray_to_array n_points x_out, carray_to_array n_points f_out)

let make_params ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~gamma ~gamma_b ~e_ex ~d_f_coeff
    ~d_s_coeff ~model ~phase ~geometry =
  let p = make params_struct in
  setf p p_Tc0 tc0;
  setf p p_d_S d_s;
  setf p p_d_F d_f;
  setf p p_xi_S xi_s;
  setf p p_xi_F xi_f;
  setf p p_gamma gamma;
  setf p p_gamma_B gamma_b;
  setf p p_E_ex e_ex;
  setf p p_D_F d_f_coeff;
  setf p p_D_S d_s_coeff;
  setf p p_model model;
  setf p p_phase phase;
  setf p p_geometry geometry;
  setf p p_geom_params (Ctypes.null);
  setf p p_spin_active (Ctypes.null);
  p

let make_depairing (ag, zeeman, orbital, spin_orbit) =
  let dp = make depairing_struct in
  setf dp dp_ag ag;
  setf dp dp_zeeman zeeman;
  setf dp dp_orbital orbital;
  setf dp dp_spin_orbit spin_orbit;
  dp

let solve_tc ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~gamma ~gamma_b ~e_ex ~d_f_coeff
    ~d_s_coeff ~model ~phase ~geometry ~depairing =
  let p = make_params ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~gamma ~gamma_b ~e_ex
      ~d_f_coeff ~d_s_coeff ~model ~phase ~geometry in
  let dp = make_depairing depairing in
  let tc_out = allocate double 0.0 in
  let rc = c_solve_tc (addr p) (addr dp) tc_out in
  check_rc rc;
  !@tc_out

let solve_tc_batch ~tc0 ~d_s ~xi_s ~xi_f ~gamma ~gamma_b ~e_ex ~d_f_coeff
    ~d_s_coeff ~model ~phase ~geometry ~depairing ~d_f_arr =
  let n = Array.length d_f_arr in
  let p = make_params ~tc0 ~d_s ~d_f:0.0 ~xi_s ~xi_f ~gamma ~gamma_b ~e_ex
      ~d_f_coeff ~d_s_coeff ~model ~phase ~geometry in
  let dp = make_depairing depairing in
  let c_df = CArray.make double n in
  Array.iteri (fun i v -> CArray.set c_df i v) d_f_arr;
  let tc_out = CArray.make double n in
  let rc = c_solve_tc_batch (addr p) (CArray.start c_df) n
      (addr dp) (CArray.start tc_out) in
  check_rc rc;
  carray_to_array n tc_out

let usadel_solve ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~t ~mode ~n_grid =
  let delta_out = CArray.make double n_grid in
  let x_out = CArray.make double n_grid in
  let rc = c_usadel_solve tc0 d_s d_f xi_s xi_f e_ex t mode n_grid
      (CArray.start delta_out) (CArray.start x_out) in
  check_rc rc;
  (carray_to_array n_grid delta_out, carray_to_array n_grid x_out)

let eilenberger_solve ~tc0 ~d_s ~d_f ~xi_s ~e_ex ~t ~n_grid =
  let f_out = CArray.make double n_grid in
  let x_out = CArray.make double n_grid in
  let rc = c_eilenberger_solve tc0 d_s d_f xi_s e_ex t n_grid
      (CArray.start f_out) (CArray.start x_out) in
  check_rc rc;
  (carray_to_array n_grid f_out, carray_to_array n_grid x_out)

let bdg_solve ~n_sites ~t_hop ~delta ~e_ex ?(mu = 0.0) () =
  let max_eig = 2 * n_sites in
  let eig_out = CArray.make double max_eig in
  let n_eig = allocate int 0 in
  let rc = c_bdg_solve n_sites t_hop delta e_ex mu
      (CArray.start eig_out) n_eig (from_voidp double null) in
  check_rc rc;
  carray_to_array (!@n_eig) eig_out

let gl_minimize ~alpha ~beta ~kappa ~nx ~ny ~dx ~mode ~h_applied =
  let n = nx * ny in
  let psi_real = CArray.make double n in
  let psi_imag = CArray.make double n in
  let rc = c_gl_minimize alpha beta kappa nx ny dx mode h_applied
      (CArray.start psi_real) (CArray.start psi_imag) in
  check_rc rc;
  (carray_to_array n psi_real, carray_to_array n psi_imag)

let josephson_cpr ~d_f ~xi_f ~e_ex ~t ?(tc0 = 9.2) ~gamma_b ~n_phases =
  let pi = Float.pi in
  let phase_arr = CArray.make double n_phases in
  for i = 0 to n_phases - 1 do
    CArray.set phase_arr i (2.0 *. pi *. Float.of_int i /. Float.of_int n_phases)
  done;
  let current_out = CArray.make double n_phases in
  let rc = c_josephson_cpr d_f xi_f e_ex t tc0 gamma_b n_phases
      (CArray.start phase_arr) (CArray.start current_out)
      (from_voidp double null) in
  check_rc rc;
  (carray_to_array n_phases phase_arr, carray_to_array n_phases current_out)

let triplet_solve ~n_layers ~thicknesses ~magnetization_angles
    ?(xi_f = 1.0) ?(xi_n = 10.0) ~t ~mode ~n_grid =
  let c_thick = CArray.make double n_layers in
  Array.iteri (fun i v -> CArray.set c_thick i v) thicknesses;
  let c_angles = CArray.make double n_layers in
  Array.iteri (fun i v -> CArray.set c_angles i v) magnetization_angles;
  let f_out = CArray.make double n_grid in
  let x_out = CArray.make double n_grid in
  let rc = c_triplet_solve n_layers
      (CArray.start c_thick) (CArray.start c_angles)
      (from_voidp double null) (from_voidp double null)
      xi_f xi_n t mode
      n_grid (CArray.start f_out) (CArray.start x_out) in
  check_rc rc;
  (carray_to_array n_grid f_out, carray_to_array n_grid x_out)
