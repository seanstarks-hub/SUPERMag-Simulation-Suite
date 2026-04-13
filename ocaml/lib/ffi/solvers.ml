(** Type-safe solver API.

    Calls [Stubs] internally. Every function validates inputs,
    converts typed enums to ints, and returns typed results. *)

open Supermag_types

let solve_tc ~(params : Params.proximity_params) ~d_f_array
    ?(depairing = Params.no_depairing) () =
  if params.tc0 <= 0.0 then
    Error "tc0 must be > 0"
  else if params.xi_s <= 0.0 then
    Error "xi_s must be > 0"
  else if params.xi_f <= 0.0 then
    Error "xi_f must be > 0"
  else if params.d_s <= 0.0 then
    Error "d_s must be > 0"
  else if Array.length d_f_array < 1 then
    Error "d_f_array must not be empty"
  else
    try
      let tc_arr = Stubs.solve_tc_batch
          ~tc0:params.tc0 ~d_s:params.d_s ~xi_s:params.xi_s ~xi_f:params.xi_f
          ~gamma:params.gamma ~gamma_b:params.gamma_b ~e_ex:params.e_ex
          ~d_f_coeff:params.d_f_coeff ~d_s_coeff:params.d_s_coeff
          ~model:(Params.model_to_int params.model)
          ~phase:(Params.phase_to_int params.phase)
          ~geometry:(Params.geometry_to_int params.geometry)
          ~depairing:(Params.depairing_to_tuple depairing)
          ~d_f_arr:d_f_array in
      Ok Result.{ d_f_values = d_f_array; tc_values = tc_arr; tc0 = params.tc0 }
    with Failure msg -> Error msg

let pair_amplitude ~d_f ~xi_f ~phase ~n_points =
  if d_f <= 0.0 then Error "d_f must be > 0"
  else if xi_f <= 0.0 then Error "xi_f must be > 0"
  else if n_points < 2 then Error "n_points must be >= 2"
  else
    try
      let (x, f) = Stubs.pair_amplitude
          ~d_f ~xi_f ~phase:(Params.phase_to_int phase) ~n_points in
      Ok Result.{ x; values = f }
    with Failure msg -> Error msg

let usadel ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~t ~mode ~n_grid =
  if n_grid < 5 then Error "n_grid must be >= 5"
  else if t <= 0.0 then Error "t must be > 0"
  else
    try
      let (delta, x) = Stubs.usadel_solve ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex
          ~t ~mode:(Params.usadel_mode_to_int mode) ~n_grid in
      Ok Result.{ x; values = delta }
    with Failure msg -> Error msg

let eilenberger ~tc0 ~d_s ~d_f ~xi_s ~e_ex ~t ~n_grid =
  if n_grid < 5 then Error "n_grid must be >= 5"
  else if t <= 0.0 then Error "t must be > 0"
  else
    try
      let (f, x) = Stubs.eilenberger_solve ~tc0 ~d_s ~d_f ~xi_s ~e_ex ~t ~n_grid in
      Ok Result.{ x; values = f }
    with Failure msg -> Error msg

let bdg ~n_sites ~t_hop ~delta ~e_ex ?mu () =
  if n_sites < 2 then Error "n_sites must be >= 2"
  else
    try
      let eigs = Stubs.bdg_solve ~n_sites ~t_hop ~delta ~e_ex ?mu () in
      Ok Result.{ eigenvalues = eigs }
    with Failure msg -> Error msg

let gl ~alpha ~beta ~kappa ~nx ~ny ~dx ~mode ~h_applied =
  if nx < 2 || ny < 2 then Error "nx and ny must be >= 2"
  else if dx <= 0.0 then Error "dx must be > 0"
  else
    try
      let (psi_r, psi_i) = Stubs.gl_minimize ~alpha ~beta ~kappa ~nx ~ny ~dx
          ~mode:(Params.gl_mode_to_int mode) ~h_applied in
      Ok Result.{ nx; ny; psi_real = psi_r; psi_imag = psi_i }
    with Failure msg -> Error msg

let josephson ~d_f ~xi_f ~e_ex ~t ~gamma_b ~n_phases =
  if n_phases < 2 then Error "n_phases must be >= 2"
  else if d_f <= 0.0 then Error "d_f must be > 0"
  else
    try
      let (phi, current) = Stubs.josephson_cpr ~d_f ~xi_f ~e_ex ~t ~gamma_b ~n_phases in
      Ok Result.{ phi; current }
    with Failure msg -> Error msg

let triplet ~n_layers ~thicknesses ~magnetization_angles ~t ~mode ~n_grid =
  if n_layers < 2 then Error "n_layers must be >= 2"
  else if Array.length thicknesses <> n_layers then
    Error "thicknesses length must equal n_layers"
  else if Array.length magnetization_angles <> n_layers then
    Error "magnetization_angles length must equal n_layers"
  else if n_grid < 5 then Error "n_grid must be >= 5"
  else if t <= 0.0 then Error "t must be > 0"
  else
    try
      let (f, x) = Stubs.triplet_solve
          ~n_layers ~thicknesses ~magnetization_angles
          ~t ~mode:(Params.triplet_mode_to_int mode) ~n_grid in
      Ok Result.{ x; values = f }
    with Failure msg -> Error msg

(* ── Depairing from physical inputs ──────────────────── *)

let depairing_from_physical ~gamma_s_mev ~h_tesla ~d_nm2ps
    ~thickness_nm ~gamma_so_mev ~t_kelvin =
  if t_kelvin <= 0.0 then Error "t_kelvin must be > 0"
  else
    try
      let (ag, zeeman, orbital, spin_orbit) =
        Stubs.depairing_from_physical ~gamma_s_mev ~h_tesla ~d_nm2ps
          ~thickness_nm ~gamma_so_mev ~t_kelvin in
      Ok Params.{ ag; zeeman; orbital; spin_orbit }
    with Failure msg -> Error msg

(* ── Optimizer wrappers ──────────────────────────────── *)

let optimize_tc ~(params : Params.proximity_params) ~d_f_lo ~d_f_hi
    ~tc_target ?(depairing = Params.no_depairing) () =
  if d_f_lo >= d_f_hi then Error "d_f_lo must be < d_f_hi"
  else if tc_target <= 0.0 then Error "tc_target must be > 0"
  else
    try
      let d_f_opt = Stubs.optimize_tc
          ~tc0:params.tc0 ~d_s:params.d_s ~xi_s:params.xi_s ~xi_f:params.xi_f
          ~gamma:params.gamma ~gamma_b:params.gamma_b ~e_ex:params.e_ex
          ~d_f_coeff:params.d_f_coeff ~d_s_coeff:params.d_s_coeff
          ~model:(Params.model_to_int params.model)
          ~phase:(Params.phase_to_int params.phase)
          ~geometry:(Params.geometry_to_int params.geometry)
          ~depairing:(Params.depairing_to_tuple depairing)
          ~d_f_lo ~d_f_hi ~tc_target in
      Ok Result.{ d_f_optimal = d_f_opt }
    with Failure msg -> Error msg

let inverse_tc ~(params : Params.proximity_params) ~tc_target ~d_f_lo ~d_f_hi
    ?(depairing = Params.no_depairing) () =
  if d_f_lo >= d_f_hi then Error "d_f_lo must be < d_f_hi"
  else if tc_target <= 0.0 then Error "tc_target must be > 0"
  else
    try
      let d_f = Stubs.inverse_tc
          ~tc0:params.tc0 ~d_s:params.d_s ~xi_s:params.xi_s ~xi_f:params.xi_f
          ~gamma:params.gamma ~gamma_b:params.gamma_b ~e_ex:params.e_ex
          ~d_f_coeff:params.d_f_coeff ~d_s_coeff:params.d_s_coeff
          ~model:(Params.model_to_int params.model)
          ~phase:(Params.phase_to_int params.phase)
          ~geometry:(Params.geometry_to_int params.geometry)
          ~depairing:(Params.depairing_to_tuple depairing)
          ~tc_target ~d_f_lo ~d_f_hi in
      Ok Result.{ d_f_optimal = d_f }
    with Failure msg -> Error msg

let fit_tc ~(params : Params.proximity_params) ~d_f_data ~tc_data
    ~(flags : Params.fit_flags) ?(depairing = Params.no_depairing) () =
  let n = Array.length d_f_data in
  if n < 2 then Error "need at least 2 data points"
  else if Array.length tc_data <> n then
    Error "d_f_data and tc_data must have same length"
  else
    try
      let (chi2, _gamma, _gamma_b, _e_ex, _xi_f) = Stubs.fit_tc
          ~tc0:params.tc0 ~d_s:params.d_s ~xi_s:params.xi_s ~xi_f:params.xi_f
          ~gamma:params.gamma ~gamma_b:params.gamma_b ~e_ex:params.e_ex
          ~d_f_coeff:params.d_f_coeff ~d_s_coeff:params.d_s_coeff
          ~model:(Params.model_to_int params.model)
          ~phase:(Params.phase_to_int params.phase)
          ~geometry:(Params.geometry_to_int params.geometry)
          ~depairing:(Params.depairing_to_tuple depairing)
          ~d_f_data ~tc_data
          ~fit_gamma:flags.fit_gamma ~fit_gamma_b:flags.fit_gamma_b
          ~fit_e_ex:flags.fit_e_ex ~fit_xi_f:flags.fit_xi_f in
      Ok Result.{ chi2 }
    with Failure msg -> Error msg
