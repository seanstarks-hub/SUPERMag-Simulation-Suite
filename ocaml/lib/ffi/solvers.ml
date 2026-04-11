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
    let tc_arr = Stubs.solve_tc_batch
        ~tc0:params.tc0 ~d_s:params.d_s ~xi_s:params.xi_s ~xi_f:params.xi_f
        ~gamma:params.gamma ~gamma_b:params.gamma_b ~e_ex:params.e_ex
        ~d_f_coeff:params.d_f_coeff
        ~model:(Params.model_to_int params.model)
        ~phase:(Params.phase_to_int params.phase)
        ~depairing:(Params.depairing_to_tuple depairing)
        ~d_f_arr:d_f_array in
    Ok Result.{ d_f_values = d_f_array; tc_values = tc_arr; tc0 = params.tc0 }

let pair_amplitude ~d_f ~xi_f ~phase ~n_points =
  if d_f <= 0.0 then Error "d_f must be > 0"
  else if xi_f <= 0.0 then Error "xi_f must be > 0"
  else if n_points < 2 then Error "n_points must be >= 2"
  else
    let (x, f) = Stubs.pair_amplitude
        ~d_f ~xi_f ~phase:(Params.phase_to_int phase) ~n_points in
    Ok Result.{ x; values = f }

let usadel ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~n_grid =
  if n_grid < 5 then Error "n_grid must be >= 5"
  else
    let (delta, x) = Stubs.usadel_solve ~tc0 ~d_s ~d_f ~xi_s ~xi_f ~e_ex ~n_grid in
    Ok Result.{ x; values = delta }

let eilenberger ~tc0 ~d_s ~d_f ~xi_s ~e_ex ~n_grid =
  if n_grid < 5 then Error "n_grid must be >= 5"
  else
    let (f, x) = Stubs.eilenberger_solve ~tc0 ~d_s ~d_f ~xi_s ~e_ex ~n_grid in
    Ok Result.{ x; values = f }

let bdg ~n_sites ~t_hop ~delta ~e_ex =
  if n_sites < 2 then Error "n_sites must be >= 2"
  else
    let eigs = Stubs.bdg_solve ~n_sites ~t_hop ~delta ~e_ex in
    Ok Result.{ eigenvalues = eigs }

let gl ~alpha ~beta ~kappa ~nx ~ny ~dx =
  if nx < 2 || ny < 2 then Error "nx and ny must be >= 2"
  else if dx <= 0.0 then Error "dx must be > 0"
  else
    let (psi_r, psi_i) = Stubs.gl_minimize ~alpha ~beta ~kappa ~nx ~ny ~dx in
    Ok Result.{ nx; ny; psi_real = psi_r; psi_imag = psi_i }

let josephson ~d_f ~xi_f ~e_ex ~t ~n_phases =
  if n_phases < 2 then Error "n_phases must be >= 2"
  else if d_f <= 0.0 then Error "d_f must be > 0"
  else
    let (phi, current) = Stubs.josephson_cpr ~d_f ~xi_f ~e_ex ~t ~n_phases in
    Ok Result.{ phi; current }

let triplet ~n_layers ~thicknesses ~magnetization_angles ~n_grid =
  if n_layers < 2 then Error "n_layers must be >= 2"
  else if Array.length thicknesses <> n_layers then
    Error "thicknesses length must equal n_layers"
  else if Array.length magnetization_angles <> n_layers then
    Error "magnetization_angles length must equal n_layers"
  else if n_grid < 5 then Error "n_grid must be >= 5"
  else
    let (f, x) = Stubs.triplet_solve
        ~n_layers ~thicknesses ~magnetization_angles ~n_grid in
    Ok Result.{ x; values = f }
