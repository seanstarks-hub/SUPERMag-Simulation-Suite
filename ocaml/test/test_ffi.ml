(** Tests for FFI stubs (ctypes bindings to C++ library). *)

let eps = 1e-10

let test_hbar () =
  let h = Supermag_ffi.Stubs.supermag_const_hbar () in
  Alcotest.(check bool) "hbar > 0" true (h > 0.0);
  Alcotest.(check bool) "hbar ~ 6.58e-16"
    true (Float.abs (h -. 6.582119514e-16) < 1e-25)

let test_kB () =
  let k = Supermag_ffi.Stubs.supermag_const_kB () in
  Alcotest.(check bool) "kB > 0" true (k > 0.0);
  Alcotest.(check bool) "kB ~ 8.617e-2"
    true (Float.abs (k -. 8.617333262e-2) < 1e-8)

let test_depairing_total () =
  let total = Supermag_ffi.Stubs.depairing_total
      ~ag:0.1 ~zeeman:0.2 ~orbital:0.05 ~spin_orbit:0.03 in
  Alcotest.(check bool) "sum ~ 0.38"
    true (Float.abs (total -. 0.38) < eps)

let test_pair_amplitude_zero () =
  let (x_arr, f_arr) = Supermag_ffi.Stubs.pair_amplitude
      ~d_f:5.0 ~xi_f:1.0 ~phase:0 ~n_points:100 in
  Alcotest.(check int) "x length" 100 (Array.length x_arr);
  Alcotest.(check int) "F length" 100 (Array.length f_arr);
  (* phase=0: F(0) = cos(0) = 1.0 *)
  Alcotest.(check bool) "F(0) ~ 1.0"
    true (Float.abs (f_arr.(0) -. 1.0) < 0.01)

let test_pair_amplitude_pi () =
  let (_x_arr, f_arr) = Supermag_ffi.Stubs.pair_amplitude
      ~d_f:5.0 ~xi_f:1.0 ~phase:1 ~n_points:100 in
  (* phase=pi: F(0) = sin(0) = 0.0 *)
  Alcotest.(check bool) "F(0) ~ 0.0"
    true (Float.abs f_arr.(0) < 0.01)

let test_bdg_eigenvalues () =
  let eigs = Supermag_ffi.Stubs.bdg_solve
      ~n_sites:10 ~t_hop:1.0 ~delta:1.5 ~e_ex:50.0 in
  Alcotest.(check int) "2*N eigenvalues" 20 (Array.length eigs);
  (* Check sorted *)
  let sorted = Array.to_list eigs |> List.sort Float.compare in
  Alcotest.(check bool) "sorted" true
    (Array.to_list eigs = sorted)

let test_josephson_output () =
  let (phi, current) = Supermag_ffi.Stubs.josephson_cpr
      ~d_f:1.0 ~xi_f:1.0 ~e_ex:100.0 ~t:4.0 ~n_phases:50 in
  Alcotest.(check int) "phi length" 50 (Array.length phi);
  Alcotest.(check int) "current length" 50 (Array.length current);
  (* Max current magnitude should be ~1 (normalized) *)
  let max_i = Array.fold_left (fun m c -> Float.max m (Float.abs c)) 0.0 current in
  Alcotest.(check bool) "max |I| ~ 1.0"
    true (Float.abs (max_i -. 1.0) < 0.3)

let test_triplet_parallel () =
  (* Parallel magnetization → no triplet component *)
  let (f_out, x_out) = Supermag_ffi.Stubs.triplet_solve
      ~n_layers:2
      ~thicknesses:[|5.0; 5.0|]
      ~magnetization_angles:[|0.0; 0.0|]
      ~n_grid:20 in
  Alcotest.(check int) "f length" 20 (Array.length f_out);
  Alcotest.(check int) "x length" 20 (Array.length x_out);
  let max_f = Array.fold_left (fun m v -> Float.max m (Float.abs v)) 0.0 f_out in
  Alcotest.(check bool) "parallel → near-zero triplet"
    true (max_f < 1e-6)

let test_usadel_output () =
  let (delta, x) = Supermag_ffi.Stubs.usadel_solve
      ~tc0:9.2 ~d_s:50.0 ~d_f:5.0 ~xi_s:38.0 ~xi_f:0.7
      ~e_ex:256.0 ~n_grid:20 in
  Alcotest.(check int) "delta length" 20 (Array.length delta);
  Alcotest.(check int) "x length" 20 (Array.length x);
  (* Delta >= 0 everywhere *)
  Array.iter (fun d ->
    Alcotest.(check bool) "delta >= 0" true (d >= 0.0)
  ) delta

let test_eilenberger_output () =
  let (f, x) = Supermag_ffi.Stubs.eilenberger_solve
      ~tc0:9.2 ~d_s:50.0 ~d_f:5.0 ~xi_s:38.0
      ~e_ex:256.0 ~n_grid:20 in
  Alcotest.(check int) "f length" 20 (Array.length f);
  Alcotest.(check int) "x length" 20 (Array.length x);
  (* f values should be bounded [0, 1] *)
  Array.iter (fun v ->
    Alcotest.(check bool) "f >= 0" true (v >= 0.0);
    Alcotest.(check bool) "f <= 1" true (v <= 1.0 +. eps)
  ) f

let test_gl_output () =
  let (psi_r, psi_i) = Supermag_ffi.Stubs.gl_minimize
      ~alpha:(-1.0) ~beta:1.0 ~kappa:1.0
      ~nx:10 ~ny:10 ~dx:0.5 in
  Alcotest.(check int) "psi_real length" 100 (Array.length psi_r);
  Alcotest.(check int) "psi_imag length" 100 (Array.length psi_i)

let () =
  Alcotest.run "FFI" [
    "constants", [
      Alcotest.test_case "hbar" `Quick test_hbar;
      Alcotest.test_case "kB"   `Quick test_kB;
    ];
    "depairing", [
      Alcotest.test_case "total" `Quick test_depairing_total;
    ];
    "proximity", [
      Alcotest.test_case "pair_amplitude zero" `Quick test_pair_amplitude_zero;
      Alcotest.test_case "pair_amplitude pi"   `Quick test_pair_amplitude_pi;
    ];
    "solvers", [
      Alcotest.test_case "bdg eigenvalues"     `Quick test_bdg_eigenvalues;
      Alcotest.test_case "josephson output"     `Quick test_josephson_output;
      Alcotest.test_case "triplet parallel"     `Quick test_triplet_parallel;
      Alcotest.test_case "usadel output"        `Quick test_usadel_output;
      Alcotest.test_case "eilenberger output"   `Quick test_eilenberger_output;
      Alcotest.test_case "gl output"            `Quick test_gl_output;
    ];
  ]
