(** Tests for the multi-parameter design optimizer. *)

open Supermag_types
open Supermag_pipeline

(* ── Nelder-Mead on a known quadratic ────────────────── *)

(** Minimize f(x,y) = (x-3)^2 + (y-7)^2.
    Global minimum at (3, 7), f = 0. *)
let test_nelder_mead_quadratic () =
  let specs = [|
    ((fun _ -> 0.0), (fun p _ -> p), 0.0, 10.0);
    ((fun _ -> 0.0), (fun p _ -> p), 0.0, 15.0);
  |] in
  let f vec =
    let x = vec.(0) and y = vec.(1) in
    Ok ((x -. 3.0) ** 2.0 +. (y -. 7.0) ** 2.0)
  in
  let x0 = [| 1.0; 1.0 |] in
  match Optimize.nelder_mead ~f ~x0 ~specs
          ~free:Optimize.no_free ~constraints:Optimize.no_constraints
          ~tol:1e-8 ~max_eval:1000 with
  | Error msg -> Alcotest.fail msg
  | Ok (best, best_val, _evals) ->
    Alcotest.(check (float 0.1)) "x near 3" 3.0 best.(0);
    Alcotest.(check (float 0.1)) "y near 7" 7.0 best.(1);
    Alcotest.(check bool) "f near 0" true (best_val < 0.01)

(* ── Constraint projection ───────────────────────────── *)

let test_project_clips () =
  let specs = [|
    ((fun _ -> 0.0), (fun p _ -> p), 1.0, 10.0);
    ((fun _ -> 0.0), (fun p _ -> p), 2.0, 8.0);
  |] in
  let vec = [| -5.0; 20.0 |] in
  let out = Optimize.project specs Optimize.no_free
      Optimize.no_constraints vec in
  Alcotest.(check (float 1e-10)) "clipped lo" 1.0 out.(0);
  Alcotest.(check (float 1e-10)) "clipped hi" 8.0 out.(1)

let test_project_passthrough () =
  let specs = [|
    ((fun _ -> 0.0), (fun p _ -> p), 0.0, 100.0);
    ((fun _ -> 0.0), (fun p _ -> p), 0.0, 100.0);
  |] in
  let vec = [| 5.0; 7.0 |] in
  let out = Optimize.project specs Optimize.no_free
      Optimize.no_constraints vec in
  Alcotest.(check (float 1e-10)) "unchanged x" 5.0 out.(0);
  Alcotest.(check (float 1e-10)) "unchanged y" 7.0 out.(1)

(* ── sensitivity_at ──────────────────────────────────── *)

let test_sensitivity_returns () =
  let params = Params.{
    tc0 = 9.2; d_s = 50.0; d_s_coeff = 1.0;
    xi_s = 38.0; xi_f = 0.7; gamma = 0.3; gamma_b = 0.0;
    e_ex = 256.0; d_f_coeff = 2.5e-4;
    model = Fominov; phase = Zero; geometry = Bilayer;
    geom_config = None; spin_active = None;
  } in
  match Optimize.sensitivity_at params ~d_f:5.0 () with
  | Error msg -> Alcotest.fail msg
  | Ok s ->
    (* d_tc_d_df should be negative (Tc decreases with d_F near minimum) or
       small near the minimum; just check it's finite *)
    Alcotest.(check bool) "dTc/ddf finite"
      true (Float.is_finite s.d_tc_d_df);
    Alcotest.(check bool) "dTc/dds finite"
      true (Float.is_finite s.d_tc_d_ds);
    Alcotest.(check bool) "dTc/dgamma finite"
      true (Float.is_finite s.d_tc_d_gamma)

(* ── Target_tc with single free d_F ──────────────────── *)

let test_optimize_target_tc_1d () =
  let params = Params.{
    tc0 = 9.2; d_s = 50.0; d_s_coeff = 1.0;
    xi_s = 38.0; xi_f = 0.7; gamma = 0.3; gamma_b = 0.0;
    e_ex = 256.0; d_f_coeff = 2.5e-4;
    model = Fominov; phase = Zero; geometry = Bilayer;
    geom_config = None; spin_active = None;
  } in
  let problem = Optimize.{
    base_params = params;
    free = { no_free with vary_d_f = Some (0.5, 20.0) };
    objective = Target_tc 5.0;
    constraints = no_constraints;
    depairing = Params.no_depairing;
    tolerance = 1e-4;
    max_evaluations = 200;
  } in
  match Optimize.optimize problem with
  | Error msg -> Alcotest.fail msg
  | Ok r ->
    (* Should find some d_F that gives Tc near 5.0 K *)
    Alcotest.(check bool) "Tc near target"
      true (Float.abs (r.tc_achieved -. 5.0) < 1.0);
    Alcotest.(check bool) "d_F in range"
      true (r.d_f_optimal >= 0.5 && r.d_f_optimal <= 20.0)

(* ── Minimize_tc with no free params ─────────────────── *)

let test_optimize_no_free () =
  let params = Params.{
    tc0 = 9.2; d_s = 50.0; d_s_coeff = 1.0;
    xi_s = 38.0; xi_f = 0.7; gamma = 0.3; gamma_b = 0.0;
    e_ex = 256.0; d_f_coeff = 2.5e-4;
    model = Fominov; phase = Zero; geometry = Bilayer;
    geom_config = None; spin_active = None;
  } in
  let problem = Optimize.{
    base_params = params;
    free = no_free;
    objective = Minimize_tc;
    constraints = no_constraints;
    depairing = Params.no_depairing;
    tolerance = 1e-4;
    max_evaluations = 10;
  } in
  match Optimize.optimize problem with
  | Error msg -> Alcotest.fail msg
  | Ok r ->
    Alcotest.(check bool) "Tc finite" true (Float.is_finite r.tc_achieved);
    Alcotest.(check int) "1 evaluation" 1 r.evaluations

(* ── free_param_specs count ──────────────────────────── *)

let test_free_param_count () =
  let free = Optimize.{
    vary_d_s = Some (10.0, 100.0);
    vary_d_f = Some (0.5, 20.0);
    vary_gamma = None;
    vary_gamma_b = None;
    vary_e_ex = None;
  } in
  let specs = Optimize.free_param_specs free in
  Alcotest.(check int) "2 free params" 2 (Array.length specs)

(* ── Test runner ─────────────────────────────────────── *)

let () =
  Alcotest.run "Optimize" [
    "nelder_mead", [
      Alcotest.test_case "quadratic" `Quick test_nelder_mead_quadratic;
    ];
    "projection", [
      Alcotest.test_case "clips bounds" `Quick test_project_clips;
      Alcotest.test_case "passthrough" `Quick test_project_passthrough;
    ];
    "sensitivity", [
      Alcotest.test_case "returns finite" `Quick test_sensitivity_returns;
    ];
    "optimize", [
      Alcotest.test_case "target_tc 1D" `Quick test_optimize_target_tc_1d;
      Alcotest.test_case "no free params" `Quick test_optimize_no_free;
    ];
    "free_params", [
      Alcotest.test_case "count" `Quick test_free_param_count;
    ];
  ]
