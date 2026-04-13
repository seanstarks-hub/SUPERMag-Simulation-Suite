(** Tests for parameter sweep logic. *)

open Supermag_pipeline

let eps = 1e-10

let test_grid_sweep_count () =
  let config = Sweep.{
    param_name = "d_F"; min_val = 0.5; max_val = 20.0;
    n_points = 50; sweep_type = Grid;
  } in
  let arr = Sweep.grid_sweep config in
  Alcotest.(check int) "n_points" 50 (Array.length arr)

let test_grid_sweep_range () =
  let config = Sweep.{
    param_name = "d_F"; min_val = 1.0; max_val = 10.0;
    n_points = 11; sweep_type = Grid;
  } in
  let arr = Sweep.grid_sweep config in
  Alcotest.(check bool) "first = min" true
    (Float.abs (arr.(0) -. 1.0) < eps);
  Alcotest.(check bool) "last = max" true
    (Float.abs (arr.(10) -. 10.0) < eps)

let test_grid_sweep_monotonic () =
  let config = Sweep.{
    param_name = "x"; min_val = 0.0; max_val = 100.0;
    n_points = 100; sweep_type = Grid;
  } in
  let arr = Sweep.grid_sweep config in
  let sorted = ref true in
  for i = 1 to Array.length arr - 1 do
    if arr.(i) < arr.(i-1) then sorted := false
  done;
  Alcotest.(check bool) "monotonically increasing" true !sorted

let test_random_sweep_count () =
  let config = Sweep.{
    param_name = "d_F"; min_val = 0.5; max_val = 20.0;
    n_points = 30; sweep_type = Random;
  } in
  let arr = Sweep.random_sweep config in
  Alcotest.(check int) "n_points" 30 (Array.length arr)

let test_random_sweep_range () =
  let config = Sweep.{
    param_name = "d_F"; min_val = 1.0; max_val = 5.0;
    n_points = 100; sweep_type = Random;
  } in
  let arr = Sweep.random_sweep config in
  Array.iter (fun v ->
    Alcotest.(check bool) "v >= min" true (v >= 1.0);
    Alcotest.(check bool) "v <= max" true (v <= 5.0);
  ) arr

let test_adaptive_sweep_count () =
  let config = Sweep.{
    param_name = "x"; min_val = 0.0; max_val = 10.0;
    n_points = 30; sweep_type = Adaptive;
  } in
  let arr = Sweep.adaptive_sweep config (fun x -> Float.sin x) in
  (* Adaptive may produce more or fewer points than requested,
     but should be in the right ballpark *)
  Alcotest.(check bool) "reasonable count" true
    (Array.length arr >= 5 && Array.length arr <= 100)

let test_grid_single_point () =
  let config = Sweep.{
    param_name = "x"; min_val = 3.0; max_val = 7.0;
    n_points = 1; sweep_type = Grid;
  } in
  let arr = Sweep.grid_sweep config in
  Alcotest.(check int) "single point" 1 (Array.length arr);
  Alcotest.(check bool) "value = min" true
    (Float.abs (arr.(0) -. 3.0) < eps)

let test_sweep_param_parsing () =
  Alcotest.(check bool) "d_F ok" true
    (Result.is_ok (Sweep.sweep_param_of_string "d_F"));
  Alcotest.(check bool) "gamma ok" true
    (Result.is_ok (Sweep.sweep_param_of_string "gamma"));
  Alcotest.(check bool) "D_S ok" true
    (Result.is_ok (Sweep.sweep_param_of_string "D_S"));
  Alcotest.(check bool) "bad err" true
    (Result.is_error (Sweep.sweep_param_of_string "banana"))

let test_sweep_param_roundtrip () =
  let params = ["d_F"; "d_S"; "D_S"; "gamma"; "gamma_B"; "E_ex"; "xi_F"; "Tc0"] in
  List.iter (fun name ->
    match Sweep.sweep_param_of_string name with
    | Error _ -> Alcotest.fail (Printf.sprintf "parse %s" name)
    | Ok p ->
      Alcotest.(check string) (Printf.sprintf "%s roundtrip" name)
        name (Sweep.sweep_param_to_string p)
  ) params

open Supermag_types

let test_geometry_parsing () =
  let cases = ["bilayer"; "trilayer"; "graded"; "domains"] in
  List.iter (fun s ->
    Alcotest.(check bool) (Printf.sprintf "geometry %s" s) true
      (Result.is_ok (Params.geometry_of_string s))
  ) cases;
  Alcotest.(check bool) "geometry bad" true
    (Result.is_error (Params.geometry_of_string "unknown"))

let test_model_parsing () =
  let cases = ["thin_s"; "fominov"; "fominov_multi"] in
  List.iter (fun s ->
    Alcotest.(check bool) (Printf.sprintf "model %s" s) true
      (Result.is_ok (Params.model_of_string s))
  ) cases;
  Alcotest.(check bool) "model bad" true
    (Result.is_error (Params.model_of_string "bogus"))

let test_grade_profile_to_int () =
  Alcotest.(check int) "Linear" 0 (Params.grade_profile_to_int Params.Linear);
  Alcotest.(check int) "Exponential" 1 (Params.grade_profile_to_int Params.Exponential);
  Alcotest.(check int) "Step" 2 (Params.grade_profile_to_int Params.Step)

let test_geometry_to_int () =
  Alcotest.(check int) "Bilayer" 0 (Params.geometry_to_int Params.Bilayer);
  Alcotest.(check int) "Trilayer" 1 (Params.geometry_to_int Params.Trilayer);
  Alcotest.(check int) "Graded" 2 (Params.geometry_to_int Params.Graded);
  Alcotest.(check int) "Domains" 3 (Params.geometry_to_int Params.Domains)

let test_material_rho () =
  let sc = Material.nb in
  Alcotest.(check bool) "Nb rho > 0" true (sc.rho > 0.0);
  Alcotest.(check bool) "Nb d_s > 0" true (sc.d_s > 0.0);
  let fm = Material.fe in
  Alcotest.(check bool) "Fe rho > 0" true (fm.rho > 0.0)

let () =
  Alcotest.run "Sweep" [
    "grid", [
      Alcotest.test_case "count"     `Quick test_grid_sweep_count;
      Alcotest.test_case "range"     `Quick test_grid_sweep_range;
      Alcotest.test_case "monotonic" `Quick test_grid_sweep_monotonic;
      Alcotest.test_case "single"    `Quick test_grid_single_point;
    ];
    "random", [
      Alcotest.test_case "count" `Quick test_random_sweep_count;
      Alcotest.test_case "range" `Quick test_random_sweep_range;
    ];
    "adaptive", [
      Alcotest.test_case "count" `Quick test_adaptive_sweep_count;
    ];
    "param", [
      Alcotest.test_case "parsing"    `Quick test_sweep_param_parsing;
      Alcotest.test_case "roundtrip"  `Quick test_sweep_param_roundtrip;
    ];
    "types", [
      Alcotest.test_case "geometry_parsing"   `Quick test_geometry_parsing;
      Alcotest.test_case "model_parsing"      `Quick test_model_parsing;
      Alcotest.test_case "grade_profile_int"  `Quick test_grade_profile_to_int;
      Alcotest.test_case "geometry_int"       `Quick test_geometry_to_int;
      Alcotest.test_case "material_rho"       `Quick test_material_rho;
    ];
  ]
