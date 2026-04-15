(** Tests for solver chaining. *)

open Supermag_pipeline
open Supermag_types

let test_chain_identity () =
  let id_step = Chain.{ name = "identity"; run = Fun.id } in
  let result = Chain.chain [id_step; id_step] 42 in
  Alcotest.(check int) "identity chain" 42 result

let test_chain_increment () =
  let inc = Chain.{ name = "inc"; run = (fun x -> x + 1) } in
  let result = Chain.chain [inc; inc; inc] 0 in
  Alcotest.(check int) "three increments" 3 result

let test_chain_compose () =
  let double = Chain.{ name = "double"; run = (fun x -> x * 2) } in
  let add_one = Chain.{ name = "add_one"; run = (fun x -> x + 1) } in
  let result = Chain.chain [double; add_one; double] 3 in
  (* 3 -> 6 -> 7 -> 14 *)
  Alcotest.(check int) "compose" 14 result

let test_chain_empty () =
  let result = Chain.chain [] 99 in
  Alcotest.(check int) "empty chain" 99 result

let test_parallel_chain () =
  let double = Chain.{ name = "double"; run = (fun x -> x * 2) } in
  let triple = Chain.{ name = "triple"; run = (fun x -> x * 3) } in
  let results = Chain.parallel_chain [double; triple] [10; 10] in
  Alcotest.(check (list int)) "parallel" [20; 30] results

let test_parallel_chain_single () =
  let inc = Chain.{ name = "inc"; run = (fun x -> x + 1) } in
  let results = Chain.parallel_chain [inc] [5] in
  Alcotest.(check (list int)) "single parallel" [6] results

(* ── tc_to_minimum tests ─────────────────────────────── *)

let make_tc d_fs tcs tc0 : Result.tc_result =
  Result.{ d_f_values = Array.of_list d_fs;
           tc_values  = Array.of_list tcs;
           tc0 }

let test_tc_to_minimum_basic () =
  (* Minimum is at index 2, d_F = 3.0, Tc = 1.0 *)
  let tc = make_tc [1.0; 2.0; 3.0; 4.0; 5.0]
                   [8.0; 5.0; 1.0; 4.0; 6.0] 9.2 in
  match Chain.tc_to_minimum tc with
  | Error e -> Alcotest.fail e
  | Ok opt ->
    Alcotest.(check (float 1e-6)) "d_f_optimal" 3.0 opt.d_f_optimal

let test_tc_to_minimum_first () =
  (* Minimum at first element *)
  let tc = make_tc [1.0; 2.0; 3.0]
                   [0.5; 5.0; 8.0] 9.2 in
  match Chain.tc_to_minimum tc with
  | Error e -> Alcotest.fail e
  | Ok opt ->
    Alcotest.(check (float 1e-6)) "d_f_optimal first" 1.0 opt.d_f_optimal

let test_tc_to_minimum_last () =
  (* Minimum at last element *)
  let tc = make_tc [1.0; 2.0; 3.0]
                   [8.0; 5.0; 0.1] 9.2 in
  match Chain.tc_to_minimum tc with
  | Error e -> Alcotest.fail e
  | Ok opt ->
    Alcotest.(check (float 1e-6)) "d_f_optimal last" 3.0 opt.d_f_optimal

let test_tc_to_minimum_empty () =
  let tc = make_tc [] [] 9.2 in
  match Chain.tc_to_minimum tc with
  | Ok _ -> Alcotest.fail "empty tc_result should return Error"
  | Error _ -> ()

(* ── tc_then_usadel error case ───────────────────────── *)

let test_tc_then_usadel_empty_array () =
  (* Empty d_f_array should return Error without calling the solver *)
  let params = Params.{
    tc0 = 9.2; d_s = 30.0; d_s_coeff = 18.0; xi_s = 38.0; xi_f = 0.7;
    gamma = 1.0; gamma_b = 0.0; e_ex = 256.0; d_f_coeff = 2.5e-4;
    model = Fominov; phase = Zero; geometry = Bilayer;
    geom_config = None; spin_active = None;
  } in
  match Chain.tc_then_usadel ~params ~d_f_array:[||] ~t:4.2 () with
  | Ok _ -> Alcotest.fail "empty d_f_array should fail"
  | Error _ -> ()

(* ── usadel_then_josephson error case ───────────────── *)

let test_usadel_then_josephson_bad_d_f () =
  match Chain.usadel_then_josephson
    ~tc0:9.2 ~d_s:30.0 ~d_f:(-1.0) ~xi_s:38.0 ~xi_f:0.7
    ~e_ex:256.0 ~t:4.2 ~gamma_b:0.0 () with
  | Ok _ -> Alcotest.fail "negative d_f should fail"
  | Error _ -> ()

let () =
  Alcotest.run "Chain" [
    "sequential", [
      Alcotest.test_case "identity"  `Quick test_chain_identity;
      Alcotest.test_case "increment" `Quick test_chain_increment;
      Alcotest.test_case "compose"   `Quick test_chain_compose;
      Alcotest.test_case "empty"     `Quick test_chain_empty;
    ];
    "parallel", [
      Alcotest.test_case "basic"  `Quick test_parallel_chain;
      Alcotest.test_case "single" `Quick test_parallel_chain_single;
    ];
    "tc_to_minimum", [
      Alcotest.test_case "basic"  `Quick test_tc_to_minimum_basic;
      Alcotest.test_case "first"  `Quick test_tc_to_minimum_first;
      Alcotest.test_case "last"   `Quick test_tc_to_minimum_last;
      Alcotest.test_case "empty"  `Quick test_tc_to_minimum_empty;
    ];
    "tc_then_usadel", [
      Alcotest.test_case "empty array" `Quick test_tc_then_usadel_empty_array;
    ];
    "usadel_then_josephson", [
      Alcotest.test_case "bad d_f" `Quick test_usadel_then_josephson_bad_d_f;
    ];
  ]
