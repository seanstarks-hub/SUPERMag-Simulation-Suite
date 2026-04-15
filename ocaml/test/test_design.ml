(** Tests for combinatorial design exploration. *)

open Supermag_types
open Supermag_pipeline

(* ── enumerate tests ─────────────────────────────────── *)

let test_enumerate_count () =
  let results = Design.enumerate_bilayers
      ~sc:Material.all_superconductors
      ~fm:Material.all_ferromagnets
      ~d_s:50.0 ~d_f_array:[|1.0; 5.0; 10.0|] () in
  (* 3 SCs × 6 FMs = 18 combinations *)
  Alcotest.(check int) "18 combos" 18 (List.length results)

let test_enumerate_subset () =
  let results = Design.enumerate_bilayers
      ~sc:[Material.nb] ~fm:[Material.fe]
      ~d_s:50.0 ~d_f_array:[|1.0; 5.0|] () in
  Alcotest.(check int) "1 combo" 1 (List.length results);
  let r = List.hd results in
  Alcotest.(check string) "sc = Nb" "Nb" r.sc_name;
  Alcotest.(check string) "fm = Fe" "Fe" r.fm_name;
  Alcotest.(check (float 0.01)) "tc0 = 9.2" 9.2 r.tc_result.tc0;
  Alcotest.(check int) "2 d_f points" 2 (Array.length r.tc_result.d_f_values)

(* ── filter tests ────────────────────────────────────── *)

let test_filter_tc_min () =
  let results = Design.enumerate_bilayers
      ~sc:Material.all_superconductors
      ~fm:[Material.fe]
      ~d_s:50.0 ~d_f_array:[|1.0|] () in
  (* Al has Tc=1.2 K, so filtering at tc_min=5.0 should exclude it *)
  let filtered = Design.filter
      { Design.tc_min = Some 5.0; tc_max = None; max_d_total = None }
      results in
  (* All remaining should have at least one Tc >= 5 *)
  List.iter (fun r ->
    let max_tc = Array.fold_left Float.max 0.0 r.Design.tc_result.tc_values in
    Alcotest.(check bool) "Tc >= 5" true (max_tc >= 5.0)
  ) filtered

let test_filter_empty () =
  let results = Design.enumerate_bilayers
      ~sc:[Material.nb] ~fm:[Material.fe]
      ~d_s:50.0 ~d_f_array:[|1.0|] () in
  (* Impossible constraint: max total thickness = 1 nm *)
  let filtered = Design.filter
      { Design.tc_min = None; tc_max = None; max_d_total = Some 1.0 }
      results in
  Alcotest.(check int) "all filtered out" 0 (List.length filtered)

let test_filter_none () =
  let results = Design.enumerate_bilayers
      ~sc:[Material.nb; Material.pb] ~fm:[Material.fe]
      ~d_s:50.0 ~d_f_array:[|1.0|] () in
  let filtered = Design.filter Design.no_constraints results in
  Alcotest.(check int) "unchanged" (List.length results) (List.length filtered)

(* ── Test runner ─────────────────────────────────────── *)

let () =
  Alcotest.run "Design" [
    "enumerate", [
      Alcotest.test_case "full count" `Quick test_enumerate_count;
      Alcotest.test_case "subset" `Quick test_enumerate_subset;
    ];
    "filter", [
      Alcotest.test_case "tc_min" `Quick test_filter_tc_min;
      Alcotest.test_case "impossible" `Quick test_filter_empty;
      Alcotest.test_case "no constraints" `Quick test_filter_none;
    ];
  ]
