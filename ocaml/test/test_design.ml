(** Tests for the Design combinatorial stack explorer. *)

open Supermag_types
open Supermag_pipeline

(* ── enumerate_bilayers tests ────────────────────────── *)

let test_enumerate_count () =
  let sc = [Material.nb; Material.pb; Material.al] in
  let fm = Material.all_ferromagnets in   (* 6 ferromagnets *)
  let stacks = Design.enumerate_bilayers ~sc ~fm ~d_s:30.0 ~d_f:5.0 in
  Alcotest.(check int) "3 SC × 6 FM = 18" 18 (List.length stacks)

let test_enumerate_materials () =
  let stacks = Design.enumerate_bilayers
    ~sc:Material.all_superconductors
    ~fm:Material.all_ferromagnets
    ~d_s:30.0 ~d_f:5.0 in
  (* Every geometry must contain exactly one S-layer and one F-layer *)
  List.iter (fun g ->
    let has_s = List.exists (function Geometry.S_layer _ -> true | _ -> false) g.Geometry.layers in
    let has_f = List.exists (function Geometry.F_layer _ -> true | _ -> false) g.Geometry.layers in
    Alcotest.(check bool) "has S-layer" true has_s;
    Alcotest.(check bool) "has F-layer" true has_f
  ) stacks

let test_enumerate_empty_sc () =
  let stacks = Design.enumerate_bilayers
    ~sc:[] ~fm:Material.all_ferromagnets ~d_s:30.0 ~d_f:5.0 in
  Alcotest.(check int) "empty SC → 0 stacks" 0 (List.length stacks)

let test_enumerate_empty_fm () =
  let stacks = Design.enumerate_bilayers
    ~sc:Material.all_superconductors ~fm:[] ~d_s:30.0 ~d_f:5.0 in
  Alcotest.(check int) "empty FM → 0 stacks" 0 (List.length stacks)

(* ── no_constraints tests ────────────────────────────── *)

let test_no_constraints () =
  let c = Design.no_constraints in
  Alcotest.(check bool) "tc_min is None" true (c.tc_min = None);
  Alcotest.(check bool) "tc_max is None" true (c.tc_max = None);
  Alcotest.(check bool) "max_total_thickness is None" true
    (c.max_total_thickness = None)

(* ── explore tests ───────────────────────────────────── *)

let test_explore_empty () =
  let result = Design.explore ~stacks:[] ~d_f_array:[|1.0; 5.0; 10.0|] () in
  Alcotest.(check int) "total_evaluated = 0" 0 result.Result.total_evaluated;
  Alcotest.(check int) "ranked is empty" 0 (List.length result.ranked)

let test_explore_description () =
  (* Each bilayer description must be non-empty *)
  let stacks = Design.enumerate_bilayers
    ~sc:[Material.nb] ~fm:[Material.fe] ~d_s:30.0 ~d_f:5.0 in
  List.iter (fun g ->
    Alcotest.(check bool) "non-empty description" true
      (String.length g.Geometry.description > 0)
  ) stacks

(* ── Test runner ─────────────────────────────────────── *)

let () =
  Alcotest.run "Design" [
    "enumerate_bilayers", [
      Alcotest.test_case "count 3×6=18" `Quick test_enumerate_count;
      Alcotest.test_case "materials valid" `Quick test_enumerate_materials;
      Alcotest.test_case "empty SC" `Quick test_enumerate_empty_sc;
      Alcotest.test_case "empty FM" `Quick test_enumerate_empty_fm;
    ];
    "no_constraints", [
      Alcotest.test_case "all None" `Quick test_no_constraints;
    ];
    "explore", [
      Alcotest.test_case "empty stacks" `Quick test_explore_empty;
      Alcotest.test_case "description non-empty" `Quick test_explore_description;
    ];
  ]
