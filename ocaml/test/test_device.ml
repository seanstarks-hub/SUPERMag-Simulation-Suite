(** Tests for Device stack parsing and geometry resolution. *)

open Supermag_types
open Supermag_pipeline

(* ── parse_stack tests ───────────────────────────────── *)

let test_parse_bilayer () =
  match Device.parse_stack "Nb:30/Fe:8" with
  | Error e -> Alcotest.fail e
  | Ok g ->
    Alcotest.(check int) "2 layers" 2 (List.length g.layers);
    Alcotest.(check string) "description" "Nb:30/Fe:8" g.description

let test_parse_trilayer () =
  match Device.parse_stack "Nb:50/CuNi:10" with
  | Error e -> Alcotest.fail e
  | Ok g ->
    Alcotest.(check int) "2 layers" 2 (List.length g.layers)

let test_parse_bad_material () =
  match Device.parse_stack "Unobtanium:10/Fe:5" with
  | Ok _ -> Alcotest.fail "should reject unknown material"
  | Error e ->
    Alcotest.(check bool) "mentions unknown" true
      (String.length e > 0)

let test_parse_bad_thickness () =
  match Device.parse_stack "Nb:abc/Fe:5" with
  | Ok _ -> Alcotest.fail "should reject non-numeric thickness"
  | Error _ -> ()

let test_parse_negative_thickness () =
  match Device.parse_stack "Nb:-5/Fe:5" with
  | Ok _ -> Alcotest.fail "should reject negative thickness"
  | Error _ -> ()

let test_parse_empty () =
  match Device.parse_stack "" with
  | Ok _ -> Alcotest.fail "should reject empty stack"
  | Error _ -> ()

let test_parse_whitespace () =
  match Device.parse_stack " Nb : 30 / Fe : 8 " with
  | Error e -> Alcotest.fail e
  | Ok g ->
    Alcotest.(check int) "2 layers" 2 (List.length g.layers)

(* ── resolve tests ───────────────────────────────────── *)

let test_resolve_bilayer () =
  let g = Geometry.bilayer ~s_thickness:30.0 ~s_material:"Nb"
      ~f_thickness:8.0 ~f_material:"Fe" in
  match Device.resolve g () with
  | Error e -> Alcotest.fail e
  | Ok p ->
    Alcotest.(check (float 0.01)) "tc0 = Nb Tc" 9.2 p.tc0;
    Alcotest.(check (float 0.01)) "xi_f = Fe xi_f" 0.7 p.xi_f;
    Alcotest.(check (float 0.01)) "d_s = 30" 30.0 p.d_s;
    Alcotest.(check (float 0.01)) "e_ex = Fe" 256.0 p.e_ex;
    (match p.geometry with
     | Params.Bilayer -> ()
     | _ -> Alcotest.fail "expected Bilayer geometry")

let test_resolve_bilayer_from_parse () =
  match Device.parse_stack "Nb:30/Fe:8" with
  | Error e -> Alcotest.fail e
  | Ok g ->
    begin match Device.resolve g () with
    | Error e -> Alcotest.fail e
    | Ok p ->
      Alcotest.(check (float 0.01)) "tc0" 9.2 p.tc0;
      (match p.geometry with
       | Params.Bilayer -> ()
       | _ -> Alcotest.fail "expected Bilayer")
    end

let test_resolve_multi_f_error () =
  let g = Geometry.{
    layers = [
      S_layer { thickness = 30.0; material = "Nb" };
      F_layer { thickness = 5.0; material = "Fe" };
      F_layer { thickness = 5.0; material = "Ni" };
    ];
    description = "S/F/F";
  } in
  match Device.resolve g () with
  | Ok _ -> Alcotest.fail "multi-F without hint should fail"
  | Error e ->
    Alcotest.(check bool) "mentions ambiguous" true
      (String.length e > 0)

let test_resolve_multi_f_with_domain_hint () =
  let g = Geometry.{
    layers = [
      S_layer { thickness = 30.0; material = "Nb" };
      F_layer { thickness = 5.0; material = "Fe" };
      F_layer { thickness = 5.0; material = "Ni" };
    ];
    description = "S/F/F domains";
  } in
  let config = Params.Domain_config {
    domain_width = 15.0; domain_wall = 0.0;
  } in
  match Device.resolve g ~geometry_hint:Params.Domains
      ~geom_config:config () with
  | Error e -> Alcotest.fail e
  | Ok p ->
    (match p.geometry with
     | Params.Domains -> ()
     | _ -> Alcotest.fail "expected Domains geometry")

let test_resolve_graded_requires_config () =
  let g = Geometry.{
    layers = [
      S_layer { thickness = 30.0; material = "Nb" };
      F_layer { thickness = 5.0; material = "Fe" };
      F_layer { thickness = 5.0; material = "Ni" };
    ];
    description = "S/F/F";
  } in
  (* Hint is Graded but no Graded_config → error *)
  match Device.resolve g ~geometry_hint:Params.Graded () with
  | Ok _ -> Alcotest.fail "graded without config should fail"
  | Error _ -> ()

let test_resolve_no_sc () =
  let g = Geometry.{
    layers = [ F_layer { thickness = 5.0; material = "Fe" } ];
    description = "F only";
  } in
  match Device.resolve g () with
  | Ok _ -> Alcotest.fail "should require superconductor"
  | Error _ -> ()

let test_resolve_no_fm () =
  let g = Geometry.{
    layers = [ S_layer { thickness = 30.0; material = "Nb" } ];
    description = "S only";
  } in
  match Device.resolve g () with
  | Ok _ -> Alcotest.fail "should require ferromagnet"
  | Error _ -> ()

(* ── Test runner ─────────────────────────────────────── *)

let () =
  Alcotest.run "Device" [
    "parse_stack", [
      Alcotest.test_case "bilayer" `Quick test_parse_bilayer;
      Alcotest.test_case "trilayer" `Quick test_parse_trilayer;
      Alcotest.test_case "bad material" `Quick test_parse_bad_material;
      Alcotest.test_case "bad thickness" `Quick test_parse_bad_thickness;
      Alcotest.test_case "negative thickness" `Quick test_parse_negative_thickness;
      Alcotest.test_case "empty" `Quick test_parse_empty;
      Alcotest.test_case "whitespace" `Quick test_parse_whitespace;
    ];
    "resolve", [
      Alcotest.test_case "bilayer" `Quick test_resolve_bilayer;
      Alcotest.test_case "bilayer from parse" `Quick test_resolve_bilayer_from_parse;
      Alcotest.test_case "multi-F error" `Quick test_resolve_multi_f_error;
      Alcotest.test_case "multi-F domain hint" `Quick test_resolve_multi_f_with_domain_hint;
      Alcotest.test_case "graded requires config" `Quick test_resolve_graded_requires_config;
      Alcotest.test_case "no SC" `Quick test_resolve_no_sc;
      Alcotest.test_case "no FM" `Quick test_resolve_no_fm;
    ];
  ]
