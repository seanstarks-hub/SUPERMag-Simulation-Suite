(** Device stack parsing and geometry resolution.

    Parses "Nb:30/Fe:8" notation into [Geometry.geometry], then
    resolves layer stacks to [Params.proximity_params] by looking up
    material properties and auto-detecting bilayer/trilayer geometry.
    Graded and domain geometries require explicit annotation. *)

open Supermag_types

(* ── Stack notation parser ───────────────────────────── *)

let parse_layer (s : string) : (Geometry.layer, string) result =
  match String.split_on_char ':' (String.trim s) with
  | [name; thick_s] ->
    (match float_of_string_opt (String.trim thick_s) with
     | None -> Error (Printf.sprintf "invalid thickness: %s" thick_s)
     | Some d when d <= 0.0 ->
       Error (Printf.sprintf "thickness must be > 0: %s" name)
     | Some d ->
       let mat = String.trim name in
       if Material.get_superconductor mat <> None then
         Ok (Geometry.S_layer { thickness = d; material = mat })
       else if Material.get_ferromagnet mat <> None then
         Ok (Geometry.F_layer { thickness = d; material = mat })
       else
         Error (Printf.sprintf "unknown material: %s" mat))
  | _ -> Error (Printf.sprintf "expected Material:thickness, got: %s" s)

let parse_stack (s : string) : (Geometry.geometry, string) result =
  let parts = String.split_on_char '/' s in
  let rec go acc = function
    | [] ->
      if acc = [] then Error "empty stack"
      else
        let layers = List.rev acc in
        let desc = String.concat "/" (List.map (fun p -> String.trim p) parts) in
        Ok Geometry.{ layers; description = desc }
    | p :: rest ->
      begin match parse_layer p with
      | Error e -> Error e
      | Ok layer -> go (layer :: acc) rest
      end
  in
  go [] parts

(* ── Geometry resolution ─────────────────────────────── *)

(** Classify layers into S, F, N, I counts for auto-detection. *)
let classify_layers (layers : Geometry.layer list) =
  List.fold_left (fun (s, f, n, i) layer ->
    match layer with
    | Geometry.S_layer _ -> (s + 1, f, n, i)
    | Geometry.F_layer _ -> (s, f + 1, n, i)
    | Geometry.N_layer _ -> (s, f, n + 1, i)
    | Geometry.I_layer _ -> (s, f, n, i + 1)
  ) (0, 0, 0, 0) layers

(** Extract the first S-layer material, or None. *)
let first_sc (layers : Geometry.layer list) =
  List.find_map (fun l ->
    match l with
    | Geometry.S_layer { material; _ } -> Material.get_superconductor material
    | _ -> None
  ) layers

(** Extract the first F-layer material, or None. *)
let first_fm (layers : Geometry.layer list) =
  List.find_map (fun l ->
    match l with
    | Geometry.F_layer { material; _ } -> Material.get_ferromagnet material
    | _ -> None
  ) layers

(** Extract S-layer thickness, or 0. *)
let s_thickness (layers : Geometry.layer list) =
  List.fold_left (fun acc l ->
    match l with
    | Geometry.S_layer { thickness; _ } -> acc +. thickness
    | _ -> acc
  ) 0.0 layers

(** Extract first N-layer for trilayer config. *)
let first_n_layer (layers : Geometry.layer list) =
  List.find_map (fun l ->
    match l with
    | Geometry.N_layer { thickness } -> Some thickness
    | _ -> None
  ) layers

let resolve (geom : Geometry.geometry)
    ?geometry_hint ?geom_config () :
    (Params.proximity_params, string) result =
  let layers = geom.layers in
  let (n_s, n_f, n_n, _n_i) = classify_layers layers in
  if n_s < 1 then Error "stack must contain at least one superconductor"
  else if n_f < 1 then Error "stack must contain at least one ferromagnet"
  else
    match first_sc layers, first_fm layers with
    | None, _ -> Error "superconductor material not found in database"
    | _, None -> Error "ferromagnet material not found in database"
    | Some sc, Some fm ->
      let d_s = s_thickness layers in
      (* Auto-detect geometry *)
      let detect () =
        if n_s = 1 && n_f = 1 && n_n = 0 then
          Ok (Params.Bilayer, None)
        else if n_s = 1 && n_f = 1 && n_n = 1 then
          begin match first_n_layer layers with
          | Some d_n ->
            Ok (Params.Trilayer,
                Some (Params.Trilayer_config {
                  d_n;
                  xi_n = 100.0;  (* default normal-metal coherence length *)
                  r_b = 0.0;     (* default: transparent interface *)
                }))
          | None -> Error "trilayer detected but N-layer thickness missing"
          end
        else
          (* Multi-F or complex: need explicit geometry *)
          begin match geometry_hint with
          | Some hint -> Ok (hint, geom_config)
          | None ->
            Error "ambiguous stack: provide ~geometry_hint for multi-F stacks"
          end
      in
      begin match detect () with
      | Error e -> Error e
      | Ok (geometry, resolved_config) ->
        (* Validate config for graded/domains *)
        let config_result = match geometry with
          | Params.Graded ->
            (match resolved_config with
             | Some (Params.Graded_config _) -> Ok resolved_config
             | _ -> Error "graded geometry requires Graded_config")
          | Params.Domains ->
            (match resolved_config with
             | Some (Params.Domain_config _) -> Ok resolved_config
             | _ -> Error "domain geometry requires Domain_config")
          | _ -> Ok resolved_config
        in
        match config_result with
        | Error e -> Error e
        | Ok validated_config ->
          Ok Params.{
            tc0 = sc.tc;
            d_s;
            d_s_coeff = sc.d_s;
            xi_s = sc.xi_s;
            xi_f = fm.xi_f;
            gamma = 1.0;       (* default coupling *)
            gamma_b = 0.0;     (* default: transparent *)
            e_ex = fm.e_ex;
            d_f_coeff = fm.d_f;
            model = Fominov;   (* default model *)
            phase = Zero;      (* default phase *)
            geometry;
            geom_config = validated_config;
            spin_active = None;
          }
      end
