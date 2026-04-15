(** Combinatorial bilayer design exploration.

    Enumerates SC×FM material combinations, runs Tc(d_F) sweeps
    in parallel via Domain.spawn, and filters results by constraints. *)

open Supermag_types

(* ── Types ───────────────────────────────────────────── *)

type constraint_t = {
  tc_min : float option;
  tc_max : float option;
  max_d_total : float option;
}

let no_constraints = { tc_min = None; tc_max = None; max_d_total = None }

type design_result = {
  sc_name : string;
  fm_name : string;
  d_s : float;
  tc_result : Result.tc_result;
}

(* ── Enumeration ─────────────────────────────────────── *)

let enumerate_bilayers ~(sc : Material.superconductor list)
    ~(fm : Material.ferromagnet list) ~d_s ~d_f_array
    ?(depairing = Params.no_depairing) () : design_result list =
  (* Build all SC×FM combinations *)
  let combos = List.concat_map (fun s ->
    List.map (fun f -> (s, f)) fm
  ) sc in
  (* Spawn one domain per combination *)
  let domains = List.map (fun (s, f) ->
    Domain.spawn (fun () ->
      let geom = Geometry.bilayer
          ~s_thickness:d_s ~s_material:s.Material.name
          ~f_thickness:1.0 ~f_material:f.Material.name in
      match Device.resolve geom () with
      | Error msg -> Error msg
      | Ok params ->
        let params = Params.{ params with d_s } in
        match Supermag_ffi.Solvers.solve_tc ~params ~d_f_array ~depairing () with
        | Error msg -> Error msg
        | Ok tc_result ->
          Ok { sc_name = s.Material.name;
               fm_name = f.Material.name;
               d_s;
               tc_result }
    )
  ) combos in
  (* Join and collect successes *)
  List.filter_map (fun d ->
    match Domain.join d with
    | Ok r -> Some r
    | Error msg ->
      Printf.eprintf "Warning: %s\n" msg;
      None
  ) domains

(* ── Filtering ───────────────────────────────────────── *)

let filter (c : constraint_t) (results : design_result list) : design_result list =
  List.filter (fun r ->
    let tc_vals = r.tc_result.tc_values in
    let n = Array.length tc_vals in
    let pass_tc_min = match c.tc_min with
      | None -> true
      | Some tmin ->
        (* At least one Tc value must meet the minimum *)
        let found = ref false in
        for i = 0 to n - 1 do
          if tc_vals.(i) >= tmin then found := true
        done;
        !found
    in
    let pass_tc_max = match c.tc_max with
      | None -> true
      | Some tmax ->
        (* All Tc values must be at or below the maximum *)
        let ok = ref true in
        for i = 0 to n - 1 do
          if tc_vals.(i) > tmax then ok := false
        done;
        !ok
    in
    let pass_d_total = match c.max_d_total with
      | None -> true
      | Some dmax ->
        let d_f_max = Array.fold_left Float.max 0.0
            r.tc_result.d_f_values in
        r.d_s +. d_f_max <= dmax
    in
    pass_tc_min && pass_tc_max && pass_d_total
  ) results

(* ── Output ──────────────────────────────────────────── *)

let to_csv (oc : out_channel) (results : design_result list) : unit =
  Printf.fprintf oc "# SUPERMag design exploration\n";
  Printf.fprintf oc "SC,FM,d_S_nm,d_F_nm,Tc_K\n";
  List.iter (fun r ->
    let n = Array.length r.tc_result.d_f_values in
    for i = 0 to n - 1 do
      Printf.fprintf oc "%s,%s,%.6f,%.6f,%.6f\n"
        r.sc_name r.fm_name r.d_s
        r.tc_result.d_f_values.(i) r.tc_result.tc_values.(i)
    done
  ) results

let to_json (oc : out_channel) (results : design_result list) : unit =
  Printf.fprintf oc "{\n  \"design_results\": [\n";
  let nr = List.length results in
  List.iteri (fun ri r ->
    let n = Array.length r.tc_result.d_f_values in
    for i = 0 to n - 1 do
      let last = ri = nr - 1 && i = n - 1 in
      let comma = if last then "" else "," in
      Printf.fprintf oc
        "    {\"SC\": \"%s\", \"FM\": \"%s\", \"d_S\": %.6f, \"d_F\": %.6f, \"Tc_K\": %.6f}%s\n"
        r.sc_name r.fm_name r.d_s
        r.tc_result.d_f_values.(i) r.tc_result.tc_values.(i) comma
    done
  ) results;
  Printf.fprintf oc "  ]\n}\n"
