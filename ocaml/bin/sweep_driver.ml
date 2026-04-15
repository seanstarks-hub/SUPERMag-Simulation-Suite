(** CLI entry point for SUPERMag parameter sweep orchestration.

    Reads sweep parameters, dispatches to C++ solvers via FFI,
    and outputs results as CSV or JSON.

    Usage:
      supermag-sweep --solver proximity --param d_F --range 0.5,20.0,50 \
                     --sc Nb --fm Fe --output results.csv
      supermag-sweep --stack "Nb:30/Fe:8" --param d_F --range 0.5,20,50
      supermag-sweep --explore --param d_F --range 0.5,20,50 --format csv
      supermag-sweep --validate buzdin_1982 *)

open Cmdliner
open Supermag_types
open Supermag_pipeline

(* ── Output formatters ───────────────────────────────── *)

let output_csv oc param_name tc0 sweep_values tc_values =
  Printf.fprintf oc "# SUPERMag sweep: param=%s, Tc0=%.6f\n" param_name tc0;
  Printf.fprintf oc "%s,Tc_K\n" param_name;
  Array.iter2 (fun x tc ->
    Printf.fprintf oc "%.6f,%.6f\n" x tc
  ) sweep_values tc_values

let output_json oc param_name tc0 sweep_values tc_values =
  Printf.fprintf oc "{\n";
  Printf.fprintf oc "  \"sweep_param\": \"%s\",\n" param_name;
  Printf.fprintf oc "  \"tc0\": %.6f,\n" tc0;
  Printf.fprintf oc "  \"results\": [\n";
  let n = Array.length sweep_values in
  for i = 0 to n - 1 do
    let comma = if i < n - 1 then "," else "" in
    Printf.fprintf oc "    {\"%s\": %.6f, \"Tc_K\": %.6f}%s\n"
      param_name sweep_values.(i) tc_values.(i) comma
  done;
  Printf.fprintf oc "  ]\n";
  Printf.fprintf oc "}\n"

(* ── Parse --range min,max,n ─────────────────────────── *)

let parse_range s =
  match String.split_on_char ',' s with
  | [a; b; c] ->
    (try Ok (Float.of_string a, Float.of_string b, int_of_string c)
     with _ -> Error (`Msg (Printf.sprintf "bad --range format: %s" s)))
  | _ -> Error (`Msg (Printf.sprintf "--range must be min,max,n — got: %s" s))

(* ── Build params from SC + FM lookup ────────────────── *)

let build_params sc_name fm_name model_str phase_str gamma gamma_b d_s_opt geom_str =
  match Material.get_superconductor sc_name with
  | None -> Error (Printf.sprintf "unknown superconductor: %s" sc_name)
  | Some sc ->
    match Material.get_ferromagnet fm_name with
    | None -> Error (Printf.sprintf "unknown ferromagnet: %s" fm_name)
    | Some fm ->
      match Params.model_of_string model_str with
      | Error e -> Error e
      | Ok model ->
        match Params.phase_of_string phase_str with
        | Error e -> Error e
        | Ok phase ->
          match Params.geometry_of_string geom_str with
          | Error e -> Error e
          | Ok geometry ->
            let d_s = match d_s_opt with
              | Some v -> v
              | None -> sc.xi_s
            in
            Ok Params.{
              tc0 = sc.tc; d_s; d_s_coeff = 1.0; xi_s = sc.xi_s; xi_f = fm.xi_f;
              gamma; gamma_b; e_ex = fm.e_ex; d_f_coeff = fm.d_f;
              model; phase; geometry;
              geom_config = None; spin_active = None;
            }

(* ── Parse --depairing gamma_s,H,D,thickness,gamma_so,T ── *)

let parse_depairing s =
  if s = "" then Ok Params.no_depairing
  else
    match String.split_on_char ',' s with
    | [gs; h; d; th; gso; t] ->
      (try
        let gs = Float.of_string gs and h = Float.of_string h
        and d = Float.of_string d and th = Float.of_string th
        and gso = Float.of_string gso and t = Float.of_string t in
        let (ag, zeeman, orbital, spin_orbit) =
          Supermag_ffi.Stubs.depairing_from_physical
            ~gamma_s_mev:gs ~h_tesla:h ~d_nm2ps:d
            ~thickness_nm:th ~gamma_so_mev:gso ~t_kelvin:t in
        Ok Params.{ ag; zeeman; orbital; spin_orbit }
      with _ -> Error (`Msg (Printf.sprintf "bad --depairing format: %s" s)))
    | _ -> Error (`Msg (Printf.sprintf "--depairing must be 6 comma-separated values: %s" s))

(* ── Explore output ──────────────────────────────────── *)

let output_explore_csv oc (result : Result.exploration_result) =
  Printf.fprintf oc "# SUPERMag exploration: %d stacks evaluated\n"
    result.total_evaluated;
  Printf.fprintf oc "rank,description,tc_min_K,d_f_at_min_nm\n";
  List.iteri (fun i dr ->
    Printf.fprintf oc "%d,%s,%.6f,%.6f\n"
      (i + 1) dr.Result.description dr.tc_min dr.d_f_at_min
  ) result.ranked

let output_explore_json oc (result : Result.exploration_result) =
  Printf.fprintf oc "{\n";
  Printf.fprintf oc "  \"total_evaluated\": %d,\n" result.total_evaluated;
  Printf.fprintf oc "  \"ranked\": [\n";
  let n = List.length result.ranked in
  List.iteri (fun i dr ->
    let comma = if i < n - 1 then "," else "" in
    Printf.fprintf oc
      "    {\"rank\": %d, \"description\": \"%s\", \"tc_min_K\": %.6f, \"d_f_at_min_nm\": %.6f}%s\n"
      (i + 1) dr.Result.description dr.tc_min dr.d_f_at_min comma
  ) result.ranked;
  Printf.fprintf oc "  ]\n}\n"

(* ── Main sweep command ──────────────────────────────── *)

let run_sweep param_str range_str sc fm model phase
    gamma gamma_b format_str output_file d_f_range_str d_s_opt geom_str
    depairing_str stack_opt explore_flag =
  match Sweep.sweep_param_of_string param_str with
  | Error msg -> Printf.eprintf "Error: %s\n" msg; 1
  | Ok sweep_param ->
    match parse_range range_str with
    | Error (`Msg msg) -> Printf.eprintf "Error: %s\n" msg; 1
    | Ok (lo, hi, n) ->
      let d_f_array = Sweep.grid_sweep Sweep.{
        param_name = "d_F"; min_val = lo; max_val = hi;
        n_points = n; sweep_type = Grid;
      } in
      match parse_depairing depairing_str with
      | Error (`Msg msg) -> Printf.eprintf "Error: %s\n" msg; 1
      | Ok depairing ->
        let oc = match output_file with
          | "" -> stdout
          | path -> open_out path
        in
        let close_oc () = if output_file <> "" then close_out oc in
        if explore_flag then begin
          (* --explore: combinatorial SC×FM sweep *)
          let stacks = Design.enumerate_bilayers
            ~sc:Material.all_superconductors
            ~fm:Material.all_ferromagnets
            ~d_s:38.0 ~d_f:lo
          in
          let result = Design.explore ~stacks ~d_f_array ~depairing () in
          if format_str = "json" then
            output_explore_json oc result
          else
            output_explore_csv oc result;
          close_oc (); 0
        end else begin
          (* Resolve params: --stack takes priority over --sc/--fm *)
          let params_result = match stack_opt with
            | Some s ->
              (match Device.parse_stack s with
               | Error e -> Error e
               | Ok geom ->
                 match Device.resolve geom () with
                 | Error e -> Error e
                 | Ok p -> Ok Params.{ p with gamma; gamma_b })
            | None ->
              build_params sc fm model phase gamma gamma_b d_s_opt geom_str
          in
          match params_result with
          | Error msg -> Printf.eprintf "Error: %s\n" msg; 1
          | Ok params ->
            let param = sweep_param in
            let sweep_values = Sweep.grid_sweep Sweep.{
              param_name = param_str; min_val = lo; max_val = hi;
              n_points = n; sweep_type = Grid;
            } in
            let d_f_arr = match d_f_range_str with
              | "" -> [|1.0|]
              | s ->
                match parse_range s with
                | Ok (lo2, hi2, n2) ->
                  Sweep.grid_sweep Sweep.{
                    param_name = "d_F"; min_val = lo2; max_val = hi2;
                    n_points = n2; sweep_type = Grid;
                  }
                | Error _ -> [|1.0|]
            in
            match Sweep.tc_parameter_sweep param sweep_values
                ~params ~d_f_array:d_f_arr ~depairing () with
            | Error msg -> Printf.eprintf "Solver error: %s\n" msg; 1
            | Ok result ->
              let param_name = Sweep.sweep_param_to_string param in
              if format_str = "json" then
                output_json oc param_name result.tc0
                  result.d_f_values result.tc_values
              else
                output_csv oc param_name result.tc0
                  result.d_f_values result.tc_values;
              close_oc (); 0
        end

(* ── Cmdliner terms ──────────────────────────────────── *)

let param_t =
  let doc = "Parameter to sweep: d_F, d_S, D_S, gamma, gamma_B, E_ex, xi_F, Tc0" in
  Arg.(required & opt (some string) None & info ["param"; "p"] ~doc ~docv:"PARAM")

let range_t =
  let doc = "Sweep range as min,max,n_points" in
  Arg.(required & opt (some string) None & info ["range"; "r"] ~doc ~docv:"RANGE")

let sc_t =
  let doc = "Superconductor name (Nb, Pb, Al)" in
  Arg.(value & opt string "Nb" & info ["sc"] ~doc ~docv:"SC")

let fm_t =
  let doc = "Ferromagnet name (Fe, Co, Ni, Py, CuNi, Cu0.43Ni0.57)" in
  Arg.(value & opt string "Fe" & info ["fm"] ~doc ~docv:"FM")

let model_t =
  let doc = "Model: thin_s, fominov, or fominov_multi" in
  Arg.(value & opt string "thin_s" & info ["model"] ~doc ~docv:"MODEL")

let phase_t =
  let doc = "Phase: zero or pi" in
  Arg.(value & opt string "zero" & info ["phase"] ~doc ~docv:"PHASE")

let gamma_t =
  let doc = "Interface transparency parameter" in
  Arg.(value & opt float 0.3 & info ["gamma"] ~doc ~docv:"GAMMA")

let gamma_b_t =
  let doc = "Interface barrier parameter (Fominov model)" in
  Arg.(value & opt float 0.0 & info ["gamma-b"] ~doc ~docv:"GAMMA_B")

let format_t =
  let doc = "Output format: csv or json" in
  Arg.(value & opt string "csv" & info ["format"; "f"] ~doc ~docv:"FMT")

let output_t =
  let doc = "Output file path (default: stdout)" in
  Arg.(value & opt string "" & info ["output"; "o"] ~doc ~docv:"FILE")

let d_f_range_t =
  let doc = "d_F range for non-d_F sweeps, as min,max,n" in
  Arg.(value & opt string "" & info ["d-f-array"] ~doc ~docv:"D_F_RANGE")

let d_s_t =
  let doc = "Superconductor layer thickness d_S (nm). Defaults to xi_S if omitted." in
  Arg.(value & opt (some float) None & info ["d-s"] ~doc ~docv:"D_S")

let geometry_t =
  let doc = "Geometry: bilayer, trilayer, graded, or domains" in
  Arg.(value & opt string "bilayer" & info ["geometry"] ~doc ~docv:"GEOM")

let depairing_t =
  let doc = "Physical depairing inputs: gamma_s_meV,H_tesla,D_nm2ps,thickness_nm,Gamma_so_meV,T_K. \
             If omitted, no depairing is applied." in
  Arg.(value & opt string "" & info ["depairing"] ~doc ~docv:"DEPAIRING")

let stack_t =
  let doc = "Device stack in Nb:30/Fe:8 notation. Takes priority over --sc/--fm if supplied." in
  Arg.(value & opt (some string) None & info ["stack"] ~doc ~docv:"STACK")

let explore_t =
  let doc = "Run combinatorial SC×FM exploration over all material combinations." in
  Arg.(value & flag & info ["explore"] ~doc)

let sweep_cmd =
  let doc = "SUPERMag parameter sweep driver" in
  let info = Cmd.info "supermag-sweep" ~doc in
  Cmd.v info
    Term.(const run_sweep $ param_t $ range_t $ sc_t $ fm_t $ model_t $ phase_t
          $ gamma_t $ gamma_b_t $ format_t $ output_t $ d_f_range_t $ d_s_t
          $ geometry_t $ depairing_t $ stack_t $ explore_t)

let () = exit (Cmd.eval sweep_cmd)
