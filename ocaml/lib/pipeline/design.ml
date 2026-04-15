(** Combinatorial device stack explorer.

    Enumerates SC × FM material combinations, evaluates Tc curves
    for each, ranks by figure of merit (e.g., deepest Tc suppression),
    and optionally filters by constraints. *)

open Supermag_types
open Supermag_ffi

(** Constraints for filtering exploration results. *)
type constraints = {
  tc_min : float option;              (** Minimum allowed Tc (K) — reject stacks below this *)
  tc_max : float option;              (** Maximum allowed Tc (K) — reject stacks above this *)
  max_total_thickness : float option; (** Maximum total stack thickness (nm) *)
}

let no_constraints = { tc_min = None; tc_max = None; max_total_thickness = None }

(** Enumerate all bilayer combinations from SC and FM material lists.
    Returns geometry list with one bilayer per SC×FM pair.
    [d_s] is the S-layer thickness used for all combinations.
    [d_f] is a nominal F-layer thickness for the geometry description. *)
let enumerate_bilayers ~(sc : Material.superconductor list)
    ~(fm : Material.ferromagnet list)
    ~d_s ~d_f : Geometry.geometry list =
  List.concat_map (fun s ->
    List.map (fun f ->
      Geometry.bilayer ~s_thickness:d_s ~sc:s ~f_thickness:d_f ~fm:f
    ) fm
  ) sc

(** Evaluate a single device stack: resolve geometry to params,
    compute Tc curve, extract minimum.
    [gamma] and [gamma_b] default to 1.0 and 0.0 respectively. *)
let evaluate_stack (geom : Geometry.geometry) ~d_f_array
    ?(gamma = 1.0) ?(gamma_b = 0.0)
    ?(depairing = Params.no_depairing) ()
    : (Result.device_result, string) result =
  match Device.resolve geom () with
  | Error e -> Error e
  | Ok params ->
    let params = Params.{ params with gamma; gamma_b } in
    match Solvers.solve_tc ~params ~d_f_array ~depairing () with
    | Error e -> Error e
    | Ok tc ->
      let n = Array.length tc.tc_values in
      if n = 0 then Error "Tc curve returned no points"
      else begin
        let idx = ref 0 in
        for i = 1 to n - 1 do
          if tc.tc_values.(i) < tc.tc_values.(!idx) then idx := i
        done;
        Ok Result.{
          description = geom.description;
          tc;
          tc_min = tc.tc_values.(!idx);
          d_f_at_min = tc.d_f_values.(!idx);
        }
      end

(** Apply constraint filters to a device result. *)
let passes_constraints (c : constraints) (geom : Geometry.geometry)
    (dr : Result.device_result) : bool =
  (match c.tc_min with
   | Some lo -> dr.tc_min >= lo
   | None -> true)
  &&
  (match c.tc_max with
   | Some hi -> dr.tc_min <= hi
   | None -> true)
  &&
  (match c.max_total_thickness with
   | Some max_t -> Geometry.total_thickness geom <= max_t
   | None -> true)

(** Run full combinatorial exploration.
    Evaluates all stacks, optionally in parallel via Domain.spawn,
    filters by constraints, ranks by deepest Tc suppression (lowest tc_min).
    [parallel]: if true, uses Domain.spawn for parallelism (default: false). *)
let explore ~(stacks : Geometry.geometry list) ~d_f_array
    ?(gamma = 1.0) ?(gamma_b = 0.0)
    ?(depairing = Params.no_depairing)
    ?(constraints = no_constraints)
    ?(parallel = false)
    () : Result.exploration_result =
  let evaluate geom =
    evaluate_stack geom ~d_f_array ~gamma ~gamma_b ~depairing ()
  in
  let results =
    if parallel then begin
      let domains = List.map (fun geom ->
        Domain.spawn (fun () -> evaluate geom)
      ) stacks in
      List.combine stacks (List.map Domain.join domains)
    end else
      List.map (fun geom -> (geom, evaluate geom)) stacks
  in
  let passing = List.filter_map (fun (geom, r) ->
    match r with
    | Ok dr when passes_constraints constraints geom dr -> Some dr
    | _ -> None
  ) results in
  let ranked = List.sort (fun a b ->
    compare a.Result.tc_min b.Result.tc_min
  ) passing in
  Result.{ ranked; total_evaluated = List.length stacks }
