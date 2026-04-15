(** Multi-parameter design optimizer with fabrication constraints.

    Composes existing solver primitives (solve_tc, usadel, josephson)
    into a closed-loop Nelder-Mead optimizer that jointly varies
    up to 5 parameters (d_S, d_F, γ, γ_B, E_ex) subject to
    fabrication constraints.

    No new C++ code — calls through the existing FFI. *)

open Supermag_types

(* ── Types ───────────────────────────────────────────── *)

(** Which parameters are free (optimized) vs. fixed.
    [Some (lo, hi)] means the parameter is free within that range. *)
type free_params = {
  vary_d_s : (float * float) option;
  vary_d_f : (float * float) option;
  vary_gamma : (float * float) option;
  vary_gamma_b : (float * float) option;
  vary_e_ex : (float * float) option;
}

let no_free = {
  vary_d_s = None; vary_d_f = None;
  vary_gamma = None; vary_gamma_b = None;
  vary_e_ex = None;
}

(** Fabrication constraints the optimizer must respect. *)
type fabrication_constraints = {
  d_s_range : (float * float) option;
  d_f_range : (float * float) option;
  max_total_thickness : float option;
  available_sc : Material.superconductor list option;
  available_fm : Material.ferromagnet list option;
}

let no_constraints = {
  d_s_range = None; d_f_range = None;
  max_total_thickness = None;
  available_sc = None; available_fm = None;
}

(** What the optimizer tries to achieve. *)
type objective =
  | Target_tc of float
  | Minimize_tc
  | Maximize_ic
  | Multi of (objective * float) list

(** Full optimization problem specification. *)
type problem = {
  base_params : Params.proximity_params;
  free : free_params;
  objective : objective;
  constraints : fabrication_constraints;
  depairing : Params.depairing;
  tolerance : float;
  max_evaluations : int;
}

(** Sensitivity: partial derivatives of Tc w.r.t. each free parameter. *)
type sensitivity = {
  d_tc_d_df : float;
  d_tc_d_ds : float;
  d_tc_d_gamma : float;
}

(** Optimization output. *)
type optimize_result = {
  optimal_params : Params.proximity_params;
  tc_achieved : float;
  d_f_optimal : float;
  evaluations : int;
  sensitivity : sensitivity option;
}

(* ── Free-parameter <-> vector mapping ───────────────── *)

(** Build the list of (getter, setter, lo, hi) for free params. *)
let free_param_specs (free : free_params) =
  let specs = ref [] in
  (match free.vary_e_ex with
   | Some (lo, hi) ->
     specs := (
       (fun (p : Params.proximity_params) -> p.e_ex),
       (fun (p : Params.proximity_params) v -> { p with e_ex = v }),
       lo, hi
     ) :: !specs
   | None -> ());
  (match free.vary_gamma_b with
   | Some (lo, hi) ->
     specs := (
       (fun (p : Params.proximity_params) -> p.gamma_b),
       (fun (p : Params.proximity_params) v -> { p with gamma_b = v }),
       lo, hi
     ) :: !specs
   | None -> ());
  (match free.vary_gamma with
   | Some (lo, hi) ->
     specs := (
       (fun (p : Params.proximity_params) -> p.gamma),
       (fun (p : Params.proximity_params) v -> { p with gamma = v }),
       lo, hi
     ) :: !specs
   | None -> ());
  (match free.vary_d_f with
   | Some (lo, hi) ->
     specs := (
       (fun _ -> 0.0),  (* d_f is in the sweep array, placeholder *)
       (fun p _ -> p),   (* applied separately *)
       lo, hi
     ) :: !specs
   | None -> ());
  (match free.vary_d_s with
   | Some (lo, hi) ->
     specs := (
       (fun (p : Params.proximity_params) -> p.d_s),
       (fun (p : Params.proximity_params) v -> { p with d_s = v }),
       lo, hi
     ) :: !specs
   | None -> ());
  Array.of_list !specs

(** Index of d_f in the specs array, if it's free. *)
let d_f_index (free : free_params) =
  let idx = ref 0 in
  let found = ref None in
  (match free.vary_d_s with Some _ -> incr idx | None -> ());
  (match free.vary_d_f with
   | Some _ -> found := Some !idx; incr idx
   | None -> ());
  ignore (free.vary_gamma, free.vary_gamma_b, free.vary_e_ex);
  !found

(** Extract the initial vector from base_params + free specs. *)
let params_to_vec (p : Params.proximity_params) (free : free_params)
    (specs : _ array) (d_f_init : float) =
  let n = Array.length specs in
  let vec = Array.make n 0.0 in
  let j = ref 0 in
  (match free.vary_d_s with
   | Some _ -> vec.(!j) <- p.d_s; incr j | None -> ());
  (match free.vary_d_f with
   | Some _ -> vec.(!j) <- d_f_init; incr j | None -> ());
  (match free.vary_gamma with
   | Some _ -> vec.(!j) <- p.gamma; incr j | None -> ());
  (match free.vary_gamma_b with
   | Some _ -> vec.(!j) <- p.gamma_b; incr j | None -> ());
  (match free.vary_e_ex with
   | Some _ -> vec.(!j) <- p.e_ex; incr j | None -> ());
  vec

(** Reconstruct params + d_f from a vector. *)
let vec_to_params (base : Params.proximity_params) (free : free_params)
    (specs : _ array) (vec : float array) : Params.proximity_params * float =
  let p = ref base in
  let d_f = ref 1.0 in
  let j = ref 0 in
  (match free.vary_d_s with
   | Some _ -> p := { !p with d_s = vec.(!j) }; incr j | None -> ());
  (match free.vary_d_f with
   | Some _ -> d_f := vec.(!j); incr j | None -> ());
  (match free.vary_gamma with
   | Some _ -> p := { !p with gamma = vec.(!j) }; incr j | None -> ());
  (match free.vary_gamma_b with
   | Some _ -> p := { !p with gamma_b = vec.(!j) }; incr j | None -> ());
  (match free.vary_e_ex with
   | Some _ -> p := { !p with e_ex = vec.(!j) }; incr j | None -> ());
  ignore specs;
  (!p, !d_f)

(* ── Constraint projection ───────────────────────────── *)

(** Clip each component of [vec] to its (lo, hi) bounds,
    then enforce max_total_thickness. *)
let project (specs : (_ * _ * float * float) array)
    (free : free_params) (constraints : fabrication_constraints)
    (vec : float array) =
  let n = Array.length vec in
  let out = Array.init n (fun i ->
    let (_, _, lo, hi) = specs.(i) in
    Float.max lo (Float.min hi vec.(i))
  ) in
  (* Enforce d_s + d_f <= max_total if both present *)
  (match constraints.max_total_thickness with
   | None -> ()
   | Some dmax ->
     let ds_idx = match free.vary_d_s with Some _ -> Some 0 | None -> None in
     let df_offset = match free.vary_d_s with Some _ -> 1 | None -> 0 in
     let df_idx = match free.vary_d_f with Some _ -> Some df_offset | None -> None in
     (match ds_idx, df_idx with
      | Some di, Some fi ->
        let total = out.(di) +. out.(fi) in
        if total > dmax then begin
          let scale = dmax /. total in
          out.(di) <- out.(di) *. scale;
          out.(fi) <- out.(fi) *. scale
        end
      | _ -> ()));
  out

(* ── Objective evaluation ────────────────────────────── *)

(** Compute the Tc at the minimum of the Tc(d_F) curve
    for a given parameter set. When d_F is free, the curve
    is evaluated at a single point; when it's fixed, a small
    sweep finds the minimum. *)
let eval_tc (params : Params.proximity_params) (d_f : float)
    (depairing : Params.depairing) (free : free_params) =
  let d_f_array = match free.vary_d_f with
    | Some (lo, hi) ->
      (* d_F is free — evaluate at the current d_F point *)
      [| d_f |]
    | None ->
      (* d_F is fixed — sweep to find minimum *)
      Array.init 50 (fun i ->
        0.5 +. Float.of_int i *. 19.5 /. 49.0)
  in
  match Supermag_ffi.Solvers.solve_tc ~params ~d_f_array ~depairing () with
  | Error msg -> Error msg
  | Ok r ->
    let tc_min = Array.fold_left Float.min Float.infinity r.tc_values in
    let idx = ref 0 in
    Array.iteri (fun i v -> if v = tc_min then idx := i) r.tc_values;
    Ok (tc_min, d_f_array.(!idx))

(** Evaluate the objective function.  Lower = better. *)
let rec eval_objective (obj : objective) (params : Params.proximity_params)
    (d_f : float) (depairing : Params.depairing) (free : free_params)
    (evals : int ref) : (float, string) result =
  incr evals;
  match obj with
  | Target_tc target ->
    (match eval_tc params d_f depairing free with
     | Error msg -> Error msg
     | Ok (tc_min, _) -> Ok (Float.abs (tc_min -. target)))
  | Minimize_tc ->
    (match eval_tc params d_f depairing free with
     | Error msg -> Error msg
     | Ok (tc_min, _) -> Ok tc_min)
  | Maximize_ic ->
    (* Usadel at the given d_F → spatial Δ(x), then Josephson CPR *)
    (match Supermag_ffi.Solvers.usadel ~tc0:params.tc0
             ~d_s:params.d_s ~d_f ~xi_s:params.xi_s ~xi_f:params.xi_f
             ~e_ex:params.e_ex ~t:(params.tc0 *. 0.5)
             ~mode:Params.Linearized ~n_grid:50 with
     | Error msg -> Error msg
     | Ok _profile ->
       (match Supermag_ffi.Solvers.josephson ~d_f ~xi_f:params.xi_f
                ~e_ex:params.e_ex ~t:(params.tc0 *. 0.5)
                ~gamma_b:params.gamma_b ~n_phases:64 with
        | Error msg -> Error msg
        | Ok cpr ->
          let max_i = Array.fold_left (fun acc v ->
            Float.max acc (Float.abs v)
          ) 0.0 cpr.current in
          Ok (-.max_i)))  (* Negate: want to maximize *)
  | Multi weighted ->
    let total = ref 0.0 in
    let err = ref None in
    List.iter (fun (sub_obj, weight) ->
      match !err with
      | Some _ -> ()
      | None ->
        match eval_objective sub_obj params d_f depairing free evals with
        | Error msg -> err := Some msg
        | Ok v -> total := !total +. weight *. v
    ) weighted;
    (match !err with Some msg -> Error msg | None -> Ok !total)

(* ── Nelder-Mead optimizer ───────────────────────────── *)

(** Initialize the simplex: one vertex at x0, others displaced
    along each axis by 5% of the axis range (or 1.0 if range=0). *)
let init_simplex (x0 : float array) (specs : (_ * _ * float * float) array) =
  let n = Array.length x0 in
  let simplex = Array.init (n + 1) (fun _ -> Array.copy x0) in
  for i = 0 to n - 1 do
    let (_, _, lo, hi) = specs.(i) in
    let delta = (hi -. lo) *. 0.05 in
    let delta = if Float.abs delta < 1e-12 then 1.0 else delta in
    simplex.(i + 1) <- Array.copy x0;
    simplex.(i + 1).(i) <- Float.min hi (x0.(i) +. delta)
  done;
  simplex

(** Centroid of all vertices except the worst (index [skip]). *)
let centroid (simplex : float array array) (skip : int) =
  let n = Array.length simplex.(0) in
  let m = Array.length simplex in
  let c = Array.make n 0.0 in
  for i = 0 to m - 1 do
    if i <> skip then
      for j = 0 to n - 1 do
        c.(j) <- c.(j) +. simplex.(i).(j)
      done
  done;
  let scale = 1.0 /. Float.of_int (m - 1) in
  Array.map (fun v -> v *. scale) c

(** Reflect point [worst] through [cent]: cent + alpha*(cent - worst) *)
let reflect (cent : float array) (worst : float array) (alpha : float) =
  Array.init (Array.length cent) (fun i ->
    cent.(i) +. alpha *. (cent.(i) -. worst.(i)))

(** Expand point: cent + gamma*(reflected - cent) *)
let expand (cent : float array) (refl : float array) (gamma : float) =
  Array.init (Array.length cent) (fun i ->
    cent.(i) +. gamma *. (refl.(i) -. cent.(i)))

(** Contract point: cent + rho*(worst - cent) *)
let contract (cent : float array) (worst : float array) (rho : float) =
  Array.init (Array.length cent) (fun i ->
    cent.(i) +. rho *. (worst.(i) -. cent.(i)))

(** Shrink all vertices toward the best. *)
let shrink (simplex : float array array) (best_idx : int) (sigma : float) =
  let best = simplex.(best_idx) in
  Array.iteri (fun i vertex ->
    if i <> best_idx then
      Array.iteri (fun j _ ->
        simplex.(i).(j) <- best.(j) +. sigma *. (vertex.(j) -. best.(j))
      ) vertex
  ) simplex

(** Main Nelder-Mead loop.
    Returns (best_vector, best_value, n_evaluations). *)
let nelder_mead ~(f : float array -> (float, string) result)
    ~(x0 : float array) ~(specs : (_ * _ * float * float) array)
    ~(free : free_params) ~(constraints : fabrication_constraints)
    ~(tol : float) ~(max_eval : int) =
  let n = Array.length x0 in
  if n = 0 then Error "no free parameters"
  else
    let simplex = init_simplex x0 specs in
    let proj v = project specs free constraints v in
    (* Project initial simplex *)
    Array.iteri (fun i v -> simplex.(i) <- proj v) simplex;
    (* Evaluate all vertices *)
    let values = Array.make (n + 1) Float.infinity in
    let evals = ref 0 in
    let err = ref None in
    Array.iteri (fun i v ->
      match !err with
      | Some _ -> ()
      | None ->
        match f v with
        | Error msg -> err := Some msg
        | Ok fv -> values.(i) <- fv; incr evals
    ) simplex;
    match !err with
    | Some msg -> Error msg
    | None ->
    (* Standard Nelder-Mead coefficients *)
    let alpha = 1.0 and gamma_nm = 2.0
    and rho = 0.5 and sigma = 0.5 in
    let converged = ref false in
    while !evals < max_eval && not !converged && !err = None do
      (* Sort indices by value *)
      let order = Array.init (n + 1) (fun i -> i) in
      Array.sort (fun a b -> Float.compare values.(a) values.(b)) order;
      let best_i = order.(0) in
      let worst_i = order.(n) in
      let second_worst_i = order.(n - 1) in
      (* Check convergence: range of values *)
      let spread = values.(worst_i) -. values.(best_i) in
      if spread < tol then
        converged := true
      else begin
        let cent = centroid simplex worst_i in
        let refl = proj (reflect cent simplex.(worst_i) alpha) in
        match f refl with
        | Error msg -> err := Some msg
        | Ok f_refl ->
          incr evals;
          if f_refl < values.(best_i) then begin
            (* Try expansion *)
            let exp = proj (expand cent refl gamma_nm) in
            match f exp with
            | Error msg -> err := Some msg
            | Ok f_exp ->
              incr evals;
              if f_exp < f_refl then begin
                simplex.(worst_i) <- exp;
                values.(worst_i) <- f_exp
              end else begin
                simplex.(worst_i) <- refl;
                values.(worst_i) <- f_refl
              end
          end else if f_refl < values.(second_worst_i) then begin
            simplex.(worst_i) <- refl;
            values.(worst_i) <- f_refl
          end else begin
            (* Contract *)
            let cont = proj (contract cent simplex.(worst_i) rho) in
            match f cont with
            | Error msg -> err := Some msg
            | Ok f_cont ->
              incr evals;
              if f_cont < values.(worst_i) then begin
                simplex.(worst_i) <- cont;
                values.(worst_i) <- f_cont
              end else begin
                (* Shrink *)
                shrink simplex best_i sigma;
                Array.iteri (fun i v ->
                  if i <> best_i then begin
                    simplex.(i) <- proj v;
                    match f simplex.(i) with
                    | Error msg -> err := Some msg
                    | Ok fv -> values.(i) <- fv; incr evals
                  end
                ) simplex
              end
          end
      end
    done;
    match !err with
    | Some msg -> Error msg
    | None ->
      (* Find best *)
      let best_i = ref 0 in
      for i = 1 to n do
        if values.(i) < values.(!best_i) then best_i := i
      done;
      Ok (simplex.(!best_i), values.(!best_i), !evals)

(* ── Public API ──────────────────────────────────────── *)

(** Single-objective optimize: find params that achieve the objective.
    Uses Nelder-Mead in the free-parameter subspace, calling
    solve_tc at each evaluation. Respects fabrication constraints
    by projecting infeasible points back to the boundary. *)
let optimize (prob : problem) : (optimize_result, string) result =
  let specs = free_param_specs prob.free in
  let n = Array.length specs in
  if n = 0 then
    (* No free params — just evaluate at base *)
    let evals = ref 0 in
    match eval_tc prob.base_params 1.0 prob.depairing prob.free with
    | Error msg -> Error msg
    | Ok (tc_min, d_f_min) ->
      Ok {
        optimal_params = prob.base_params;
        tc_achieved = tc_min;
        d_f_optimal = d_f_min;
        evaluations = 1;
        sensitivity = None;
      }
  else
    (* For 1D Target_tc with only d_F free, delegate to golden-section *)
    let use_golden = n = 1 && prob.free.vary_d_f <> None in
    match prob.objective, use_golden with
    | Target_tc target, true ->
      let (lo, hi) = match prob.free.vary_d_f with
        | Some r -> r | None -> assert false in
      (match Supermag_ffi.Solvers.optimize_tc ~params:prob.base_params
               ~d_f_lo:lo ~d_f_hi:hi ~tc_target:target
               ~depairing:prob.depairing () with
       | Error msg -> Error msg
       | Ok opt ->
         (* Evaluate Tc at the optimum *)
         (match eval_tc prob.base_params opt.d_f_optimal prob.depairing prob.free with
          | Error msg -> Error msg
          | Ok (tc_val, _) ->
            Ok {
              optimal_params = prob.base_params;
              tc_achieved = tc_val;
              d_f_optimal = opt.d_f_optimal;
              evaluations = 1;  (* golden-section internal count unknown *)
              sensitivity = None;
            }))
    | _ ->
      (* Multi-dimensional Nelder-Mead *)
      let d_f_init = match prob.free.vary_d_f with
        | Some (lo, hi) -> (lo +. hi) /. 2.0
        | None -> 5.0
      in
      let x0 = params_to_vec prob.base_params prob.free specs d_f_init in
      let evals = ref 0 in
      let f vec =
        let (params, d_f) = vec_to_params prob.base_params prob.free specs vec in
        eval_objective prob.objective params d_f prob.depairing prob.free evals
      in
      match nelder_mead ~f ~x0 ~specs ~free:prob.free
              ~constraints:prob.constraints
              ~tol:prob.tolerance ~max_eval:prob.max_evaluations with
      | Error msg -> Error msg
      | Ok (best_vec, _best_val, n_evals) ->
        let (opt_params, d_f_opt) =
          vec_to_params prob.base_params prob.free specs best_vec in
        (match eval_tc opt_params d_f_opt prob.depairing prob.free with
         | Error msg -> Error msg
         | Ok (tc_val, d_f_at_min) ->
           let sens = match sensitivity_at opt_params ~d_f:d_f_at_min
                              ~depairing:prob.depairing () with
             | Ok s -> Some s
             | Error _ -> None
           in
           Ok {
             optimal_params = opt_params;
             tc_achieved = tc_val;
             d_f_optimal = d_f_at_min;
             evaluations = n_evals;
             sensitivity = sens;
           })

(** Sensitivity analysis at a given operating point.
    Central finite differences: ∂Tc/∂p ≈ (Tc(p+h) - Tc(p-h)) / (2h). *)
and sensitivity_at (params : Params.proximity_params) ~(d_f : float)
    ?(depairing = Params.no_depairing) ?(step = 0.01) () =
  let eval_at p df =
    let d_f_arr = [| df |] in
    match Supermag_ffi.Solvers.solve_tc ~params:p ~d_f_array:d_f_arr ~depairing () with
    | Error msg -> Error msg
    | Ok r -> Ok r.tc_values.(0)
  in
  let h_df = Float.max 0.01 (d_f *. step) in
  let h_ds = Float.max 0.01 (params.d_s *. step) in
  let h_g = Float.max 0.001 (params.gamma *. step) in
  match eval_at params (d_f +. h_df), eval_at params (d_f -. h_df) with
  | Ok tp, Ok tm ->
    let dtc_ddf = (tp -. tm) /. (2.0 *. h_df) in
    let p_ds_p = { params with d_s = params.d_s +. h_ds } in
    let p_ds_m = { params with d_s = params.d_s -. h_ds } in
    (match eval_at p_ds_p d_f, eval_at p_ds_m d_f with
     | Ok tp2, Ok tm2 ->
       let dtc_dds = (tp2 -. tm2) /. (2.0 *. h_ds) in
       let p_g_p = { params with gamma = params.gamma +. h_g } in
       let p_g_m = { params with gamma = Float.max 0.0 (params.gamma -. h_g) } in
       (match eval_at p_g_p d_f, eval_at p_g_m d_f with
        | Ok tp3, Ok tm3 ->
          let dtc_dg = (tp3 -. tm3) /. (2.0 *. h_g) in
          Ok { d_tc_d_df = dtc_ddf; d_tc_d_ds = dtc_dds; d_tc_d_gamma = dtc_dg }
        | Error msg, _ | _, Error msg -> Error msg)
     | Error msg, _ | _, Error msg -> Error msg)
  | Error msg, _ | _, Error msg -> Error msg

(** Design-for-manufacturing: find the parameter set that
    hits the objective AND is most robust to fabrication tolerances.
    Modifies the objective with a sensitivity penalty:
    f_robust(p) = f(p) + λ · Σ_i |∂Tc/∂p_i| · σ_i *)
let robust_optimize (prob : problem) ~(tolerances : float array)
    : (optimize_result, string) result =
  let specs = free_param_specs prob.free in
  let n = Array.length specs in
  if n = 0 then optimize prob
  else if Array.length tolerances <> n then
    Error (Printf.sprintf
             "tolerances length (%d) must match number of free params (%d)"
             (Array.length tolerances) n)
  else
    let lambda = 0.5 in  (* Penalty weight *)
    let d_f_init = match prob.free.vary_d_f with
      | Some (lo, hi) -> (lo +. hi) /. 2.0
      | None -> 5.0
    in
    let x0 = params_to_vec prob.base_params prob.free specs d_f_init in
    let evals = ref 0 in
    let f vec =
      let (params, d_f) = vec_to_params prob.base_params prob.free specs vec in
      match eval_objective prob.objective params d_f prob.depairing prob.free evals with
      | Error msg -> Error msg
      | Ok base_cost ->
        (* Add sensitivity penalty: penalize designs where Tc is sensitive
           to parameter drift within fabrication tolerances *)
        let penalty = ref 0.0 in
        let pen_err = ref None in
        for i = 0 to n - 1 do
          if !pen_err = None && tolerances.(i) > 0.0 then begin
            let h = tolerances.(i) *. 0.1 in  (* Small step for finite diff *)
            let v_plus = Array.copy vec in
            let v_minus = Array.copy vec in
            let (_, _, lo, hi) = specs.(i) in
            v_plus.(i) <- Float.min hi (vec.(i) +. h);
            v_minus.(i) <- Float.max lo (vec.(i) -. h);
            let (p_plus, df_plus) =
              vec_to_params prob.base_params prob.free specs v_plus in
            let (p_minus, df_minus) =
              vec_to_params prob.base_params prob.free specs v_minus in
            match eval_tc p_plus df_plus prob.depairing prob.free,
                  eval_tc p_minus df_minus prob.depairing prob.free with
            | Ok (tc_p, _), Ok (tc_m, _) ->
              let grad = Float.abs ((tc_p -. tc_m) /. (2.0 *. h)) in
              penalty := !penalty +. grad *. tolerances.(i)
            | Error msg, _ | _, Error msg ->
              pen_err := Some msg
          end
        done;
        (match !pen_err with
         | Some msg -> Error msg
         | None -> Ok (base_cost +. lambda *. !penalty))
    in
    match nelder_mead ~f ~x0 ~specs ~free:prob.free
            ~constraints:prob.constraints
            ~tol:prob.tolerance ~max_eval:prob.max_evaluations with
    | Error msg -> Error msg
    | Ok (best_vec, _best_val, n_evals) ->
      let (opt_params, d_f_opt) =
        vec_to_params prob.base_params prob.free specs best_vec in
      (match eval_tc opt_params d_f_opt prob.depairing prob.free with
       | Error msg -> Error msg
       | Ok (tc_val, d_f_at_min) ->
         let sens = match sensitivity_at opt_params ~d_f:d_f_at_min
                            ~depairing:prob.depairing () with
           | Ok s -> Some s
           | Error _ -> None
         in
         Ok {
           optimal_params = opt_params;
           tc_achieved = tc_val;
           d_f_optimal = d_f_at_min;
           evaluations = n_evals;
           sensitivity = sens;
         })
