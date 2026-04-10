(** Geometry type: layer stacks, thicknesses, interfaces. *)

type layer =
  | S_layer of { thickness : float; material : string }
  | F_layer of { thickness : float; material : string }
  | I_layer of { thickness : float }  (** Insulator *)
  | N_layer of { thickness : float }  (** Normal metal *)

type geometry = {
  layers : layer list;
  description : string;
}

let bilayer ~s_thickness ~s_material ~f_thickness ~f_material : geometry =
  { layers = [
      S_layer { thickness = s_thickness; material = s_material };
      F_layer { thickness = f_thickness; material = f_material };
    ];
    description = Printf.sprintf "S(%s)/F(%s) bilayer" s_material f_material;
  }

let total_thickness (g : geometry) : float =
  List.fold_left (fun acc layer ->
    acc +. match layer with
    | S_layer l -> l.thickness
    | F_layer l -> l.thickness
    | I_layer l -> l.thickness
    | N_layer l -> l.thickness
  ) 0.0 g.layers
