(** Geometry type: layer stacks, thicknesses, interfaces. *)

type layer =
  | S_layer of { thickness : float; sc : Material.superconductor }
  | F_layer of { thickness : float; fm : Material.ferromagnet }
  | I_layer of { thickness : float }  (** Insulator *)
  | N_layer of { thickness : float }  (** Normal metal *)

type geometry = {
  layers : layer list;
  description : string;
}

let bilayer ~s_thickness ~(sc : Material.superconductor)
    ~f_thickness ~(fm : Material.ferromagnet) : geometry =
  { layers = [
      S_layer { thickness = s_thickness; sc };
      F_layer { thickness = f_thickness; fm };
    ];
    description = Printf.sprintf "%s:%.4g/%s:%.4g" sc.name s_thickness fm.name f_thickness;
  }

let trilayer ~s_thickness ~(sc : Material.superconductor)
    ~n_thickness ~f_thickness ~(fm : Material.ferromagnet) : geometry =
  { layers = [
      S_layer { thickness = s_thickness; sc };
      N_layer { thickness = n_thickness };
      F_layer { thickness = f_thickness; fm };
    ];
    description = Printf.sprintf "%s:%.4g/N:%.4g/%s:%.4g"
      sc.name s_thickness n_thickness fm.name f_thickness;
  }

let total_thickness (g : geometry) : float =
  List.fold_left (fun acc layer ->
    acc +. match layer with
    | S_layer l -> l.thickness
    | F_layer l -> l.thickness
    | I_layer l -> l.thickness
    | N_layer l -> l.thickness
  ) 0.0 g.layers
