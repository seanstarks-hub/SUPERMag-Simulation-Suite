(** Material type: superconductors and ferromagnets with physical parameters. *)

type superconductor = {
  name : string;
  tc : float;          (** Critical temperature (K) *)
  xi_s : float;        (** Coherence length (nm) *)
  lambda_l : float;    (** London penetration depth (nm) *)
  delta_0 : float;     (** BCS gap at T=0 (meV) *)
  rho : float;         (** Normal-state resistivity (Ohm·nm) *)
  d_s : float;         (** Diffusion coefficient (nm^2/ps) *)
}

type ferromagnet = {
  name : string;
  e_ex : float;        (** Exchange energy (meV) *)
  xi_f : float;        (** Coherence length in F (nm) *)
  d_f : float;         (** Diffusion coefficient (m^2/s) *)
  rho : float;         (** Normal-state resistivity (Ohm·nm) *)
}

type material = SC of superconductor | FM of ferromagnet

(* ── Superconductors ─────────────────────────────────── *)

let nb : superconductor =
  { name = "Nb"; tc = 9.2; xi_s = 38.0; lambda_l = 39.0; delta_0 = 1.55;
    rho = 150.0; d_s = 18.0 }

let pb : superconductor =
  { name = "Pb"; tc = 7.2; xi_s = 83.0; lambda_l = 37.0; delta_0 = 1.35;
    rho = 220.0; d_s = 60.0 }

let al : superconductor =
  { name = "Al"; tc = 1.2; xi_s = 1600.0; lambda_l = 16.0; delta_0 = 0.18;
    rho = 27.0; d_s = 4500.0 }

let all_superconductors = [nb; pb; al]

(* ── Ferromagnets ────────────────────────────────────── *)

let fe : ferromagnet =
  { name = "Fe"; e_ex = 256.0; xi_f = 0.7; d_f = 2.5e-4; rho = 100.0 }

let co : ferromagnet =
  { name = "Co"; e_ex = 309.0; xi_f = 0.5; d_f = 1.8e-4; rho = 63.0 }

let ni : ferromagnet =
  { name = "Ni"; e_ex = 75.0; xi_f = 2.3; d_f = 5.0e-4; rho = 69.0 }

let py : ferromagnet =
  { name = "Py"; e_ex = 20.0; xi_f = 5.0; d_f = 3.0e-4; rho = 400.0 }

let cuni : ferromagnet =
  { name = "CuNi"; e_ex = 5.0; xi_f = 10.0; d_f = 4.0e-4; rho = 350.0 }

let cu043ni057 : ferromagnet =
  { name = "Cu0.43Ni0.57"; e_ex = 11.2; xi_f = 4.2; d_f = 4.0e-4; rho = 500.0 }

let all_ferromagnets = [fe; co; ni; py; cuni; cu043ni057]

(* ── Lookup ──────────────────────────────────────────── *)

let get_superconductor name =
  List.find_opt (fun s -> s.name = name) all_superconductors

let get_ferromagnet name =
  List.find_opt (fun f -> f.name = name) all_ferromagnets
