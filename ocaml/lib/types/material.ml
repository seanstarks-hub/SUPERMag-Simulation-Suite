(** Material type: superconductors and ferromagnets with physical parameters. *)

type superconductor = {
  name : string;
  tc : float;          (** Critical temperature (K) *)
  xi_s : float;        (** Coherence length (nm) *)
  lambda_l : float;    (** London penetration depth (nm) *)
  delta_0 : float;     (** BCS gap at T=0 (meV) *)
}

type ferromagnet = {
  name : string;
  e_ex : float;        (** Exchange energy (meV) *)
  xi_f : float;        (** Coherence length in F (nm) *)
  d_f : float;         (** Diffusion coefficient (m^2/s) *)
}

type material = SC of superconductor | FM of ferromagnet

let nb : superconductor =
  { name = "Nb"; tc = 9.2; xi_s = 38.0; lambda_l = 39.0; delta_0 = 1.55 }

let fe : ferromagnet =
  { name = "Fe"; e_ex = 256.0; xi_f = 0.7; d_f = 2.5e-4 }
