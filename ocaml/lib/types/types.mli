(** Shared type signatures for SUPERMag OCaml types. *)

type superconductor = {
  name : string;
  tc : float;
  xi_s : float;
  lambda_l : float;
  delta_0 : float;
}

type ferromagnet = {
  name : string;
  e_ex : float;
  xi_f : float;
  d_f : float;
}

type material = SC of superconductor | FM of ferromagnet
