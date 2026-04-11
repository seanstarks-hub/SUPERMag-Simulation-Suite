(** Tests for solver chaining. *)

open Supermag_pipeline

let test_chain_identity () =
  let id_step = Chain.{ name = "identity"; run = Fun.id } in
  let result = Chain.chain [id_step; id_step] 42 in
  Alcotest.(check int) "identity chain" 42 result

let test_chain_increment () =
  let inc = Chain.{ name = "inc"; run = (fun x -> x + 1) } in
  let result = Chain.chain [inc; inc; inc] 0 in
  Alcotest.(check int) "three increments" 3 result

let test_chain_compose () =
  let double = Chain.{ name = "double"; run = (fun x -> x * 2) } in
  let add_one = Chain.{ name = "add_one"; run = (fun x -> x + 1) } in
  let result = Chain.chain [double; add_one; double] 3 in
  (* 3 -> 6 -> 7 -> 14 *)
  Alcotest.(check int) "compose" 14 result

let test_chain_empty () =
  let result = Chain.chain [] 99 in
  Alcotest.(check int) "empty chain" 99 result

let test_parallel_chain () =
  let double = Chain.{ name = "double"; run = (fun x -> x * 2) } in
  let triple = Chain.{ name = "triple"; run = (fun x -> x * 3) } in
  let results = Chain.parallel_chain [double; triple] [10; 10] in
  Alcotest.(check (list int)) "parallel" [20; 30] results

let test_parallel_chain_single () =
  let inc = Chain.{ name = "inc"; run = (fun x -> x + 1) } in
  let results = Chain.parallel_chain [inc] [5] in
  Alcotest.(check (list int)) "single parallel" [6] results

let () =
  Alcotest.run "Chain" [
    "sequential", [
      Alcotest.test_case "identity"  `Quick test_chain_identity;
      Alcotest.test_case "increment" `Quick test_chain_increment;
      Alcotest.test_case "compose"   `Quick test_chain_compose;
      Alcotest.test_case "empty"     `Quick test_chain_empty;
    ];
    "parallel", [
      Alcotest.test_case "basic"  `Quick test_parallel_chain;
      Alcotest.test_case "single" `Quick test_parallel_chain_single;
    ];
  ]
