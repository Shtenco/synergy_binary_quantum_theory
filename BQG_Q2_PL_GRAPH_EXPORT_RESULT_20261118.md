# q=2 PL graph export — implementation result (2026-11-18)

## What was implemented

Added `scripts/q2_pl_graph_export.py`. It starts from the same 16-cell boundary
construction used by `bcqg_global_manifold_gate.py`, then applies the same
global barycentric subdivision. It exports the 1-skeleton at levels 0, 1 and 2.

At refinement level n>0, each new graph vertex is canonically identified with
one non-empty simplex of the preceding PL complex. Two such vertices are
adjacent exactly when their parent simplices are strictly comparable by
inclusion. The exporter records this parent-simplex map, the full adjacency,
degree extrema, connected-component count, and a SHA-256 hash of canonical
vertex/edge data.

Added `tests/test_q2_pl_graph_export.py` and a dedicated workflow
`.github/workflows/q2-pl-graph-export.yml`.

## Independent combinatorial count check

The generated construction is expected to reproduce the independently
published PL refinement counts:

| PL level | graph vertices | graph edges | tetrahedra in complex |
|---:|---:|---:|---:|
| 0 | 8 | 24 | 16 |
| 1 | 80 | 464 | 384 |
| 2 | 1696 | 10912 | 9216 |

A separate in-session Python reconstruction of the same face-poset algorithm
reproduced all six vertex/edge counts exactly. GitHub Actions run
`44269971966` completed successfully: all 3 exporter tests passed, the
three-level JSON export completed with `passed: true`, and the 70,352-byte
`q2-pl-graph-export` artifact was uploaded. The run was triggered at commit
`fdab36fc27bc04801e98a252d5ae172b185c7631`; the subsequent report-only commit
did not alter the exporter or tests.

## Critical distinction — what this does not solve

This is the 1-skeleton of the selected PL completion and its barycentric
refinements. It is not the separate causal edge-rewrite graph whose scalar
vertex count obeys

`N_g=(4\,8^g+10)/7`.

The causal rewrite currently does not export its own canonical adjacency and
gluing maps. Therefore the exporter closes a real, useful PL-graph data gap,
but it does not yet provide the common graph identity required to claim that
the Q2 causal-volume exponent and QFI spectral dimension were measured on the
same microscopic graph.

Likewise, the PL gate does not supply physical reduced density matrices
`rho_e` or edge generators `K_e`. The exporter deliberately writes
`quantum_state_assignment: null` and `qfi_edge_conductances: null`; it does
not fabricate QFI values. A QFI heat-kernel calculation on this graph remains
blocked until those physical state/operator assignments are derived.

## Decision

- **Closed in this increment:** deterministic adjacency construction for the
  canonical 16-cell PL completion and two existing barycentric refinements;
  explicit refinement labels; graph hashes; deterministic tests and CI entry.
- **Still open:** canonical adjacency export for the distinct causal rewrite;
  proof that the microscopic dynamics selects the chosen PL gluing; physical
  edge-state/generator construction; QFI spectral plateau and topology
  controls on the same graph family.
- **No new physical dimension claim:** the exact causal-volume fixed point
  remains a result about the frozen rewrite, and the PL graph export does not
  upgrade it to a QFI spectral-dimension derivation.
