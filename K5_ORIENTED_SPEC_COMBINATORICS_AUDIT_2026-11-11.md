# K5 oriented-spec combinatorics audit — 2026-11-11

## Purpose

This note records a new source-level regression gate for the orientation generator used by the production K5 node-H construction. It is deliberately separate from operator-adjoint tests: a correct primitive adjoint does not establish that the full set of oriented terms is enumerated with the intended signs and multiplicities.

## Definition under test

For a vertex (v), let its four ordered neighbors be (n_0,n_1,n_2,n_3). For each omitted position (r\in\{0,1,2,3\}), retain the other three neighbors in their inherited order as ((a,b,c)). The generator must emit exactly these three cyclic orderings:

- (+(v,a,b,c)), ((v,b,c,a)), ((v,c,a,b)) with common sign ((-1)^r).

Thus there are exactly (4\times3=12) oriented specifications per vertex.

## Committed test

Added `test_all_twelve_oriented_specs_have_exact_combinatorics` to `scripts/test_k5_primitive_dense_adjoint.py`.

The test:
1. constructs the expected sequence independently from the documented omitted-neighbor/cyclic-order rule;
2. compares the exact ordered output at vertex 0;
3. checks count and uniqueness;
4. repeats the independent construction and exact comparison for all five K5 vertices.

This is a combinatorial generator test. It does not prove covariance of the resulting quantum operator, adjointness of the coherent T sum, or HDA closure.

## Evidence status

- Test source committed on `research/consolidated-bqg-audit-20261111`.
- Test method inventory increases from 39 to 40 across six workflow scripts.
- GitHub Actions execution is pending; no PASS is claimed until the job completes.
- The research branch remains separate from `main`; no merge or rebase is performed.

## Next operator-level gates

1. Matrix-element adjointness of the complete 32-term T sum, not only one monomial.
2. Matrix-element Hermiticity of assembled (H_v=(T_v+T_v^\dagger)/2), with the full 12-spec signed sum.
3. Spin-wall stability with explicit accounting of transitions clipped at each finite wall.
4. Only after these gates: full HH composition and leakage Gram (L^\dagger L), retaining off-diagonal interference.

## Full coherent T-sum test added — 2026-11-11

The production test file now includes `test_full_32_term_T_sum_matches_reverse_matrix_element`. Unlike the five-hit single-monomial test, this evaluates all 32 signed monomials, coherently sums their contributions, and checks the reverse matrix element from the complete adjoint sum. It uses one deterministic largest-magnitude nonzero output element from a fixed all-spin-1/all-(K=0) input at doubled-spin wall 2.

This is a materially stronger local acceptance gate, but not a full matrix adjoint proof. Its scope remains one matrix element, one input, one orientation triple, and one finite wall. The production primitive/orientation script now has five methods; total workflow source inventory is 41 methods. CI execution remains unconfirmed.


## CI and next-scope clarification — 2026-11-11

The full 32-term test is committed in code commit `ba4dfc48db4f1a2ea5caea8ef7a3160fc63bdd00`; its Actions run #28 (`42944740425`) was queued at the last observation. Runs #26 and #27 were also queued. No PASS is claimed.

The test checks one selected matrix element of the complete coherent T sum. It is not a matrix-wide adjoint certificate. The next target is an explicit finite-basis matrix comparison over multiple input columns and orientation triples, with spin-wall clipping recorded rather than hidden.

## Local oriented-volume self-adjointness gate added

The volume suite now includes a direct test of the vertex-oriented tensor action under the Frobenius inner product. It checks all five vertices and all four spin triples, with deterministic complex test tensors, and therefore exercises the production epsilon conjugations before and after the local spectral volume map.

This expands the volume test file to seven methods and the six-script source inventory to 42. The test is committed in `50c807d03c8616dce639487afc83704ec19c5b36`. It remains unexecuted until the corresponding Actions job completes.


## Matrix-wide T/T† continuation checkpoint — 2026-11-11

The latest production audit now records 42 source-inventoried test methods across six workflow scripts. The 32-term T check is one selected matrix element only. The next required gate is a rectangular matrix comparison over multiple orthonormal input columns, with coherent no-prune accumulation of all 32 signed terms and independently generated reverse sequences. Forward outputs outside the declared codomain and reverse outputs outside the declared domain must be counted and reported; no silent projection or thresholding is permitted.

The assembled node operator `H_v=(T_v+T_v†)/2` is a separate subsequent gate over all 12 signed orientation specifications. A successful H test cannot substitute for checking the constituent T/T† relation. Actions runs #26–#29 were all observed queued with no job conclusion; no CI PASS is claimed. The branch is 109 commits ahead and 2 behind `main`, with no merge/rebase performed.
