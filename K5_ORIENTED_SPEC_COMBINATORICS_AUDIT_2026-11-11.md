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
