# C1 no-magnitude-pruning derivation audit — 2026-11-27

## Purpose

This note records the isolated depth-6 C1 producer path created for a controlled comparison against the recovered thresholded producer. It is an implementation audit, not a C1 proof or a rank certificate.

## Isolated dependency path

The shard runner at `research/depth6_no_prune/github_mixed_master_shard.py` imports the isolated generic master-row module. The path is:

`github_mixed_master_shard.work`
→ `bqg_depth6_generic_master_row.master_map(..., tol=0.0)`
→ `reduced_s3_node / reduced_s2_node`
→ `bqg_block_apply_invariants.H_node_apply_ultrafast(..., tol=0.0)`
→ `initial_reduced_block(..., tol=0.0)`
→ `oriented_T_ultrafast(..., tol=0.0)`
→ `apply_radial_collapsed_block / triangle_forward_collapsed / triangle_adjoint_collapsed`
→ `local_two_hit_outputs`
→ `bqg_reduced_transfer_hamiltonian_no_prune.local_hit_outputs(..., tol=0.0)`.

The selector side is isolated at `bqg_depth6_generic_jucys_selector.py`; its `expand(..., tol=0.0)` now drops only exactly-zero coordinate segments. The Clebsch–Gordan recoupling module is isolated at `bqg_covariant_volume_recoupling.py` and retains every coefficient except exact structural zeros.

The recovered source directory is kept unchanged and is not imported as the no-prune operator implementation. Shared imports from the recovered tree are limited to helpers used on this path that do not perform magnitude pruning: spin/permutation structure, allowed labels, exact CG evaluation, and the oriented-spec enumeration. The no-prune regression test asserts the resolved module paths for the master, block operator, transfer kernel, recoupling module and selector.

## Threshold policy

- Operator amplitudes and sparse map entries are not removed by a positive magnitude threshold on the selected ultrafast path.
- Intermediate sparse accumulation receives `tol=0.0`; exact-zero results may be removed using elementwise nonzero checks, not a norm test.
- Gauss projection receives `tol=0.0`; no entrywise positive cutoff is applied, and retained-support detection is elementwise.
- The zero-tolerance branch avoids `np.linalg.norm(x) > 0` as a support test: squaring a value such as `1e-200` can underflow the norm to zero even though the amplitude is nonzero. A dedicated regression now checks this case in both block and master-row helpers.
- Selector expansion removes exact-zero coordinate segments only.
- The selector's QR/SVD rank decision still uses a machine-precision-scaled rank tolerance and validates the expected multiplicity. This is a numerical rank diagnostic, not an operator-amplitude pruning threshold. The recovered `exact_mult` still validates floating character traces against `1e-8`/`1e-7` rounding tolerances; this is a discrete multiplicity-estimation dependency and remains a numerical risk, not an amplitude filter.
- Floating-point underflow, roundoff, conditioning and finite-precision rank decisions remain. This is not exact rational/algebraic arithmetic.
- The copied block module also retains legacy factorized/fast entry points with their historical thresholds. The shard runner does not call those entry points; it calls only `H_node_apply_ultrafast`. Do not describe the whole module as globally threshold-free.

## Regression coverage

`tests/test_bqg_depth6_no_prune_mode.py` covers:
1. preservation of a sub-threshold input row when `tol=0`;
2. order-dependence of the legacy thresholded sparse accumulator versus the zero-tolerance path;
3. preservation through zero-tolerance Gauss projection;
4. propagation of zero tolerance through the node wrapper;
5. import resolution to the isolated derived modules;
6. shard-runner wiring and explicit `master_map(..., tol=0.0)`;
7. retention of sub-threshold CG products;
8. retention of tiny selector coordinates while omitting an exactly-zero segment;
9. retention of `1e-200` amplitudes even when their squared Euclidean norm underflows to zero.

The prior CI run on commit `e7e22fc7088df017f6b526dc3f869fa08a91b011` reported the Neumaier/no-prune tests, Python compilation, core-scope audit, truth-ledger validation and finite-frontier validation as successful. That run failed later at the public README LaTeX gate; the README delimiter issue was corrected in a subsequent commit. A new CI run is required to certify the current head.

## What this does not establish

- It does not prove that the isolated derived operator is mathematically identical to the intended canonical master operator.
- It does not establish chunking invariance, source-to-derived operator identity, or full dependency closure beyond the audited call path.
- It does not establish exact linearity in finite precision.
- It does not provide a complete all-orbit execution, rank result, nullspace result or non-vacuum injectivity certificate.
- It does not close C1. C1 remains OPEN pending full producer equivalence, complete finite rank verification in every S5 sector, and independent replayable evidence.

## Status

`NO_MAGNITUDE_PRUNING_PATH_IMPLEMENTED_FOR_THE_SELECTED_ULTRAFAST_PRODUCER`

`FINITE_FLOAT64_REGRESSION_TESTS_PENDING_CURRENT_HEAD_CI`

`C1=OPEN`
