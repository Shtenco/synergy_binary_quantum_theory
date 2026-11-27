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

The selector side is isolated at `bqg_depth6_generic_jucys_selector.py`; its `expand(..., tol=0.0)` drops only exactly-zero coordinate segments. The Clebsch–Gordan recoupling module is isolated at `bqg_covariant_volume_recoupling.py`. A second hidden threshold was found in the S5 intertwiner action used by the symmetry-basis projector: the recovered `act_on_key` defaulted to `tol=1e-12`, and the recovered depth-4 transform builder calls it without an explicit tolerance. The derived `bqg_s5_intertwiner_action_gate.py` now defaults both `act_on_key` and `act_on_state` to zero, so `bqg_depth4_symmetric_master_builder._transform_matrix` resolves the no-prune action through normal module import resolution.

The recovered source directory remains byte-for-byte unchanged. The no-prune runner imports the isolated derived master, block operator, transfer, recoupling, selector and S5 action modules; it reuses recovered structural helpers such as spin/permutation layout, allowed labels, exact CG evaluation and oriented-spec enumeration. The derived S5 action is used by the recovered depth-6 basis module and depth-4 transform builder. Regression tests assert all these module paths.

## Threshold policy

- Operator amplitudes and sparse map entries are not removed by a positive magnitude threshold on the selected ultrafast path.
- Intermediate sparse accumulation receives `tol=0.0`; exact-zero results may be removed using elementwise nonzero checks, not a norm test.
- Gauss projection receives `tol=0.0`; no entrywise positive cutoff is applied, and retained-support detection is elementwise.
- The zero-tolerance branch avoids `np.linalg.norm(x) > 0` as a support test: squaring a value such as `1e-200` can underflow the norm to zero even though the amplitude is nonzero. A dedicated regression now checks this case in both block and master-row helpers.
- Selector expansion removes exact-zero coordinate segments only.
- S5 intertwiner action used by the selector's symmetry-basis projectors defaults to `tol=0.0`; a regression checks a five-vertex product of local recoupling coefficients that yields a nonzero `1e-15` global coefficient.
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
9. retention of `1e-200` amplitudes even when their squared Euclidean norm underflows to zero;
10. resolution of the S5 action module to the isolated derived copy;
11. retention of a `1e-15` S5 recoupling product through both `act_on_key` and the actual depth-4 `_transform_matrix` caller that omits an explicit tolerance; the recovered `1e-12` default would have deleted it;
12. end-to-end superposition linearity of the non-vacuum `master_map` output on a seeded two-input test.

On commit `38756a19da67c091c7a0248e4c4d0d3098e04ebf`, the canonical-core job reported the Neumaier/no-prune tests, Python compilation, core-scope audit, truth-ledger validation, finite-frontier validation, public status surface, LaTeX normalization, README local links and multiple mathematical control steps as successful. The full workflow was still in progress, with HDA and higher-shell jobs outstanding. The S5 no-prune action module and its regression were validated in the next core run. The end-to-end superposition-linearity regression was added after that run and still requires a fresh CI result.

## What this does not establish

- It does not prove that the isolated derived operator is mathematically identical to the intended canonical master operator.
- It does not establish chunking invariance, source-to-derived operator identity, or full dependency closure beyond the audited call path.
- It does not establish exact linearity in finite precision.
- It does not provide a complete all-orbit execution, rank result, nullspace result or non-vacuum injectivity certificate.
- It does not close C1. C1 remains OPEN pending full producer equivalence, complete finite rank verification in every S5 sector, and independent replayable evidence.

## Status

`NO_MAGNITUDE_PRUNING_PATH_IMPLEMENTED_FOR_THE_SELECTED_ULTRAFAST_PRODUCER`

`FINITE_FLOAT64_REGRESSION_TESTS_PASS_PRIOR_HEAD; MASTER_MAP_SUPERPOSITION_REGRESSION_PENDING_CURRENT_HEAD_CI`

`C1=OPEN`
