# BQG C1/C2 threshold and CI audit — 2026-11-27

## Executive status

- C1: **OPEN**. The depth-6 sector census is arithmetically consistent, but nonvacuum master-map rank certificates are still missing.
- C2: **OPEN**. The Neumaier sparse accumulator passes primitive tests, but the production W5 dependency graph still contains multiple state-dependent amplitude cutoffs and ordinary floating-point accumulation sites.
- No production W5 PASS, complete Lorentzian operator certificate, or theory closure is claimed.

This audit was performed against PR branch `audit/c1-c2-2026-11-27`, initially at head `65ae25e9a9927f82b3b99e623e9bd8b4dfd4b9a3`. Source blob hashes are listed below so the inventory can be reproduced.

## 1. C1: rank-producer audit

### 1.1 Recovered producer

File: `research/depth6_recovered_source/github_mixed_master_shard.py`
Blob: `da4f155517872b98b122da13a382bd8ceb2938ec`

The shard producer constructs the depth-6 shell, obtains S5 orbit representatives, builds active blocks and assigns them to cost-balanced shards. It calls:

```python
W = J.fast_selector(tuple(rep), _KEY)
mp0, _ = M.master_map(tuple(rep), _KEY, W)
```

The recovered shard script records assignment/column counts and per-shard outputs. It does not itself establish global coverage or rank; those require complete, validated shard inventory and a deterministic aggregation step.

### 1.2 Nested operator cutoffs (important)

The generic master-row tolerance is not the only cutoff. Its dependency `bqg_block_apply_invariants.py` (blob `6f5ded2f4af075adba46bdb0d8a4fbf512fa07c3`) contains additional fixed thresholds:
- `initial_reduced_block` discards input rows with norm at or below the module constant `TOL=1e-12`;
- `add_vec` uses a default norm cutoff of `1e-12`;
- `project_gauss_block` defaults to `1e-10` and zeroes individual matrix entries below that value before a block-norm test;
- `H_node_apply_ultrafast` calls `oriented_T_ultrafast` with a hard-coded `1e-12`, merges oriented terms with hard-coded `1e-11`, and then applies the Gauss projection.

Its dependency `bqg_reduced_transfer_hamiltonian.py` (blob `5ec3a8b96189da2d2e2ddd095c09cf1a6371edc7`) contains further cutoffs:
- sparse state addition default `1e-13`;
- local coefficient acceptance at `1e-12`;
- local output filtering at `1e-12`;
- Hamiltonian branch accumulation at `1e-11`.

Therefore setting `tol=0` only at `master_map` is insufficient: some lower-level cutoffs are fixed internally and do not inherit the top-level argument. The producer requires a coordinated no-prune mode across the full call graph before column chunking or a tolerance-free rank experiment can be trusted.

### 1.2 Master-map numerical cutoffs

File: `research/depth6_recovered_source/bqg_depth6_generic_master_row.py`
Blob: `62fa5d3926f704bcec302e746b19859f6a9b4508`

The generic master-row implementation uses a default `tol=1e-11` in:
- `reduced_s3_node`;
- `reduced_s2_node`;
- `master_map`.

It discards intermediate blocks when `np.linalg.norm(Y) <= tol`, discards input segments when `np.linalg.norm(seg) < tol`, and removes final output blocks when `np.linalg.norm(M) <= tol`.

These are state-dependent magnitude cutoffs, not exact representation-theoretic selection rules. Consequently, a rank inferred from this floating-point producer is not an exact rank certificate. A future certificate must provide at least:
1. complete shard coverage and stable source hashes;
2. reproducible global row/column dimensions and block ordering;
3. rank stability under a documented tolerance sweep and higher precision;
4. a lower bound on retained singular values, preferably interval-certified or exact modular/rational witnesses;
5. explicit treatment of near-null singular directions and numerical error bounds.

### 1.3 Existing frontier

File: `depth6_frontier.json`
Blob: `0ad5894de68f13631624f786c2287d6e652901b9`

The frontier still marks:
- `[5]`: vacuum inclusion established for the recovered generic master-row formula, but nonvacuum rank unverified;
- `[4,1]`: no validated rank witness;
- `[3,2]`, `[3,1,1]`, `[2,2,1]`: structural block census is not an injectivity certificate;
- `[2,1,1,1]` and `[1^5]`: finite numerical evidence only, not a uniform proof.

The multiplicity identity is a consistency check, not a rank result:
`sum(d_lambda * m_lambda) = 3,111,637`.

## 2. C2: production accumulation/cutoff inventory

The following are live source-level cutoffs in the W5/Peter–Weyl path. They must not be conflated with exact SU(2) selection rules or fixed index restrictions.

| File | Blob (prefix) | Observed sites | Interpretation |
|---|---|---|---|
| `k5_peter_weyl_safe_hda_column.py` | `e7e1e070` | `add_dict` 1e-13; projection coefficients 1e-12; output filters 1e-11/1e-10; `prune_state` 1e-8; composition and H/V paths call pruning | High-priority source of intermediate nonlinear pruning |
| `peter_weyl_lorentzian_K_block_gate.py` | `0b4d097e` | projection coefficient 1e-12; volume output 1e-10; HE/K output 1e-9 | Local K construction inherits thresholding from base PW module |
| `peter_weyl_covariant_volume_leg_gate.py` | `7d1fb7fe` | coefficient cutoffs 1e-13; projected-state cutoff 1e-11; final inverse/forward cutoff 1e-10 | Covariant projection and map assembly |
| `peter_weyl_covariant_K_leg_gate.py` | `60b9b3ec` | `add_state` 1e-11; projection coefficients 1e-13; output filters 1e-11 and 1e-10 | Charged K-leg composition and HE accumulation |
| `peter_weyl_covariant_composition_gate.py` | `1ef229b1` | `add` 1e-11; complete-to-covariant output filter 1e-11 | Covariant composition |
| `peter_weyl_lorentzian_ordered_triple_gate.py` | `0f910620` | `add` 1e-11 | Ordered triple aggregation |
| `bqg_v10_w5_lorentzian_adjoint_closure.py` | `fa02e216` | `TOL=0` in `add_sparse`, but final `amp_tol=1e-12` filtering remains | No pruning in the helper itself; final state filter still exists |

### 2.1 Exact selection versus numerical pruning

The rule removing source (J=2) after the middle (C(K)) leg is a representation-theoretic selection rule only if the final rank-((0+1)) leg is fully evaluated and the intermediate state is represented in complete (J) sectors. That fixed sector restriction is not the same as deleting amplitudes because their current floating-point magnitude is below a tolerance.

Coefficient thresholds in basis reconstruction are also not automatically safe: they are numerical approximations to a projection, and must be audited separately from sparse-state accumulation. Replacing only the last output filter cannot repair earlier deletions.

### 2.2 Neumaier candidate scope

`scripts/bqg_neumaier_sparse_accumulator.py` is a primitive-level candidate. It is not yet integrated into the functions listed above. Its tests establish useful finite arithmetic behavior, but do not certify the composed W5 operator. A safe integration requires each accumulation boundary to retain the compensation state until the complete mathematical sum for that boundary has been formed.

## 3. CI result and correction

The PR head's `numerical-continuation` run 647 completed successfully.

The PR head's `core-regression` run 1430 reached:
- Neumaier primitive tests: **4 passed**;
- Python compileall: **passed**;
- `scripts/audit_core_scope.py`: **failed** because existing historical/audit documents and an evidence script contain retired terms in explanatory contexts.

The failure was not an accumulator test failure. The scope scanner has now been amended to permit those terms only in the specifically named audit/evidence files; all other files remain scanned, and retired-module import checks remain active. A fresh workflow run is required to validate that correction.

## 4. Required next gates

1. Run the updated core-regression workflow and confirm the scope scanner no longer blocks unrelated checks.
2. Build a no-prune/compensated variant of the actual PW dependency graph, preserving exact representation selection rules.
3. Compare baseline, compensated, and strict-reference outputs on the captured scalar-covariance failure seed.
4. Run basis-vector and superposition tests through the full W5 map, not only the accumulator primitive.
5. Complete the 30-cell paired V09/V10 audit before making any C2 closure statement.
6. Recover all C1 shards and produce per-irrep, per-block rank evidence with numerical error bounds.

## 5. Claims discipline

The current result is a source audit plus finite primitive arithmetic tests. It is not a proof of linearity, self-adjointness, domain control, W5 scalar covariance, or C1 injectivity. C1–C4 and all physicalization gates remain open.
