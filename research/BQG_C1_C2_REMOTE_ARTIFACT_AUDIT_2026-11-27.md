# Remote artifact audit — BQG C1/C2 — 2026-11-27

This note records additional GitHub Actions evidence recovered after the initial continuation audit. It does not change the formal gate statuses.

## C1 — depth-6 master-rank pipeline is not closed

The repository contains a real per-block computation path:

- `research/depth6_recovered_source/github_mixed_master_shard.py` builds the depth-6 shell, computes S5 orbit multiplicities, calls the recovered `master_map`, and writes per-orbit block matrices as pickle files.
- `.github/workflows/bqg-depth6-mixed-stage-a.yml` expects 112 shards for each of `[3,2]` and `[3,1,1]`, then a fail-closed aggregator requiring all 2,755/130,903 and 2,719/153,455 blocks/columns respectively.
- `.github/workflows/bqg-depth6-mixed-stage-b-final.yml` expects 160 shards for `[2,2,1]` and then verifies the combined frontier.

However, the latest Stage-A run found in the inspected run history was run **#9**, ID `42821076917`, on 2026-11-10. Its paginated job list has 112 `[3,2]` shard jobs:

- 16 succeeded;
- 93 failed;
- 3 were cancelled.

The successful jobs did not produce a complete shard set. Several failed/cancelled jobs show runner shutdown signals and exit code 143 after long-running shard computation. The associated Stage-A keymeta salvage run **#8**, ID `43286958315`, failed closed because the source Stage-A run was not a successful same-repository main-branch run. The corresponding Stage-B run **#10**, ID `43286958231`, was skipped. The old raw artifacts had short retention and are no longer available.

The earlier Stage-A run ID `36899125190` also failed and has no artifacts. Therefore, the 112-shard source path exists, but no complete distributed rank certificate was recovered. The checked-in `DEPTH6_5_ORBIT_BLOCK_CENSUS_20261126.json` explicitly has `rank_claim: null` and states that certified per-block ranks remain open.

There is a second qualification: the recovered shard producer is a floating-point, tolerance-filtered map. `github_mixed_master_shard.py` calls `master_map(..., tol=1e-11)`; the generic row implementation also prunes intermediate reduced S3/S2 coordinates and final blocks by norm thresholds (with a separate `1e-12` coordinate cutoff). Thus even a future complete shard aggregate would first be a finite numerical rank witness, not an exact-arithmetic proof. A robust certificate must include tolerance sweeps, residual bounds, and a lower singular-value margin that remains separated from the pruning/error budget.

**C1 remains OPEN.** Do not infer a rank result from the structural census, partial shard completions, or the separate 32D logical calculation below.

## C2 — actual W5 tolerance ablation

Recovered successful GitHub Actions run:

- Workflow: `W5 scalar tolerance ablation`, run #2, ID `42820210450`.
- Artifact: `w5-scalar-tolerance-ablation-k0-node0-j10`, SHA-256 `6237661d422f02698bce72d9aba6dc46fc5ea1571011f062499992e1eb07d111`.
- Input: node 0, source K=0, Jmax2=10.
- One branch-generation pass; 768 generated sequences and 98,304 generated branches.

On this specific finite input:

| Variant | Support | Norm | Allowed J=0 weight | Wrong-sector amplitude norm |
|---|---:|---:|---:|---:|
| Production projection/accumulation tolerance 1e-11 | 122 | 2.5770296247418205 | 6.6410816867969675 | 0 |
| Zero projection and accumulation tolerance | 1020 | 2.577029624741821 | 6.641081686796968 | 1.9580e-15 |

The full-map relative L2 difference against the `math.fsum` reference is approximately `2e-15` for the tested variants. Reversing the oriented-spec order in the zero-tolerance grouped run did not materially change that one result.

**Important interpretation:** this ablation shows a large support-set difference (122 vs 1,020) but negligible amplitude-level L2 difference for this one captured seed. It neither proves that threshold pruning is harmless for all W5 inputs nor certifies linearity of the complete operator. The artifact also records a remaining internal basis coefficient cutoff of `1e-13`; the zero-tolerance variant is therefore not a globally unpruned symbolic evaluation.

### Source-level threshold inventory

The following production paths still have amplitude-dependent pruning:

- `scripts/peter_weyl_covariant_K_leg_gate.py`: `add_state` deletes a key when a partial sum is at or below `1e-11` (lines 53–59); projection paths also use `1e-13`, `tol`, and final amplitude filters (lines 89–134, 181–185, 236).
- `scripts/peter_weyl_covariant_volume_leg_gate.py`: charged/Gauss projection uses coefficient and output amplitude filters (lines 98–151, 175).
- `scripts/peter_weyl_covariant_composition_gate.py`: thresholded `add` with deletion of existing keys (lines 45–54), plus final projection filter.
- `scripts/peter_weyl_lorentzian_ordered_triple_gate.py`: same intermediate thresholded add/delete pattern (lines 48–57).
- `scripts/bqg_v10_w5_lorentzian_adjoint_closure.py`: `add_sparse` defaults to `TOL=0`, but the final forward/adjoint columns are filtered by user-provided `amp_tol` (default `1e-12`, lines 238, 295–296).

Do not remove exact structural selection rules. In particular, V10's `prune_scalar_channel` retains only J=0 or J=1 because the final rank-(0+1) leg cannot couple other J sectors to source J=0; that is a fixed representation-theoretic rule, not a magnitude cutoff.

**C2 remains OPEN.** The Neumaier accumulator remains a candidate only. Integration must cover every accumulation and projection site, followed by basis/superposition/scaling tests and the full 30-cell paired W5 V09/V10 audit.

## Separate finite result — not a C1/C2 closure

A separate complete 32-column finite Euclidean higher-shell artifact was recovered from core-regression run #1419 (ID `45823637097`). I independently rebuilt its 32×32 Gram matrices from the saved sparse column states:

- 32/32 columns present; first-state supports 64–120 and second-state supports 1,253–3,000.
- (K) rank 32; eigenvalue range 4.3060759870–13.3527813527.
- (M=P H^4 P-K^2) eigenvalue range 47.9777767497–186.9023442232.
- (Lambda) eigenvalue range 10.6357598783–15.0599276660.
- Rebuilt (K) relative difference from the saved matrix: (1.92\times10^{-17}).
- Rebuilt (Lambda) relative difference: (3.33\times10^{-15}).

This verifies the saved finite 32D Euclidean artifact internally. It is not the depth-6 graph-changing master map, not the Lorentzian W5 operator, and not physicalization evidence.

## Current disposition

- C1: `OPEN` — complete per-block rank certificates not recovered.
- C2: `OPEN` — source-level nonlinear pruning risk; one actual finite W5 ablation is amplitude-stable for its single seed; production integration and paired W5 audit missing.
- C3/C4: unchanged, `OPEN`.
