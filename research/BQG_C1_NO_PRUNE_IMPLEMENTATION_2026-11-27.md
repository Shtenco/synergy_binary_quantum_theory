# BQG C1: isolated no-magnitude-pruning path (2026-11-27)

## Purpose and status

This is a derived research execution path for investigating the effect of intermediate amplitude deletion in the depth-6 mixed master map. It is not a proof of injectivity, a rank certificate, a production run, or a closure of C1.

The byte-for-byte recovered source snapshot under `research/depth6_recovered_source/` is preserved. The derived implementation lives separately under `research/depth6_no_prune/`.

## Active call path

```text
github_mixed_master_shard.py
  -> bqg_depth6_generic_jucys_selector.py
       -> expand(..., tol=0)
  -> bqg_depth6_generic_master_row.py
       -> master_map(..., tol=0)
       -> reduced_s3_node / reduced_s2_node / raw_to_s2_coords (tol=0)
  -> bqg_block_apply_invariants.py
       -> H_node_apply_ultrafast(..., tol=0)
       -> radial / triangle collapsed local maps
  -> bqg_reduced_transfer_hamiltonian_no_prune.py
       -> local_hit_outputs(..., tol=0)
  -> bqg_covariant_volume_recoupling.py
       -> coupled_tensor: exact-zero checks only
```

The shard runner explicitly calls `master_map(..., tol=0.0)` and records `operator_mode=FLOAT64_NO_MAGNITUDE_PRUNING`.

## Removed intermediate magnitude filters on the active path

- Initial reduced rows: the recovered `TOL=1e-12` is passed as zero.
- Sparse vector accumulation and scaled merges: `tol=0`; only exact-zero results are deleted.
- Gauss projection: entrywise and norm filters receive zero.
- Local one-hit CG overlaps: `local_hit_outputs` receives zero.
- Volume-local coefficients and radial local-matrix entries: filters receive zero.
- Two-hit and shared-color local outputs: filters receive zero.
- Jucys expansion: the recovered `1e-13` segment cutoff is parameterized and defaults to zero in the derived selector.
- Coupled-tensor reconstruction: the recovered `1e-15` CG-factor cutoffs are replaced by exact-zero checks in the derived recoupling module.

The recovered `bqg_safe_boundary_master_gate.py` is not modified. The active path uses its CG matrix, epsilon, axis-application, edge-weight, spin-layout and oriented-spec helpers; the legacy thresholded factorized-map functions in that module are not called by `H_node_apply_ultrafast`.

## Numerical limits that remain

1. Arithmetic is still complex float64 with ordinary NumPy accumulation. This path disables magnitude pruning; it does not provide exact arithmetic or compensated accumulation. The Neumaier sparse accumulator remains a separate candidate and is not integrated here.
2. The selector still uses a machine-epsilon singular-value rank threshold and a residual acceptance tolerance (derived default (10^{-10})). These are basis/rank-generation controls, not sparse amplitude filters, and do not constitute a rank proof.
3. The S3 basis uses the existing spectral projector range selection (`eigenvalue > 0.5`) with a residual/orthogonality check. Its finite-dimensional basis construction is not an injectivity certificate.
4. The radial scalar-consistency guard at (10^{-10}) raises on a failed identity; it does not delete amplitudes.
5. The supplied selector/shard architecture currently covers four keys: `32`, `311`, `221`, and `2111`. It does not cover all seven S5 sectors.

## Regression coverage

`tests/test_bqg_depth6_no_prune_mode.py` checks:

- preservation of sub-threshold initial rows;
- cancellation terms surviving sparse accumulation when pruning is disabled;
- sub-threshold entries surviving Gauss projection;
- zero-tolerance propagation through the node wrapper;
- isolation of derived operator, transfer, recoupling and selector modules;
- explicit zero-tolerance routing through the shard worker;
- retention of sub-threshold Clebsch-Gordan reconstruction coefficients;
- retention of sub-threshold Jucys expansion coordinates.

These are primitive and wiring tests. A complete depth-6 shard run, all-sector block rank analysis, and independent rank verification remain outstanding.

## Status

- No-prune primitive/wiring tests: CI-gated.
- Full no-prune depth-6 production dependency graph: not executed.
- C1 vacuum-only kernel: **OPEN**.
- C2 Lorentzian/W5 operator closure: **OPEN**.
- `core_theory_closed_declared=false`.
- `experimentally_confirmed=false`.
