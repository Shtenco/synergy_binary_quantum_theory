# Peter-Weyl physical projector: no-prune integration audit (2026-11-11)

## Scope

This change integrates the previously separate finite Peter-Weyl master-projector pilot into the consolidated research branch. It does not change `main`, the regulated K5 operator definition, or the physical interpretation of the model.

## Code changes

- Added `scripts/peter_weyl_master_projector_pilot.py`, based on the pilot from `research/post-hda-physicalization-v2`.
- Added explicit `--mode legacy|no_prune`; the selected mode is passed from `constraint_columns` to both `PW.apply_H_cached_state` and `PW.prune_state`.
- In `no_prune`, the producer's no-prune path is used and the outer `prune_state` call is a no-op.
- `master_from_gram` now returns the raw assembled matrix. Its Frobenius and relative Hermiticity defects are recorded before any Hermitian projection.
- The Hermitian part is formed only after diagnostics, for the eigensolver. A defect above `--hermitian-rtol` causes the internal acceptance check to fail.
- The zero sector uses `abs(lambda) <= zero_threshold`; positive semidefiniteness is tested separately against `psd_tolerance`. Negative eigenvalues outside the zero window are reported separately and are not counted as positive gaps.
- Tolerances and raw matrix inputs are checked for finiteness.
- Added 12 focused tests in `scripts/test_peter_weyl_master_projector_integrity.py`.
- Added both scripts to push and pull-request path filters and added the test command to `.github/workflows/bqg-sparse-pruning.yml`.

## Local algebraic control

An independent synthetic two-node/two-column calculation compared the sparse cross-Gram contraction

`M_ab = sum_vw G_vw <H_v a | H_w b>`

with a separately assembled dense matrix expression

`M = sum_vw G_vw H_v^dagger H_w`.

For the deterministic complex-valued fixture, the Frobenius difference was exactly 0 in float64, the raw Hermiticity defect was 0, and the Hermitian spectrum was approximately `[11.58082101, 14.31292899]`.

This is a synthetic algebraic control of the contraction/index convention, not a run of the production K5 physical projector.

## Reproduction commands

Legacy comparison:

`python scripts/peter_weyl_master_projector_pilot.py --domain full --jmax2 5 --mode legacy`

Threshold-free float64 finite-carrier calculation:

`python scripts/peter_weyl_master_projector_pilot.py --domain full --jmax2 5 --mode no_prune`

The full 32-state no-prune calculation has **not** been executed in this audit. Its output must not be described as PASS or as a physical no-go until it has completed and its artifacts have been independently inspected. No-prune removes amplitude thresholds; it does not provide exact arithmetic or eliminate floating-point roundoff.

## Current scientific boundary

Even a successful no-prune finite-carrier result is a regulated finite-dimensional Galerkin diagnostic. It is not a continuum rigging map, a continuum no-go theorem, or a complete quantum HDA closure proof. The next quantitative comparison is to assemble legacy and no-prune constraint columns on the same declared carrier and compute

`||C_no_prune - C_legacy||_2`

alongside `sigma_min(C_legacy)`, with all matrix dimensions, support hashes, spin cutoff, and numerical tolerances recorded.

## Validation status

- The 12 new integrity tests are committed but have not yet been observed executing in CI.
- The workflow runs created for this branch were queued at the last poll; no PASS is claimed.
- `main` remains unchanged.
