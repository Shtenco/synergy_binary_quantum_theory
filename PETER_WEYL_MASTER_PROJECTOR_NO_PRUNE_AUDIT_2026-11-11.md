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
- Added 15 focused tests in `scripts/test_peter_weyl_master_projector_integrity.py`.
- Added both scripts to push and pull-request path filters and added the test command to `.github/workflows/bqg-sparse-pruning.yml`.

## Follow-up integrity correction (2026-11-11 09:27 UTC)

A second-pass audit found a hidden coefficient cutoff in the zero-sector
annihilation verifier: coefficients with magnitude below `1e-14` were
skipped while reconstructing `H_v|psi>`. This could under-report the
constraint residual for a zero eigenvector with small but nonzero components.
The verifier has been factored into `constraint_annihilation_norm` and
now skips only coefficients that are exactly `0j`. Contributions are
accumulated with `math.fsum` separately for real and imaginary parts
before the Euclidean norm is formed.

Two regression tests exercise a `1e-15` coefficient multiplying a
`1e15` column amplitude (the residual must remain 1, not disappear
behind a threshold). This correction removes the last explicit amplitude
threshold in the zero-sector residual reconstruction; it does not establish
exact arithmetic.

## Producer-path audit after the verifier correction

The no-prune path was traced through the production K5 column producer,
including `apply_primitive_oriented`, `apply_T_oriented_basis`,
`apply_T_cached_state`, `apply_H_cached_state`, and
`prune_state`. In no-prune mode these layers route through
`coherent_aggregate` / `aggregate_mappings`, which use
deterministic real/imaginary `math.fsum` accumulation and do not apply
magnitude thresholds. The local channel/output guards are exact comparisons
against `0j`, not epsilon cutoffs. The declared finite spin regulator
`Jmax2` remains active by design.

The audit therefore distinguishes three things: (1) legacy amplitude pruning,
(2) no-prune coherent float64 aggregation with exact-zero elimination, and
(3) the finite spin cutoff. “No-prune” means no amplitude-threshold pruning;
it does not mean infinite spin, exact symbolic arithmetic, or absence of
floating-point cancellation/roundoff.

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

## Legacy/no-prune constraint-map perturbation harness

Added `scripts/peter_weyl_constraint_stability.py`. It constructs the
direct-sum map

`C psi = (H_0 psi, ..., H_4 psi)`

for the same finite carrier and spin cutoff in both pruning modes. It then
forms only the small input-space Gram matrices (C^\dagger C), using sparse
output dictionaries, and evaluates:

- `sigma_min(C_legacy)` (including numerical nullity);
- `sigma_min(C_no_prune)`;
- `||C_no_prune-C_legacy||_2`;
- support counts and deterministic SHA-256 digests of the stacked columns;
- the sufficient Weyl rank-stability inequality.

The script avoids materializing a potentially enormous dense output-space
matrix. Its small Gram matrices are assembled independently in both
conjugate directions, so the raw Hermiticity defect is not imposed to zero by
construction. Sparse inner products now use compensated real/imaginary
`math.fsum` sums, and the legacy operator cache is cleared and garbage-collected
before constructing no-prune columns to reduce peak memory.

Important numerical limitation: singular values are estimated from
`eigvalsh(C^dagger C)`, which squares the condition number. The script
therefore uses a roundoff-scaled PSD tolerance (`100 * machine epsilon`),
reports whether the declared Gram-eigenvalue rank threshold clears that
resolution floor, and refuses to certify the Weyl sufficient condition when
the legacy rank threshold does not resolve roundoff. A near-null result still
needs an independent high-precision or stable QR/SVD cross-check.

Added five synthetic dense-oracle tests in
`scripts/test_peter_weyl_constraint_stability.py`, including direct
comparison of the sparse perturbation norm with dense SVD. Both new files
are included in push and pull-request workflow path filters, and the tests
are part of the CI job. A local independent synthetic control of the new comparison algebra gave
exact agreement between the sparse Gram construction and dense matrix
multiplication (Frobenius difference 0); for the two-column perturbation
fixture, both sparse-Gram and dense-SVD operator norms were exactly 0.2.
This validates only the synthetic contraction/perturbation implementation,
not the production K5 column calculation.

The stability suite contains five tests. The combined integrity/stability
test count is therefore 20. CI is wired for both push and pull-request paths;
the GitHub Actions runs remain queued, so these tests are not yet claimed as
CI-passed.

The production full-32 comparison has not yet run.

Reproduction:

`python scripts/peter_weyl_constraint_stability.py --domain low --jmax2 5`

`python scripts/peter_weyl_constraint_stability.py --domain full --jmax2 5`

The low-domain run is a six-state smoke calculation; only the full-domain run
addresses the complete 32-state carrier.

## Current scientific boundary

Even a successful no-prune finite-carrier result is a regulated finite-dimensional Galerkin diagnostic. It is not a continuum rigging map, a continuum no-go theorem, or a complete quantum HDA closure proof. The next quantitative comparison is to assemble legacy and no-prune constraint columns on the same declared carrier and compute

`||C_no_prune - C_legacy||_2`

alongside `sigma_min(C_legacy)`, with all matrix dimensions, support hashes, spin cutoff, and numerical tolerances recorded.

## Validation status

- The 12 new integrity tests are committed but have not yet been observed executing in CI.
- The workflow runs created for this branch were queued at the last poll; no PASS is claimed.
- `main` remains unchanged.


## Follow-up numerical certificate hardening (2026-11-11)

The research branch now contains an isolated streaming singular-value implementation
in `scripts/peter_weyl_streaming_svd.py`. It accepts row blocks of the original
constraint map (C), applies reduced Householder QR to each block, merges compact
(R) factors in a deterministic binary tree, and performs a direct SVD of the final
small (R). It does not form (C^\dagger C).

The stability comparison now reconstructs deterministic row blocks from the sparse
columns and computes direct singular values for the legacy map, no-prune map, and
their represented-float difference. Gram spectra remain available as diagnostics,
but their rank threshold is explicitly a relative **Gram eigenvalue** threshold.
The direct SVD has a separate relative **singular-value** threshold. This avoids
silently treating those two tolerances as interchangeable.

The scale floor `max(1.0, ...)` was removed from the Gram spectral classifier.
The master-projector pilot now derives zero and PSD tolerances from the actual
spectral scale, with a zero-scale branch. Its raw Hermiticity diagnostic uses a
scaled sum-of-squares norm rather than a unit/tiny denominator floor.

Constraint annihilation now records both absolute and Frobenius-relative residuals:
[
r_v^{\rm rel}(\psi)=\frac{\|H_v\psi\|_2}
{\|H_v\|_F\,\|\psi\|_2}.
]
The relative residual is used for acceptance; absolute residuals remain in the
certificate for auditability. Norm evaluation uses scaled sum-of-squares to avoid
naive square overflow/underflow.

The master-projector pilot also snapshots and validates the producer's raw volume
Hermiticity registry. An empty registry, non-finite entry, negative norm/defect, or
defect above the declared tolerance fails closed. Non-finite diagnostics are
serialized as null rather than non-standard JSON NaN/Infinity tokens.

### Numerical controls performed locally

- 11 TSQR unit tests passed against dense SVD oracles, including complex inputs,
  rank-deficient maps, full right-nullspace output, row permutation, block-size
  changes, controlled singular values down to (10^{-14}), scale changes from
  (10^{-150}) to (10^{150}), and explicit dynamic-range-underflow detection.
- A separate sparse-row assembly probe passed for duplicate output keys and the
  empty-support zero operator.
- Separate numerical probes passed for scale-invariant Gram thresholds,
  scale-invariant Hermiticity ratios, relative annihilation residuals, large
  complex norms, and fail-closed volume-defect handling.

These are local synthetic/helper tests, not a full Peter-Weyl production run and
not GitHub Actions results. The new direct-SVD comparison has **not** yet been
run on the full 32-dimensional no-prune carrier. No upstream operator-construction
error bound has been established, so the result remains numerical for the
represented float64 maps and cannot certify exact rank of the symbolic operator.
