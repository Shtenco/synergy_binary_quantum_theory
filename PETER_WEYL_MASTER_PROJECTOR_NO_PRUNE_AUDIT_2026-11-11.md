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

## Follow-up re-audit: direct TSQR/SVD and current test inventory (2026-11-11 10:09 UTC)

The current branch version supersedes earlier test-count statements in this
report. Counting test methods in the checked-in files at the current branch
head gives:

- `test_peter_weyl_master_projector_integrity.py`: 22 tests;
- `test_peter_weyl_constraint_stability.py`: 10 tests;
- `test_peter_weyl_streaming_svd.py`: 11 tests;
- total across these three focused suites: **43 tests**.

The workflow runs all three suites. This count is a static inventory of test
methods, not a statement that CI has passed.

### Independent implementation review

The streaming solver uses binary-tree Householder TSQR and a final SVD of the
compact R factor. The implementation avoids forming the Gram matrix for the
direct singular-value path. The row-block adapter deterministically merges
the sparse output supports and inserts an explicit zero row for an empty
operator, allowing the zero-map nullspace to be represented.

The TSQR implementation's current `scale_separation_loss` flag detects
components that underflow to exact zero during common rescaling/normalization.
It does not estimate general QR backward error or flag every sub-ulp loss of
relative information. Consequently the solver remains a floating-point
numerical diagnostic; its own output explicitly does not certify exact rank
or bound upstream operator-construction error. A future hardening pass should
add a conservative near-machine-precision scale-separation warning and a
regression fixture where a nonzero block is below machine epsilon relative to
the dominant block.

### Current repository and CI snapshot

At this audit point the research branch head is
`2ba2c8e87b85de99012d1905835d7a9d7f8207d1`. Comparison with `main`
reports 259 commits ahead and 2 behind, with merge base
`23669209c83fe0aa2741612e5b2b119562ee6d7d`. The `main` head remains
`254f19fdab41d32b920bf454f987aa7c0fdba34d`; no merge or main-branch write
was performed.

The latest observed workflow run is
`42999838772` for commit
`2ba2c8e87b85de99012d1905835d7a9d7f8207d1`. At the time of this check it
was still `queued`, and its `sparse-pruning-tests` job was also `queued`.
No CI PASS is claimed. Several earlier pushes also created queued runs; the
latest run is the relevant branch-head check.

The full 32-state production legacy/no-prune comparison has still not been
executed. No physical-kernel conclusion is upgraded by the presence of the
new TSQR implementation or by synthetic dense-oracle tests.

## Follow-up: subnormal Gram threshold regression (2026-11-11 10:10 UTC)

A review of the normalized Gram-rank classifier found that the absolute
quantity `rtol * spectral_scale` can underflow to zero for subnormal Gram
matrices even when the normalized eigenvalues remain well resolved. The
classifier now compares normalized eigenvalues directly with the relative
threshold and reports both normalized and absolute threshold quantities.

Added a regression test using a Gram matrix with scale (10^{-320}). It
requires the absolute threshold to underflow to zero while the normalized
spectrum still classifies rank one. The focused suite inventory is now:

- master-projector integrity: 22 tests;
- constraint stability: 11 tests;
- streaming TSQR/SVD: 11 tests;
- **44 tests across these three focused suites**.

This test is committed but has not yet been observed passing in GitHub Actions.
The test-only change is covered by the workflow's path filters. No production
32-state run has been performed.

## Follow-up: conservative TSQR scale-resolution gate and expanded producer-ledger audit (2026-11-11)

A further source-level pass traced the no-prune path through the production
K5 implementation. The amplitude-threshold branches in
`apply_primitive_oriented`, `apply_T_oriented_basis`,
`apply_T_cached_state`, `apply_H_cached_state`, and `prune_state`
are bypassed in `mode="no_prune"`; aggregation is routed through the
compensated sparse aggregation helpers. The remaining explicit cutoff in
this path is the declared finite spin regulator `Jmax2`, not an amplitude
threshold. The legacy-only threshold expressions remain in the legacy
branches. The five-tetrahedron input vector also uses exact-nonzero component
selection in no-prune mode.

The volume-operator audit was strengthened independently: each projector run
clears the volume cache and raw-defect registry, then requires the number of
unique registry records to match the cache-miss count. It also validates
spin-triple shape/type and that every doubled spin is within the declared
`Jmax2`. Malformed records, malformed keys, missing coverage, or out-of-cutoff
entries fail closed. This closes a weakness where a nonempty but incomplete
registry could previously pass the Hermiticity gate.

### TSQR scale-resolution warning

The streaming TSQR implementation now conservatively flags two cases even
when no represented component underflows to exact zero:

1. a nonzero normalized input component is below machine epsilon relative to
   the common input scale;
2. a previously accumulated compact factor is rescaled by a factor below
   machine epsilon when a later block raises the global scale.

The warning deliberately over-approximates possible loss: an orthogonal
sub-ulp component can survive in a special matrix geometry, and the regression
fixture demonstrates exactly that the numerical rank may still be reported as
two while the scale-resolution status is inconclusive. The warning is not a
proof that the singular value was destroyed; it blocks treating the result as
a robust certificate without a backward-error analysis.

Two new regression cases exercise (i) a sub-ulp component within one block
and (ii) a small block arriving before a much larger block. A local synthetic
reproduction of both cases retained the represented singular values
`[1, 1e-20]` and rank two, while correctly raising the conservative warning.
A separate 113-by-9 complex random dense-oracle reproduction agreed with dense
SVD to maximum absolute singular-value difference about (2.2\times10^{-14})
and did not raise the warning. These are local synthetic probes, not a full
suite execution or production K5 calculation.

### Current focused test inventory

At the latest source snapshot checked during this pass:

- master-projector integrity: 27 test methods;
- legacy/no-prune constraint stability: 12 test methods;
- streaming TSQR/SVD: 13 test methods;
- total across these three focused suites: **52 test methods**.

The workflow invokes all three suites. The count is static source inspection;
the corresponding GitHub Actions run is still queued and no PASS is claimed.

### Current branch and scientific status

The research branch head at the end of this pass is
`88038f1ff8fa89a11e15f1454c9a2615eea518d4`; `main` remains
`254f19fdab41d32b920bf454f987aa7c0fdba34d`. The branch remains isolated and
has not been merged. The latest observed Actions run for the test commit was
queued, with its job also queued.

The full 32-state production legacy/no-prune comparison has not been executed.
No physical-kernel, exact-rank, or continuum conclusion is established by
these additional numerical controls.

## Follow-up: volume symmetrization error and spectral-resolution ledger (2026-11-11 13:45 UTC)

The raw volume producer now records a quantitative bound for the Hermitian
projection actually used by the implementation. With

\[
Q_H=\frac{Q+Q^\dagger}{2},
\]

the identity

\[
Q-Q_H=\frac{Q-Q^\dagger}{2}
\]

gives

\[
\|Q-Q_H\|_F=\frac12\|Q-Q^\dagger\|_F,
\qquad
\|Q-Q_H\|_2\leq\frac12\|Q-Q^\dagger\|_F.
\]

The producer records the absolute Frobenius displacement, the same quantity
as a valid upper bound on the spectral-norm displacement, and its ratio to
\(\|Q\|_F\). The projector audit independently checks consistency among
the absolute and relative records and rejects malformed or inconsistent
certificates. This is a rigorous bound for the change from raw Q to its
Hermitian projection Q_H, not an error bound for the nonlinear volume
map V=sqrt(|Q_H|).

The Hermitian spectrum now also receives a scale-normalized numerical
resolution diagnostic: spectral scale, relative resolution threshold,
numerical nullity at that threshold, normalized gap to zero, minimum
resolved absolute eigenvalue, and the maximum scalar derivative
1/(2 sqrt(|lambda|)) over resolved eigenvalues. The classification is
performed in normalized units so the absolute threshold may underflow
without changing the relative classification. Tests include a subnormal
scale probe.

The derivative diagnostic is deliberately labelled as a scalar
eigenvalue-map sensitivity indicator. It is not an operator-Lipschitz bound,
does not control eigenvector rotation, and is not a certified eigensolver
backward-error estimate. No stability claim for the full matrix-valued
sqrt(|Q_H|) map is made from this number alone.

### Current focused test inventory (static source count)

- master-projector integrity: 28 test methods;
- legacy/no-prune constraint stability: 12 test methods;
- streaming TSQR/direct-SVD: 13 test methods;
- production K5 raw-volume operator: 12 test methods;
- total across these four focused suites: **65 test methods**.

The workflow path filters and steps include all four suites. These are
committed test methods, not an execution receipt. No local execution of the
current branch snapshot and no successful GitHub Actions run are claimed.

### Current branch state at this update

The latest code/test commit before this documentation update was
8bc3faaf07ccd22b0845f0ea64b9f85fd54b5301. The research branch was
293 commits ahead of main and two commits behind it. The latest commit
had no reported combined-status checks and no associated workflow runs in
the connector response. Earlier queued runs remain non-successful,
non-failing pending work. The 32-state production comparison remains
unexecuted. Main was not changed or merged.

## Follow-up: physical-kernel subspace geometry (2026-11-11 14:23 UTC)

The legacy/no-prune comparator now goes beyond rank stability. It extracts the
right-kernel basis from the full right-singular-vector matrix returned by the
streaming direct SVD, forms

\[
P_L=Z_LZ_L^\dagger,\qquad P_N=Z_NZ_N^\dagger,
\]

and records

\[
d_{\rm ker}=\|P_L-P_N\|_2.
\]

For equal nonzero nullities, this is the sine of the largest principal angle
between the two finite-dimensional numerical kernels. The report also checks
kernel-basis orthonormality and projector Hermiticity/idempotence. If nullities
differ, it explicitly returns `NULLITY_MISMATCH` and suppresses the
principal-angle field. If either direct SVD carries the conservative
TSQR scale-separation warning, the geometric comparison is marked
`INCONCLUSIVE_SCALE_SEPARATION`; a computed projector distance is retained
only as a diagnostic, not a certificate.

Four regression methods were added:
- analytically rotated one-dimensional kernels, with the projector distance
  compared against \(|\sin\theta|\);
- identical kernels, requiring zero projector distance;
- unequal nullities, requiring no principal-angle claim;
- scale-separation warning, requiring an inconclusive status.

An independent local NumPy reproduction of the analytic geometry gave
\(0.17213833761259534\) against the expected
\(\sin(0.173)=0.17213833761259542\); identical subspaces gave zero, and
the nullity-mismatch projector distance was one. This is a synthetic
mathematical sanity check, not execution of the checked-in test suite or a
production Peter-Weyl calculation.

### Current focused test inventory (static source count)

- master-projector integrity: 29 methods;
- legacy/no-prune constraint stability: 16 methods;
- streaming TSQR/direct-SVD: 13 methods;
- production K5 raw-volume operator: 14 methods;
- total: **72 test methods**.

All four suites are included in the workflow definition. These counts are
source inventory only; the latest Actions runs remain queued and no CI PASS
is claimed.

### Current repository and scientific boundary

At this update, the research branch head is
`0de647f4cfc80cb9b85e87361d8a33b748c1c46e`; comparison with main reports
297 commits ahead and two behind, with merge base
`23669209c83fe0aa2741612e5b2b119562ee6d7d`. Main remains
`254f19fdab41d32b920bf454f987aa7c0fdba34d` and was not changed or merged.

The latest observed Actions run is `43037834137` for the current branch
head; its workflow job is queued and has no conclusion. Runs for the two
kernel-geometry commits are also queued. No full 32-state production
legacy/no-prune comparison has been executed. No exact operator rank,
physical-kernel equivalence, continuum limit, or physical prediction is
established by these changes.

## Follow-up: independent SVD oracle and dimension-aware spectral ledger (2026-11-11 14:24 UTC)

The production K5 volume regression suite now contains a second, independent
construction of the absolute Hermitian operator. For each declared small-spin
case it forms Q_H, computes an SVD

\[
Q_H=U\Sigma W^\dagger,
\]

and constructs

\[
|Q_H|_{\rm SVD}=U\Sigma U^\dagger,
\qquad
V_{\rm SVD}=U\sqrt{\Sigma}U^\dagger.
\]

For Hermitian Q_H, this is the positive square root of Q_H^dagger Q_H and
therefore independently checks both V^2=|Q_H| and V itself against the
producer's Hermitian-eigendecomposition path. The test also exercises the
all-zero (0,0,0) spin triple and requires a finite, zero-valued raw-volume
ledger entry.

The master-projector audit now cross-checks the recorded scalar sensitivity
against 1/(2 sqrt(lambda_min_resolved)), checks the normalized spectral gap
against lambda_min_resolved / spectral_scale, verifies the reported threshold,
and constrains the numerical nullity by the exact tensor-product dimension
(s1+1)(s2+1)(s3+1). The zero-spectrum case must report the full matrix
dimension as numerical nullity. Boolean values are rejected where integer
spin indices or coverage limits are expected, preventing Python's bool-is-int
subclassing from admitting malformed metadata.

### Updated focused test inventory (static source count)

- master-projector integrity: 31 test methods;
- legacy/no-prune constraint stability: 16 test methods;
- streaming TSQR/direct-SVD: 13 test methods;
- production K5 raw-volume operator: 14 test methods;
- total across these four focused suites: **74 test methods**.

These counts were obtained by source inspection. The tests have been added to
the branch; they have not been executed in this environment. The current
commit has no reported combined-status entries and no associated workflow
runs in the connector response. No CI PASS is claimed.

### Branch state

The latest confirmed commit at this update is
370f36562a078583fe48261ed65f3feb2157fa2d. The research branch is 304
commits ahead of main and two commits behind it, with merge base
23669209c83fe0aa2741612e5b2b119562ee6d7d. Main remains unchanged and no
merge was performed. The full 32-state production legacy/no-prune run and
the low-domain pilot execution remain unperformed; source changes alone do
not establish a physical kernel or a continuum result.

## Follow-up: independent volume oracle and CI queue inspection (2026-11-11 14:25 UTC)

A second independent oracle is now committed in
scripts/test_k5_volume_raw_operator.py. It constructs |Q_H| and sqrt(|Q_H|)
from the left singular vectors and singular values of an SVD of Q_H, then
compares both the squared production volume and the production volume itself.
A zero-spin (0,0,0) case checks that the raw operator, volume, and ledger are
all exactly zero and that the spectral nullity equals the one-dimensional
carrier size.

The integrity gate now also rejects boolean spin indices and boolean
coverage/cutoff arguments. This is intentional: Python treats bool as a
subclass of int, but True is not a valid doubled spin or call count.

### Current focused test inventory (static source count)

- master-projector integrity: 31 test methods;
- legacy/no-prune constraint stability: 16 test methods;
- streaming TSQR/direct-SVD: 13 test methods;
- production K5 raw-volume operator: 14 test methods;
- total: **74 test methods**.

### CI state: queue saturation observed

A direct Actions API read for this branch returned 122 workflow runs total;
the latest 100 returned runs were all still in status queued. The latest
runs include:
- run 122, ID 43037988475, head 370f36562a078583fe48261ed65f3feb2157fa2d;
- run 121, ID 43037980024, head ed066c8d0b57d26db5369a3a220d0e6dff4c23e4.

Their jobs are also queued and have no steps yet. This confirms a queue
backlog rather than a test failure. To avoid worsening that backlog, no
manual workflow dispatch or rerun is initiated. No test execution or CI PASS
is claimed.

The current branch head before this documentation update is
370f36562a078583fe48261ed65f3feb2157fa2d; comparison reports 304 commits
ahead of main and two behind. Main is unchanged and unmerged. The low-domain
and 32-state production calculations remain pending actual execution.

## Follow-up: spectral-threshold sensitivity of the physical kernel (2026-11-11 14:36 UTC)

A review of the comparator's two rank tolerances found a semantic mismatch
that must remain explicit in every interpretation:

- `singular_rtol` is applied to singular values of C;
- `rank_rtol` in the Gram diagnostic is applied to eigenvalues of
  C^dagger C, so its equivalent singular-value threshold is
  `sqrt(rank_rtol)`.

Thus the defaults `singular_rtol=1e-10` and `rank_rtol=1e-10`
do not describe the same numerical rank decision: the Gram threshold
corresponds to a singular-value relative threshold of 1e-5. This is not
an algebraic contradiction; it is a tolerance-semantics difference.
Reports must not compare the two rank fields as though they were computed
with a common cutoff.

The stability comparator now includes a declared singular-value threshold
grid (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4). For each legacy and no-prune
map it records:
- numerical rank and nullity at each grid point;
- largest retained and largest discarded normalized singular values;
- distances of those boundary values from the threshold;
- sampled rank plateaus;
- legacy/no-prune kernel-projector distance at each identical threshold;
- explicit nullity mismatch and TSQR scale-resolution status.

This is a grid-based sensitivity diagnostic, not a proof for thresholds
between grid points and not an exact-rank certificate. It also does not
supply the missing upstream operator-construction error bound.

Three regression methods were added: known rank transitions at singular
scales 1, 1e-6, 1e-12; rejection of unsorted/non-finite threshold grids;
and principal-angle comparison across a shared threshold grid.

Independent local NumPy oracle checks returned ranks [3, 2, 1] for the
matrix with singular values [1, 1e-6, 1e-12] at thresholds
[1e-13, 1e-9, 1e-5]. A rotated one-dimensional kernel at angle 0.1 gave
projector distance 0.09983341664682824 versus sin(0.1) =
0.09983341664682815. These are independent synthetic checks, not execution
of the checked-in suite or a production Peter-Weyl calculation.

### Current focused test inventory (static source count)

- master-projector integrity: 31 methods;
- legacy/no-prune constraint stability: 19 methods;
- streaming TSQR/direct-SVD: 13 methods;
- production K5 raw-volume operator: 14 methods;
- total: **77 test methods**.

### Current repository and CI state

The research branch head is
`4854eb2d7f9c029da209259d2006cb4aea5264c0`; comparison with main reports
308 commits ahead and two behind, with merge base
`23669209c83fe0aa2741612e5b2b119562ee6d7d`. Main remains
`254f19fdab41d32b920bf454f987aa7c0fdba34d` and was not changed or merged.

The latest observed Actions run is `43040087118` for the threshold-test
commit; its job remains queued with no conclusion. The preceding implementation
run `43040079658` is also queued. No CI PASS is claimed. The queue backlog
continues; no manual reruns or workflow dispatches were started.

The low-domain and full 32-state production legacy/no-prune comparisons have
not been executed. No physical-kernel equivalence, exact rank, regulator
stability, or continuum theorem is established.


## Follow-up: singular-threshold and subspace-geometry hardening (2026-11-11 14:38 UTC)

The threshold-sensitivity helper computes numerical ranks from direct singular
values of the represented constraint map C, normalizes by its largest singular
value, and compares right-kernel projectors at an identical declared threshold
grid. For equal-dimensional kernels, the spectral norm of the projector
difference equals the sine of the largest principal angle. Grid plateaus remain
sampling diagnostics only; they do not certify exact rank or exclude a
transition between sampled thresholds.

The input-integrity gate now rejects non-positive or boolean dimensions,
boolean threshold entries, and full right-singular-vector matrices that are
not orthonormal within a dimension-scaled floating-point tolerance. This is
necessary because projector-distance and principal-angle claims rely on an
orthonormal singular-vector basis; finite entries and correct array shape alone
are insufficient.

The threshold-stability suite now contains 20 test methods. Static source
inventory across the four focused suites is therefore 78 methods:
31 master-projector integrity, 20 constraint stability, 13 streaming
TSQR/direct-SVD, and 14 raw-volume operator tests. These are source counts,
not execution receipts. The newly added method and validation have not been
executed in this environment.

This hardening is based on branch head b31517b0c94f297813562c930a152b2ede346e0e. The branch was observed as
309 commits ahead of main and two behind immediately before this patch.
Main remains untouched; no merge was performed. CI was still queued at the
last observed API snapshot, and no manual rerun or dispatch was initiated.


## Follow-up: fail-closed validation of kernel SVD certificates (2026-11-11 14:49 UTC)

A further audit found that the threshold-sweep path checked full right-singular-vector
unitarity, but the primary `kernel_projector_from_direct_svd` path did not.
That left a validation asymmetry: a malformed full-V matrix could reach the
kernel-projector comparison without first proving that it was a unitary SVD
basis. The projector's own idempotence and Hermiticity checks are not a
substitute for validating the full singular-vector frame.

The kernel-projector constructor now fails closed unless all of the following
hold:

- dimension is a positive integer and not a boolean;
- rank is an integer in range and not a boolean;
- the singular spectrum is a finite, nonnegative, descending vector of the
  declared dimension;
- the record contains a finite nonnegative `singular_value_rtol`;
- the declared rank exactly matches the thresholded singular spectrum;
- the full complex right-singular-vector matrix is finite and unitary within
  the same dimension-scaled floating-point tolerance used by the threshold
  sweep.

The constructor now records the full-frame unitarity defect in addition to
the kernel-basis and projector geometry defects. Three negative regression
tests were added for a non-unitary frame, a rank/spectrum inconsistency, and
boolean rank/dimension values.

Static source inventory now counts 81 methods across the four focused suites:
31 master-projector integrity, 23 constraint stability, 13 streaming
TSQR/direct-SVD, and 14 raw-volume operator tests. This is a source count only;
the tests have not been executed locally in this environment.

The code and test commits are:
- `97e14833cecd1015aebbeaf8a9fd00349a4aa73f` — validate the full SVD certificate;
- `f908871db55b460126d9daff87b5ec9aa108edff` — add negative regression tests.

The latest confirmed branch head is `f908871db55b460126d9daff87b5ec9aa108edff`.
It is 312 commits ahead of main and two behind; merge base remains
`23669209c83fe0aa2741612e5b2b119562ee6d7d`. Main remains
`254f19fdab41d32b920bf454f987aa7c0fdba34d`; no merge was performed.

GitHub Actions run `43042013274` (test commit) and run `43041993637`
(implementation commit) are both queued with no conclusion. The branch has
127 recorded workflow runs in the latest API response, confirming that the
queue backlog persists. No rerun or manual dispatch was initiated, and no CI
PASS is claimed. The low-domain and full 32-state production calculations
remain unexecuted; no physical-kernel equivalence, exact rank, regulator
stability, or continuum theorem is established.


## Follow-up: regulator and tolerance argument validation (2026-11-11 14:50 UTC)

The stability runner accepts a doubled-spin cutoff `jmax2` that is passed
into the elementary spin-raising/lowering branches. Since both supported
carrier domains start from links with doubled spin 1 (j=1/2), a cutoff below
1 does not contain the declared input carrier and should not silently run as
though it were a valid regulator. The runner now rejects non-integer,
boolean, or sub-1 `jmax2` values before constructing any columns.

The same early validation now rejects boolean rank tolerances (Python's
`bool` is an `int` subclass), nonnumeric/nonfinite/negative tolerance
values, and nonpositive or boolean streaming block sizes. One negative test
checks all six invalid argument classes and confirms the failures occur
before the expensive Peter-Weyl calculation is entered.

The constraint-stability suite now contains 24 test methods; the four-suite
static source inventory is **82 methods** (31 + 24 + 13 + 14). This remains
a static count, not an execution result.

Commits:
- `5ca5062104752a4e23a026f7fd55f4e8d44883a1` — validate regulator and numerical tolerances;
- `af67604a36c6220dc9356c8110ede1c32260b854` — add invalid-input regression test.

The research branch remains isolated from main. GitHub Actions continues to
queue push-triggered runs; no CI result is inferred from queued jobs. No
production low-domain or 32-state calculation is claimed.


## Follow-up: finite spin-cutoff sensitivity sweep (2026-11-11 14:52 UTC)

A first regulator-sensitivity driver is now available in
`scripts/peter_weyl_constraint_stability.py`:

```bash
python scripts/peter_weyl_constraint_stability.py --domain low --jmax2-sweep 1,3,5
```

For each requested doubled-spin cutoff it records the legacy and no-prune
numerical ranks/nullities, largest and smallest singular values, smallest
retained singular value, normalized singular spectrum, and sparse-column
digest. Between adjacent cutoffs it compares the no-prune right-kernel
projectors and the normalized singular-spectrum shapes. The driver clears
the T-column cache and runs garbage collection between cutoff points so
transient cached constraint columns do not accumulate across the sweep.

Two synthetic regression tests were added:
- adjacent-cutoff kernel comparison reproduces the analytic
  `sin(theta)` projector distance for a rotated one-dimensional kernel;
- unordered, duplicate, boolean, and sub-1 cutoff entries are rejected
  before the expensive calculation is invoked.

Important scope boundary: the driver does **not** retain the full sparse
constraint maps across cutoffs, and therefore does not calculate
`||C_J-C_(J+1)||_2`. Normalized spectrum differences are shape diagnostics,
not operator-norm convergence bounds. Even a stable kernel projector across
the sampled cutoffs would be evidence only for the finite sampled regulators,
not regulator independence or a continuum theorem. The next stronger step
must preserve a common row-key embedding and directly estimate the
inter-cutoff operator difference, with a separate error ledger for floating
point and truncation.

The constraint-stability suite now has 26 test methods; the four-suite
static source inventory is **84 methods** (31 + 26 + 13 + 14). These are
source counts only. The new cutoff sweep and its tests have not been run
against the production Peter-Weyl columns in this environment.


## Follow-up: scale-loss gates across kernel and regulator diagnostics (2026-11-11)

The streaming TSQR routine flags cases where common-scale normalization or
restoration may erase singular modes. The audit found that downstream
consumers could still build projectors or report threshold ranks from the
physical-scale singular-value array before honoring that flag. In an
underflow case this could turn a lost nonzero mode into an apparent exact
zero mode. The Weyl summary could also display a strict inequality computed
from the same lossy spectrum despite separately marking the run unstable.

The constraint-stability layer now fails closed before kernel-projector
construction and before threshold-grid rank classification whenever either
input map reports scale-separation loss. The Weyl numerical values and
full-rank flags become unavailable in that case. The finite-cutoff sweep
suppresses singular-spectrum summaries and adjacent-cutoff comparisons if
either cutoff loses scale separation. This prevents a numerical information
loss warning from being converted into a false zero-mode or convergence
statement.

One synthetic regression method was added. Static test-method inventory is
85 methods across the four focused suites (31 master-projector integrity,
27 constraint stability, 13 streaming TSQR/direct-SVD, 14 raw-volume
operator tests). This is source inspection only; no test execution is claimed.

This change is based on branch head 9dd4f9eb73922d4355bf3648575b8f6b4a4c4479. Main remains untouched and no
merge is performed. It does not establish regulator convergence or
physical-kernel equivalence.


## Follow-up: direct adjacent-regulator operator-difference norm (2026-11-11 14:56 UTC)

The initial spin-cutoff sweep compared kernel projectors and normalized
singular-spectrum shapes, but did not measure the operator change between
successive cutoffs. That is now strengthened.

For each adjacent pair, the runner temporarily retains the sparse constraint
columns and forms the difference in a common row-label embedding:
`Delta C = C_(J_next) - C_J`. The ambient row set is the union of the
canonical direct-sum labels `(node, Peter-Weyl output-state)`; a row absent
at one cutoff is zero-extended. Since the state key carries the spin and
intertwiner labels, identical labels refer to the same canonical output
basis vector. The difference is then passed to the direct streaming
TSQR/SVD, and its largest singular value is reported as the represented
float64 operator 2-norm. Both no-prune and legacy cutoff differences are
reported, together with the previous kernel-projector and normalized-spectrum
diagnostics.

The comparison helper validates the carrier dimension, complete input-column
index set, non-boolean finite nonnegative singular tolerance, and positive
block size. If TSQR reports scale-separation loss, the operator norm and
difference rank are suppressed rather than inferred from a possibly lossy
spectrum. The sweep retains only the previous cutoff's sparse maps while
building the next cutoff, avoiding accumulation of all cutoff maps in memory.

Synthetic tests now cover:
- a common row union with a newly appearing output state and an analytic
  operator-difference norm of 0.5;
- the adjacent-cutoff sweep wiring, including operator norm and kernel angle;
- suppression of the norm and rank when scale-separation loss is flagged;
- rejection of invalid tolerance, block-size, and incomplete-column inputs.

This is a meaningful upgrade from shape-only comparison, but it is not yet
a convergence theorem. The norm is exact only for the represented float64
matrices up to the numerical SVD error; it does not bound errors in
recoupling coefficients, finite-spin truncation of the underlying theory,
or the correctness of the canonical basis embedding itself. Those require
independent producer-level and representation-level error bounds. In
particular, a small adjacent difference over a few cutoffs is not enough to
prove a Cauchy sequence or a regulator-independent physical kernel.

The constraint-stability suite now contains 30 test methods. The static
inventory across the four focused suites is **88 methods**: 31 master-projector
integrity, 30 constraint stability, 13 streaming TSQR/direct-SVD, and 14
raw-volume operator tests. These counts are source inspection only; the
updated tests and production cutoff sweep have not been executed here.

Commits for this stage:
- `6d9e6daf68077e1eca1345dfbcd48e36052e74a9` — add adjacent-cutoff
  common-row operator norm;
- `184631ab25bac0b836b03fe102b7ec9404fd357a` — wire operator norm,
  kernel-projector and normalized-spectrum comparisons into the sweep;
- `555fbb7d0bee58af0f9fc2e60582c2dd6078a461` — test common-row embedding
  and operator norm;
- `2fdfab5d8a27d42ccda7dacbd781bd595735de4c` — validate comparison inputs;
- `a7956fcd508074b4380286fbd328955e67190c0a` — test invalid controls.

The branch head at this report update is `a7956fcd508074b4380286fbd328955e67190c0a`,
326 commits ahead and two behind main. Main remains unchanged and no merge
was performed. GitHub Actions remains queued; no CI PASS or local test
execution is claimed. The full production low-domain and 32-state runs
remain outstanding.


## Follow-up: common-row adjacent-cutoff norm and fail-closed rank metadata (2026-11-11 14:58 UTC)

The regulator sweep now computes the direct operator 2-norm of each adjacent
cutoff difference for both pruning modes. The maps are compared in the union
of their canonical direct-sum output labels, with absent rows zero-extended.
The largest singular value of the difference map is the reported
`operator_difference_norm_2`. A synthetic dense-oracle case with one newly
appearing row has the expected norm 0.5 and is covered by a regression test.

The sweep now reports, for every adjacent cutoff pair:
- no-prune and legacy `||C_(J_next)-C_J||_2`;
- difference-map numerical rank at the declared singular tolerance;
- scale-separation warning and a stable digest of the represented difference;
- no-prune kernel-projector distance;
- normalized singular-spectrum shape difference when resolvable.

The comparison helper rejects incomplete input-column maps, boolean/nonfinite
or negative tolerances, and invalid block sizes. A scale-separation warning
suppresses the operator norm and difference rank.

A further audit of the fail-closed path found that the nested direct-SVD
record still exposed the normalized SVD's integer rank even when the scale
warning made that rank inconclusive. The record now sets `rank: null` under
that warning and preserves the internal count only under the explicitly
untrusted field `untrusted_normalized_scale_rank`. A synthetic regression
test verifies this behavior.

The constraint-stability suite now has 31 test methods. The four-suite static
inventory is **89 methods** (31 master-projector integrity, 31 constraint
stability, 13 streaming TSQR/direct-SVD, 14 raw-volume operator). These are
source counts only. No checked-in tests or production cutoff sweep have been
executed in this environment.

Commits:
- `6d9e6daf68077e1eca1345dfbcd48e36052e74a9` — adjacent-cutoff difference
  operator SVD;
- `184631ab25bac0b836b03fe102b7ec9404fd357a` — sweep integration;
- `555fbb7d0bee58af0f9fc2e60582c2dd6078a461` — common-row norm regression;
- `2fdfab5d8a27d42ccda7dacbd781bd595735de4c` — helper input validation;
- `a7956fcd508074b4380286fbd328955e67190c0a` — validation tests;
- `022dca65f648396ef0ff1552183ab2b7f766aeea` — suppress inconclusive rank;
- `00b7b43eaec1cd5b429cbf7b7a600eeeef595a7b` — rank-suppression regression.

The current branch head at this report update is
`00b7b43eaec1cd5b429cbf7b7a600eeeef595a7b`, 330 commits ahead and two
behind main. Main remains unchanged; no merge was performed. CI has not
completed, and no PASS is claimed. The operator norm is a finite float64
comparison under the stated common-row embedding; it does not yet bound
recoupling, truncation, or continuum-limit errors.


## Follow-up: independent dense SVD oracle for the adjacent-cutoff norm (2026-11-11 15:00 UTC)

The synthetic adjacent-cutoff operator-norm regression now includes a second,
independently assembled dense oracle. It explicitly constructs the union of
row labels, forms the dense matrix for `C_current-C_previous`, computes its
singular values with NumPy's dense SVD, and compares the complete available
singular spectrum and leading operator norm against the streaming TSQR/SVD
path. The regression uses one-row streaming blocks to exercise the TSQR
merge path rather than a single dense block.

The independent oracle is intentionally a test-only small-matrix check; it
does not materialize production cutoff matrices. It verifies the sparse
row-union assembly and streaming norm on a controlled case, but it is not a
production calculation and does not bound upstream recoupling or truncation
errors.

The constraint-stability suite remains at 31 test methods (the oracle was
added within an existing method). The static total across the four focused
suites remains **89 methods**. Source inspection only; no suite execution is
claimed.

Latest test commit:
`09fea78072ebf7d09d3d4e773ebc22277fba1b6e` — add independent dense SVD
oracle to the adjacent-cutoff comparison regression.

At this update the research branch head is
`09fea78072ebf7d09d3d4e773ebc22277fba1b6e`, 341 commits ahead and two
behind main. Main remains `254f19fdab41d32b920bf454f987aa7c0fdba34d`;
merge base remains `23669209c83fe0aa2741612e5b2b119562ee6d7d`. No merge
was performed. Latest observed workflow runs remain queued with no
conclusion; no CI PASS or local execution is claimed.


## Follow-up: independent dense SVD oracle for adjacent-cutoff difference (2026-11-11)

The existing common-row adjacent-cutoff regression now performs a second,
independent assembly of the represented difference matrix. It constructs
the explicit union of nonzero row labels, materializes the small dense
matrix for the synthetic fixture, and compares the complete singular-value
list against the streaming TSQR/SVD result. The TSQR comparison uses a
one-row block size to exercise the streaming merge path rather than only a
single dense block.

This strengthens the test beyond checking one hand-computed operator norm:
it checks the nonzero singular spectrum of a rectangular row-embedded
difference. It remains a synthetic finite-dimensional regression; it does
not execute the production Peter-Weyl column producer or certify upstream
recoupling coefficients.

The change modifies an existing test method, so the static method inventory
is unchanged: 31 constraint-stability methods and 89 methods across the
four focused suites. The updated test has not been executed in this
environment. GitHub Actions remains queued; no CI PASS is claimed.


## Follow-up: distinguish full row union from nonzero difference support (2026-11-11)

The adjacent-cutoff comparison now reports three separate row counts:
the full union of labels appearing in either cutoff map, the intersection
of labels appearing in both maps, and the nonzero row support remaining
after coefficient-wise subtraction. This distinction matters because rows
whose coefficients are unchanged belong to the common embedding but cancel
exactly from the sparse difference operator. The synthetic fixture now
asserts union count 3, common count 2, and nonzero difference support 2.

For that fixture the represented difference is equivalent, up to row
ordering, to diag(0.5, 0.25). An independent local NumPy calculation gives
singular values (0.5, 0.25) and spectral norm 0.5. This checks the analytic
fixture only; it is not execution of the repository test or production
Peter-Weyl calculation.

The static inventory remains 31 constraint-stability methods and 89 methods
across the four focused suites. The test is committed but has not been run.
At this update the research branch is 345 commits ahead and two behind main;
main remains untouched. The latest observed Actions runs are queued with no
conclusion. No CI PASS or production execution is claimed.


## Follow-up: positive-metric kernel theorem and spectral bounds (2026-11-11 15:0x UTC)

For the direct-sum constraint map
`C psi = (H_0 psi, ..., H_(n-1) psi)`, the weighted master operator is
`M_G = C^dagger (G tensor I) C`. If the node metric `G` is Hermitian
strictly positive definite, then
`lambda_min(G) I <= G <= lambda_max(G) I`, and congruence by `C` gives

`lambda_min(G) C^dagger C <= M_G <= lambda_max(G) C^dagger C`.

Consequently, by the min-max principle, the ordered eigenvalues obey
`lambda_min(G) sigma_i(C)^2 <= lambda_i(M_G)
<= lambda_max(G) sigma_i(C)^2`, and
`ker(M_G) = ker(C) = intersection_v ker(H_v)`.
The metric can change the nonzero spectrum and condition number, but not the
exact common kernel. If `G` is only positive semidefinite, the lower bound
vanishes and this kernel-equivalence conclusion need not hold.

The constraint-stability module now implements
`positive_metric_spectral_bounds`. It:
- fails closed when direct TSQR reports scale-separation loss;
- validates finite, descending singular values and a finite Hermitian metric;
- rejects metrics that are not strictly positive definite;
- reports lower/upper eigenvalue bounds, metric condition number when
  representable, and lower-bound underflow;
- labels the result as an exact finite-dimensional matrix inequality, not
  an error bound for the upstream construction.

The result is integrated into `run()` for three deterministic positive
metrics (identity, diagonal, and random positive) and both legacy/no-prune
constraint maps. The random metric seed is recorded as `20261111`.

Three synthetic regressions were added:
1. Dense complex two-node constraint map: the actual eigenvalues of
   `C^dagger(G tensor I)C` lie within the predicted bounds.
2. The weighted master kernel projector agrees with the direct kernel
   projector for a rank-deficient C and strictly positive G.
3. Semidefinite/non-Hermitian metrics are rejected, and scale-separation
   warnings suppress the bounds.

An independent local NumPy calculation of the dense fixture gives
constraint singular values approximately
`(5.20810480, 1.69577249, 2.56e-16)`, metric eigenvalues
`(1.23095842, 2.16904158)`, and weighted-master eigenvalues
`(50.6686448, 5.73135518, 9.27e-16)`. All three ordered eigenvalues
satisfy the corresponding lower and upper bounds. This is a synthetic
oracle calculation, not execution of the repository test suite or a
production BQG calculation.

The constraint-stability suite now contains **34 test methods** by source
inspection. The four-suite static inventory is **92 methods** (31
master-projector integrity, 34 constraint stability, 13 streaming TSQR/SVD,
14 raw-volume operator). None of these counts implies test execution.

Commits for this stage:
- `636d5a0fbf51c6795adc79d9e561ed5f2d7a1450` — add positive-metric
  spectral-bound helper;
- `770c2e58a1f0827bea26a835d35b48c09890e4ba` — add dense-oracle,
  kernel-invariance, and metric-validation tests;
- `9f9c0b56f045522f72534c9f3fe2c408a84563b2` — integrate three positive
  metric bounds into the constraint audit;
- `dda2d4fe252df9bc38c5148639a09a0d9aea2dbd` — guard condition-number and
  lower-bound underflow reporting.

This theorem concerns the finite represented matrices. It does not establish
that the chosen finite-spin constraint maps converge to the intended
continuum operators, nor that the discretized metric is a faithful
approximation to a continuum rigging construction.
