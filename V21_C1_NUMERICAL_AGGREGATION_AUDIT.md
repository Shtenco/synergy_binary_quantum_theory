# V21 C1 numerical aggregation audit

Status: code-level comparison of V21 staged aggregation with the historical V2 ordered-triple implementation. No physical computation was run by this audit.

## 1. Historical pruning points

The historical ordered-triple path uses `TOL = 1e-11` and `add(dst, src, scale, tol)` when combining the eight auxiliary `(i,j,k)` contributions. The add operation removes a key whenever the current partial sum has magnitude less than or equal to the threshold. A later auxiliary contribution can then re-create the key from zero.

The production V2 path is `pl_lorentzian_group_worker_v2.py` (or the one-term V2 worker), which wraps the shared group worker with the tetrahedral backend. The base group worker sets `TOL = W.TOL = 1e-10` and passes that value explicitly to `LGA.project_scalar_gauss`; this overrides the module default `LGA.TOL = 1e-11`. The V21 collector therefore uses `W.TOL = 1e-10` for the matched production V2 comparison. The inner ordered-triple auxiliary accumulator still uses `LP.TOL = 1e-11`.

This ordered thresholding is not algebraically equivalent to summing all eight contributions first and pruning once. Example: two contributions of magnitude `0.75e-11` to the same key are each dropped by the historical incremental accumulator, while an unpruned sum retains `1.5e-11`.

## 2. Two comparison lanes

The V21 collector now calculates both:

- **Unpruned lane:** stable complex summation using separate real/imaginary `math.fsum`; no support threshold during auxiliary aggregation; then scalar Gauss projection at `LGA.TOL`.
- **Legacy-aligned lane:** ordered incremental accumulation reproducing the historical `add(..., tol=LP.TOL)` behavior; then scalar Gauss projection at `LGA.TOL`.

The four independent V2 references are gated against the legacy-aligned lane at relative delta below `1e-10`. The legacy-aligned lane reproduces `LP.TOL` auxiliary accumulation followed by the production V2 `W.TOL` scalar projection. Unpruned-versus-V2 deltas are separately reported. Both lanes must satisfy the nonscalar rejection gate.

The final V21 physical pair is assembled from the unpruned lane, with stable summation and no tolerance pruning. Its `VERIFIED_EXACT_ZERO_PAIR` label means zero in the resulting floating-point complex sum only; it is not a symbolic proof of exact cancellation.

## 3. Why this matters

If only the legacy-aligned lane is examined, a sub-threshold nonzero residue can disappear by construction. If only the unpruned lane is compared to a legacy V2 artifact, threshold-induced differences can be misdiagnosed as an operator mismatch. Reporting both lanes separates these two questions.

## 4. Remaining caveats

- The two lanes still inherit floating-point arithmetic inside the underlying SU(2) operator kernels.
- The unpruned lane is not arbitrary precision and does not establish a rigorous interval bound.
- The independent V2 artifact metadata must remain externally pinned and semantically matched.
- The 32-path computation, four V2 comparisons, and all structural tests still require actual execution on a working runner.
- A successful software comparison is not a proof of the full Hamiltonian constraint algebra, physical TT closure, or Einstein dynamics.

## 5. Underflow and overflow handling

State and pair L2 norms use iterative `math.hypot`, not `sqrt(sum(abs(z)**2))`. This prevents a nonzero amplitude such as `1e-200` from being mislabeled as zero because its square underflows. Pair zero classification is based on whether the assembled amplitude map is empty after exact floating-point zero removal, not whether a squared norm happens to evaluate to zero.

Relative V2 comparison no longer uses an arbitrary `1e-15` denominator floor, which could conceal a complete mismatch between tiny states. If any comparison norm or delta overflows to a non-finite value, the relative delta is set to infinity so the acceptance comparison fails closed. Projection and final pair norms are also required to be finite.

## 6. Production V2 workflow and backend identity

The first pass inspected the legacy standalone `pl_lorentzian_triple_worker.py` and initially mistook it for the canonical independent V2 path. That was the wrong workflow target. The repository's production V2 workflows use `pl_lorentzian_group_worker_v2.py` and `pl_lorentzian_single_worker_v2.py`; both already wrap their base computation in `install_tetrahedral_volume_backend()`. Therefore the production V2 path already uses the same tetrahedral backend as V21.

The repair now adds explicit source-commit and backend-SHA provenance to those production V2 workers and requires the V2 collector and V21 collector to verify those values. The legacy standalone worker is not a substitute for the 48-term production V2 artifacts. V2 references generated before the provenance fields were added remain historical/unattested until regenerated.

The production V2 worker's projection tolerance is `W.TOL = 1e-10`, passed explicitly into `LGA.project_scalar_gauss`; it is not the module default `LGA.TOL = 1e-11`. V21 uses the same explicit production V2 tolerance for its matched comparison lane.
## 7. Strict JSON integer metadata

A source-level adversarial review found a Python/JSON type pitfall: bool is a subclass of int, and equality treats True == 1 and False == 0. Membership checks such as value in (0, 1) and integer coercion such as int(manifest["source_node"]) can therefore accept malformed boolean selectors.

The V21 pair collector now rejects booleans and non-integer values in all three stage path metadata records for source_node, omit, cycle, Jmax2, and each of the three auxiliary indices. It also requires the three stage metadata records to describe exactly the same mode, order, and auxiliary tuple. Top-level manifest selectors are type-checked before use rather than coerced.

Four path-metadata and three manifest-selector regressions were added (tests 50–56). A further guard now requires checkpoint state_support to be a nonnegative JSON integer rather than accepting bool/float equality with a support count; test 57 covers this. The stage-to-stage predecessor gate also requires a positive non-boolean integer support count before a worker may consume the prior checkpoint (test 58). The repository suite now contains 58 tests. These tests have not yet executed in CI: GitHub Actions runs 42470000116 (attempts 1 and 2), 42470158108, 42470338726, and 42470819556 failed before runner assignment (runner_id=0, empty step lists); the log endpoints are unavailable. This is an infrastructure-level no-execution result, not a test assertion failure.

The independent local synthetic harness remains separate: its 17/17 PASS result does not validate the production collector or replace the 58 repository tests.
## 8. V2 48-term artifact binding and orbit validation

The independent V2 path also needed a stricter artifact boundary. The V2 single-term and six-term group workers now write the SHA-256 of each exact NPZ beside its JSON metadata. The V2 48-term preflight requires that digest to match the sibling NPZ before the historical reducer can load it. It also validates the NPZ key set, fixed spin/K shapes, integer index dtypes, complex finite amplitude vector, unique state keys, and exact integer support count before the reducer can deserialize it.

Each V2 term record is now checked for literal boolean passed=true, an exact integer index in 0..23, exact integer source_node=0, an exact integer orientation coefficient in {-1,+1}, Jmax=3.5, and a nonempty checks map whose values are all literal boolean true. The collector verifies the complete unique forward/adjoint 24+24 orbit, recomputes each PL orientation coefficient and ordered target tuple from the canonical worker specification, and enforces one source commit and one backend SHA across all 48 terms.

The V21 collector additionally requires each independent V2 metadata npz_sha256 to equal the externally pinned NPZ SHA-256. This closes the metadata-to-state binding seam for the four references used in V21 comparison.

Ten more regressions were added (tests 59–68): boolean passed/index/source/coefficient rejection, NPZ digest mismatch and schema rejection, duplicate-key rejection, support-count typing, canonical V2 term acceptance, and V21 metadata-to-NPZ digest binding. The repository test inventory at this stage was 68 tests. These tests remain unexecuted because GitHub Actions jobs continue to terminate before runner assignment.

These checks validate provenance, path identity, and file integrity. They do not regenerate the four V2 references, run the 48 physical V2 terms, or complete the 32-path V4 calculation.
## 9. Exact V2 acceptance-check inventory and strict reader metadata

A further source review found that checking only whether a V2 `checks` mapping was nonempty and all its supplied values were `true` was insufficient: a malformed artifact could omit a required acceptance criterion and still pass. The V2 production worker emits a fixed nine-key inventory. Both the 48-term preflight and the V21 independent-reference reader now require the exact key set and literal `true` for every value.

The V21 reader also rejects string-valued `Jmax` (only a finite numeric 3.5 is accepted) and requires `gauss_support` to be a non-boolean integer exactly matching the number of NPZ amplitudes. Four regressions were added (tests 69–72), including rejection of incomplete check maps at both collector boundaries. The repository suite now contains 72 tests.

A second review found a malformed indentation in the V2 six-term group worker's metadata rewrite block. This syntax error was corrected in `scripts/pl_lorentzian_group_worker_v2.py`; the workflow's `py_compile` step includes that module. However, GitHub Actions run 42495070894 also ended before exposing any job steps, so the compile and test commands are still not verified as executed.

The local archive was rerun after updating its isolated V2 guard copy: 32/32 synthetic tests pass under Python 3.13.5. This remains a separate synthetic result, not a production-module import or the 72-test repository suite.
## 10. V2 diagnostic completeness and numeric type safety

Review of the actual V2 group worker found another fail-open path: missing leakage fields were read with `dict.get(key, 0.0)`, so absent diagnostics could silently look like perfect zero leakage. The production evaluator now calls the dependency-light `scripts/v2_term_guards.py` validator. It requires all four operation-specific fields (`CV_complete_basis_leakage`, `CK_outer_complete_basis_leakage`, `CK_internal_volume_sector_leakage`, `CK_complete_charge_basis_leakage`), rejects negative/non-finite values, and validates both projection norms as finite and nonnegative.

Numeric metadata is type-checked before conversion: bool and numeric strings are rejected rather than silently treated as 1/0 or parsed numbers. The same strictness applies to `accepted2` and `rejected2`. The helper was separated from the heavy group-worker imports so the regression tests can import the guard without loading the full Lorentzian engine.

Six new repository tests were added in tests 73–78: missing diagnostics, negative leakage, non-finite charge diagnostic, negative/non-finite projection norms, boolean diagnostic, numeric-string diagnostic, and boolean projection norm (the negative/non-finite projection cases are combined in one test). Repository suite count is 78.

The CI workflow now compiles both `v2_term_guards.py` and `pl_lorentzian_group_worker.py`, in addition to the V2 workers and collectors. It uses `ubuntu-latest` rather than Windows to reduce hosted-runner cost. Nevertheless, the latest Actions jobs still fail without any step or log records, so neither compilation nor the 78 repository tests are confirmed as executed.

The local synthetic archive now contains 40/40 passing checks, including five diagnostic-guard tests and three additional numeric-type edge cases. These are isolated copies of the small validator, not production module tests.
