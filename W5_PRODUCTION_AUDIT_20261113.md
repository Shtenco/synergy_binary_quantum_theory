# W5 production audit — 2026-11-13

## Scope and provenance

This audit targets the actual V09/V10/V11 W5 Peter–Weyl production path, not a replacement toy operator.

The source snapshot was obtained from the W5 source-staging artifact associated with PR #81:

- source commit: `e7f6b983e2990b1c268d9164f7f2e438f54e9e6e`
- source tree: `60e32745e0fb630964714876775e8cde20e1adab`
- 760 manifest-listed files passed byte-count and SHA-256 verification.
- The staged runner recorded Python 3.12.16, NumPy 2.5.4, SymPy 1.14.0, mpmath 1.3.0. The production audit workflows use Python 3.11 and currently install unpinned NumPy/SymPy; pinning the numerical environment is still required for bit-for-bit reproducibility.

## 1. L1 sparse-accumulator audit

GitHub Actions run `43444754151` completed successfully. The local reproduction of `tests.test_bqg_production_linearity` reports 9 tests, 6 passing and 3 expected failures.

- The strict reference `SparseLinearAccumulator` passes exact linearity/additivity, complex homogeneity, and order-independence tests.
- Legacy production helpers `COMP.add` and `CK.add_state` still prune a running partial sum. The expected failures deliberately preserve this known sub-threshold nonlinearity.
- The strict reference accumulator is **not yet wired** into the active V09/V10/V11 production call graph.

A passing L1 job therefore does not certify production linearity.

## 2. L2 real-production-primitive audit

GitHub Actions run `43444729954` completed successfully in 92.8 seconds: 12 tests, 9 ordinary passes and 3 expected failures.

The safe-domain complex homogeneity and independent-basis additivity checks pass for the real production primitives. The three expected failures characterize threshold-boundary additivity violations in the volume, `H_E`, and full covariant `C(V)` paths.

These are diagnostic expected failures, not acceptance of full linearity. A fixed matrix-level sparsification policy or a proven error-controlled accumulator must replace state-dependent intermediate pruning before a production linearity claim is made.

## 3. Full V10 W5 closure: reproducible blocker

GitHub Actions run `43444229133` executed the five W5 seeds and all five K5 source nodes at `Jmax2=10,12,14`, `max_states=500`, `max_depth=8`, `amp_tol=1e-12`, production pruning, cache on.

**All three cutoff jobs failed at the same guard before producing closure JSON.**

The stack ends in `peter_weyl_covariant_K_composition_cached_gate.restore_scalar_charge_M`, called by the sine-ordered H_E path. The reported transition is:

```text
node 0
input J2=0, canonical M2=0
output J2=2, M2=0
RuntimeError: scalar H_E violated canonical J,M preservation
```

This invalidates the current highest-weight-only scalar Wigner–Eckart reduction for this charged-sector path unless the underlying representation/selection rule is re-derived. Do not suppress the guard to make CI green.

It is not yet a no-go theorem for BQG or proof that the full physical Hamiltonian is non-scalar. The charged label may describe only the local four-leg sector while an external holonomy index carries additional SU(2) charge. The correct total representation and tensor rank must be derived.

## 4. Production trace resolves the scale of the apparent SU(2) violation

The production artifact from run `43446220696` was downloaded and inspected directly (artifact SHA-256: `38405a69db9bf0de71b2f214f89ac1a436e19b8bf9a121c94a47c76237fd83d9`).

For the first direct K call, K seed 0, node 0, Jmax2=10, input `J2=0,M2=0`:

- complete sine-`H_E` output support: 172;
- correct-sector weight: `6.641081547822338`;
- wrong-sector weight: `2.9364927759944046e-20`;
- wrong-sector fraction: approximately `4.42e-21`;
- relative wrong-sector amplitude: approximately `6.65e-11`.

The passive trace contains 780 individual sequence rows. Individual fixed-index sequences have large wrong-sector contributions, but after the full oriented-spec sum each completed spec is in the expected sector; the final accumulated wrong-sector residue is at floating-point scale. The current exact-key guard fails on the mere existence of this residue even though the downstream physical acceptance condition already allows a wrong-sector squared-weight fraction below `1e-18`.

This is a concrete **guard-contract mismatch**. It is not evidence that SU(2) is physically broken. It also does not yet prove that pruning is innocent: the paired shadow job is still required.

### Norm-bounded restoration patch

Commit `040d038dfd62c0627efa2cd01ac7212b912279b3` changes the cached wrapper to:

- compute per-charged-node total and wrong-(J,M) weights before cleanup;
- enforce the existing `1e-18` squared-weight fraction contract;
- raise on material leakage;
- discard only rows in a sector already certified below that bound;
- never relabel a wrong-J or wrong-M row as a correct row;
- retain per-call leakage diagnostics in the runtime memoization output.

Seven regression tests were added in `tests/test_peter_weyl_scalar_sector_restore.py`. The first CI attempt exposed a missing `math` import (runs `43448564400` and `43448588477`); the import was added in commit `7e4caa7a68374cac9c5f0b467118f8cf54e85d91`. The rerun result is pending.

A second independent gate, `scripts/peter_weyl_w5_direct_all_M_su2_gate.py`, now evaluates the actual sine-ordered (H_E) on all three magnetic components of a fixed charged (J=1) seed and measures ([J_+,H_E]) and ([J_-,H_E]). It runs both production and instrumented-shadow modes, bypasses memoization, and keeps the representation test independent of the highest-weight cache. Its CI result is pending.

The norm-bounded restore regression suite completed successfully in both matrix jobs of run `43448695465`: 7 tests passed in each job after adding the missing `math` import. The heavier context-probe step was then cancelled in the shadow job by runner shutdown; the production probe was still running at the latest check. A later direct all-M workflow run is queued/running; no covariance conclusion is claimed yet.

## 5. Discriminating diagnostics launched

Two diagnostics were launched on this isolated branch:

- Run `43445958758`: one-node W5 closure in `shadow_no_prune` mode, with all currently instrumented pruning sites bypassed. It was still executing at the latest status check; no artifact/result is available. The ledger itself still reports incomplete threshold-site coverage, so a passing shadow run would not yet prove zero loss.
- Run `43446220696`: paired `production` and `shadow_no_prune` executions of the existing `peter_weyl_w5_scalar_failure_context_probe.py`, restricted to the K=0 seed and source node 0. The `shadow_no_prune` matrix job was cancelled after the runner received a shutdown signal before producing an artifact; this is not a numerical result. The `production` matrix job was still executing at the latest status check. The probe preserves the guard and captures the exact input key, operator stage, per-sequence output J,M weights, support, and largest complex amplitudes.

The production probe is the next artifact to inspect. If the wrong-sector weight appears before pruning or persists in a completed shadow run, the representation/canonicalization or operator construction must be corrected. If it appears only after thresholded accumulation, the exact sequence trace will identify the first asymmetric pruning site. No conclusion is drawn from the cancelled shadow job.

## 6. Transitive pruning coverage inventory added

The prior W5 runner hard-coded `lower_level_pruning_ledger.coverage.complete=false`. The V10 import closure and every threshold-looking branch are now being independently inventoried by `scripts/bqg_pruning_coverage_inventory.py`. It walks the transitive local Python import graph, lists threshold comparisons with source lines, verifies each amplitude-pruning site uses the shared `PRUNE`/loss-ledger path, rejects unclassified sites and non-audit LRU caches, and separately classifies exact SU(2) rules, spectral zero-mode conventions, support caps, and residual-acceptance thresholds.

Seven existing modules already contain instrumentation for the important intermediate sites: primitive intertwiner projection, product coefficients, cached T-state outputs, complete charge projection, volume matrix/state outputs, Euclidean sine outputs, and Lorentzian local-volume columns. The new inventory makes that coverage auditable rather than inferred from a hand-maintained list.

The standalone CI run `43450087875` is queued. Only after it passes will the V10 runtime ledger be allowed to report static threshold-site coverage complete. Cache-on/off ledger equality, spin-wall/state/depth-cap accounting, and full complex support convergence remain separate open controls; even a complete static inventory is not a zero-loss theorem.

## 7. Required repair sequence

1. Inspect both artifacts from run `43446220696` and the shadow result from `43445958758`.
2. Derive the correct SU(2) action on the full charged state, including external holonomy/matrix indices. Decide whether the operator is rank-0 on the total coupled charge or a nontrivial tensor on the local charge label.
3. Build a finite full-M reference for one fixed W5 seed and node. Compare every magnetic block against the proposed reduced Wigner–Eckart implementation. Keep the existing guard active.
4. Only after the representation is correct, address the threshold-induced nonlinearity. Inventory every pruning site transitively; `lower_level_pruning_coverage_complete=false` remains a hard blocker.
5. Integrate coherent accumulation before a frozen, globally declared output-basis mask. Prove convergence and error bounds; do not replace intermediate pruning with an unbounded dictionary and call it solved.
6. Rerun V10 at Jmax2 10/12/14 with cache on/off and production/shadow comparisons. Require complete closure, no spin-wall omissions, no cap hits, zero unexplained scalar-sector leakage, direct-adjoint residual convergence, and full complex support comparison.
7. Only then rebuild V11 common carrier and the W5 joint-Gram candidate. The current W5 cone certificate is a reduced spectral lower bound, not a measured Lorentzian cancellation or physical kernel.

## Verdict

- L1 strict reference accumulator: **PASS**.
- L2 safe-domain production primitive controls: **PASS**, with three expected threshold-boundary failures.
- Legacy production linearity: **NOT CERTIFIED**.
- V10 full W5 closure at Jmax2 10/12/14: **BLOCKED** by charged-sector J,M preservation assertion.
- Full no-pruning coverage: **INCOMPLETE**.
- Physical Lorentzian W5 operator, joint kernel, and physical causal cone: **OPEN**.

No physical claim is promoted by these finite computational controls.


## 8. Follow-up: close the static cache-coverage blind spot (2026-11-13)

The first end-to-end inventory execution exposed a real control-flow defect: cache decorators were inspected only after the function's threshold-comparison scan, so a cached function with no threshold comparison could bypass the cache audit through an early `continue`. The scan now inspects decorators for every function before that early exit. A dedicated fixture test constructs exactly such a function and requires the inventory to fail when it finds an ordinary, unclassified LRU cache.

The next production run also identified two cases requiring semantic classification rather than blanket exclusions:

- `bqg_pruning_audit.cached` is the private `functools.lru_cache` implementation inside `audit_lru_cache`. The wrapper stores the numerical result together with the ordered pruning-event trace, and replays that trace on cache hits. It is now reported separately as `audit_lru_cache_internal_value_and_trace_cache`; dedicated cold-hit, warm-hit replay, and cache-bypass regression tests were added.
- `bqg_sparse_linear_accumulator.sparsify_matrix_columns` applies a mask only after matrix columns have been assembled, so it is not state-dependent intermediate amplitude pruning. Its `tol < 0` guard is classified separately as input validation; `abs(z) > tol` is classified as a post-assembly fixed operator mask. **This classification does not certify the approximation:** an operator-norm/error bound for the discarded matrix remains required before this mask can support a quantitative production claim.

GitHub Actions run `43492745470` completed successfully on commit `26e1bd75bc8244c6feb74e1e35173f74be3ec57d`:

- transitive local import closure: 19 modules;
- threshold-like comparisons inventoried: 78;
- unclassified sites: 0;
- unclassified ordinary LRU caches: 0;
- the internal audit-aware cache was explicitly identified;
- all four then-current inventory regression tests passed.

This is a **static inventory PASS only**. It establishes coverage of the declared V10 local import closure and classification of discovered sites; it does not establish zero runtime loss, a bound for the post-assembly matrix mask, convergence, W5 closure, or physical SU(2) covariance. The strengthened follow-up run `43492874165` completed successfully on code commit `3fe0e9e83499d85be78af7c989df99a9b285973f`: all four inventory tests and all eleven pruning-audit collector tests passed. The run reports 78 threshold-like comparisons, zero unclassified sites, zero ordinary unclassified LRU caches, and one explicitly classified audit-aware internal cache. The subsequent precision refinement separated tolerance validation from the fixed matrix mask and is covered by that run's four inventory tests.


## 9. Runtime W5 ledger now reports static inventory without overstating completeness

The W5 closure runner now imports and executes the transitive inventory and embeds a compact result under `lower_level_pruning_ledger.coverage.static_inventory`. It deliberately keeps the legacy overall `coverage.complete=false` and top-level `lower_level_pruning_coverage_complete=false`: static classification is not a runtime zero-loss certificate.

The new status is `STATIC_THRESHOLD_INVENTORY_PASS_RUNTIME_LOSS_CONTROLS_PENDING` when the inventory is complete. The remaining runtime controls are explicitly enumerated:

1. A certified operator-norm or induced-norm error bound for the post-assembly fixed matrix mask.
2. Cache-on/cache-off pruning-ledger equality across the full composition chain.
3. Separate, complete accounting for spin-wall and state/depth-cap omissions.
4. Full complex-support convergence across spin cutoffs and amplitude-threshold ladders.

The `shadow_no_prune` result remains marked incomplete. Its bypass removes only wired nonzero amplitude thresholds; it does not remove the post-assembly matrix mask or support caps. This prevents a static inventory PASS from being misreported as a full no-pruning reference.

This integration is committed, but its new runtime output is awaiting the W5 closure and shadow diagnostic workflows triggered by the latest code changes. No closure or physical covariance result is inferred from the static inventory tests.


### Latest inventory after runtime-ledger integration

Run `43493076092` completed successfully on inventory commit `743a8825afdf76e0e3f73a1a0e2f65551cf2c9eb`, after the W5 closure runner began importing the inventory module. The closure now includes the auditor in its own declared dependency set:

- transitive modules: 20;
- threshold-like sites: 78;
- unclassified sites: 0;
- ordinary unclassified LRU caches: 0;
- inventory tests: 5 passed;
- pruning collector/cache replay tests: 11 passed.

This confirms the static inventory remains complete when the inventory itself is part of the W5 dependency closure, and confirms the ledger-builder's separation between static classification and overall runtime completeness. The W5 cutoff-closure run triggered by the same integration is still pending; its numerical closure result is not inferred from this inventory PASS.


## 10. Certified error bound for the standalone post-assembly matrix-mask utility

The standalone helper `bqg_sparse_linear_accumulator.sparsify_matrix_columns_with_error_bound` now returns both the masked columns and an explicit loss report. For the discarded matrix (E=A-A_\tau), it computes the Frobenius norm and records the rigorous inequality

`||A-A_tau||_2 <= ||A-A_tau||_F`.

The report contains the absolute spectral-norm error upper bound, relative Frobenius error, retained/discarded nonzero counts, structural-zero counts, per-column discarded norms, and an explicit scope limitation. The implementation uses iterative `math.hypot` accumulation to avoid directly squaring large finite magnitudes.

The new test suite checks:
- the exact discarded Frobenius norm on a complex-valued fixture;
- the induced output-error inequality for a test coefficient vector;
- consistency with the legacy mask API;
- empty/zero matrix behavior;
- rejection of invalid tolerances and non-finite matrix entries.

**Scope qualification:** the module documentation explicitly states that this sparse accumulator is not wired into the V09/V10/V11 production path. The W5 closure imports its `SparseLinearAccumulator` for coherent summation, but the W5 closure runner does not call `sparsify_matrix_columns` or its new certified companion. Therefore the bound API is available and tested, but the runtime ledger correctly reports `fixed_operator_mask_bound_api_available=true`, `fixed_operator_mask_bound_applied_in_production=false`, and `fixed_operator_mask_error_bound_complete=false`. This avoids treating a tested utility as evidence that a production operator has been sparsified with a certified error.

GitHub Actions run `43498074138` completed successfully on the corresponding code/workflow commit `5426d0f2dc0419ac536b9f7557f510476e210b87`:
- 20 modules in the static dependency inventory;
- 80 threshold-like sites;
- zero unclassified sites;
- 5 inventory tests passed;
- 11 pruning-audit/cache tests passed;
- 5 sparse-matrix error-bound tests passed.

The final small test-only follow-up is run `43498098817` on branch HEAD `2d5639b9e57c8446ca867e54a3635a4297a5b815`; it is being checked separately. Static inventory PASS and helper tests do not establish the pending W5 closure, full no-pruning reference, cache-on/cache-off equality, or physical SU(2) covariance.


### 10.1. Final verification of the certified matrix-mask helper

A first overflow-regression fixture used two entries of (10^{308}), but the stable `hypot` accumulator correctly represented their combined norm as a finite value. The test expectation—not the implementation—was wrong. The fixture was tightened to two finite entries of (1.7\times10^{308}), whose combined Frobenius norm exceeds the finite floating-point range. The helper now rejects this case rather than emitting an infinite or misleading finite certificate.

The corrected CI run `43498364038` completed successfully on the branch HEAD commit `4aad457123c99d533c64cffbc2da22815eee3c5d`:

- 20 modules;
- 80 threshold-like sites;
- zero unclassified sites;
- 5 inventory tests passed;
- 11 pruning collector/cache tests passed;
- 6 matrix sparsification/error-bound tests passed, including the overflow rejection.

The separate production L1 audit `43498254862` also completed successfully on the implementation commit `0744524e459e42d7b3907effd5196c02d259fd6c`; both module compilation and the L1 characterization/strict-reference test step passed.

The helper's bound is now numerically regression-tested, but remains a standalone utility. No claim is made that the actual W5 production closure currently applies this matrix mask or uses its error report.


## 11. Dedicated cache-on/cache-off equivalence runner

Added `scripts/bqg_v10_w5_cache_equivalence.py` and workflow `.github/workflows/bqg-w5-cache-equivalence.yml`. The runner launches two isolated fresh Python processes with identical W5 parameters, changing only the context-local cache mode (`on` versus `off`). It refuses to claim equivalence unless both processes produce complete finite closure matrices.

The comparison checks exact support and complex coefficients for:
- the forward W5 matrix;
- the independently ordered direct-adjoint matrix;
- the ordered lower-level pruning-event ledger (excluding only cache-policy metadata);
- the coefficient threshold and closure-cap loss ledger;
- the coherent node-sum ledger;
- spin-wall and closure diagnostics;
- the spectral sensitivity ledger.

A timeout or capped closure yields `INCOMPLETE_CACHE_EQUIVALENCE`, not PASS. Any exact matrix or ledger mismatch yields `FAIL_CACHE_EQUIVALENCE`. The workflow is currently running as Actions run `43498621358`; no result is claimed until its uploaded JSON and completed job are inspected.

The full-node closure study and shadow no-prune diagnostic remain separate gates. Their currently queued/running jobs are not interchangeable with the cache comparison and do not yet establish a no-pruning reference or physical SU(2) result.


### 11.1. Workflow execution status and non-results

The first cache-equivalence run `43498621358` compiled the comparison runner successfully and entered the two-process matrix comparison step; it has not yet produced a completed report. The follow-up workflow revision that adds explicit regression tests for incomplete-run summaries is queued as run `43498807251`.

The shadow diagnostic run `43498254640` ended with its numerical step marked `cancelled`; its upload step was skipped, so no diagnostic JSON exists for inspection. This is **not** a numerical or physical FAIL. The full-node closure run `43498254658` remains pending and has not yet exposed matrix jobs. Consequently, neither the full closure nor the no-pruning reference has a confirmed numerical result at this point.


### 10.2. Explicitly scoped W5 seed subsets

The W5 closure runner now accepts `--seed-K` as an explicit subset of the canonical labels `(0, 2, 4, 6, 8)`. Omitting the option retains the complete five-seed production run. Every result records both `seed_count` and `seed_Ks`; bounded one-seed diagnostics therefore cannot be confused with the full W5 seed family.

The shadow smoke workflow now requests only `--seed-K 0`, with `max_states=1) and `max_depth=0`, and asserts that the result reports exactly `[0]`. This is deliberately a scoped execution diagnostic, not a five-seed W5 closure or a physical validation.

The closure workflow additionally runs the dedicated algebra/closure-control unit-test module before launching its full five-seed, five-source-node matrix over (J_{max 2}=10,12,14). This ensures the new seed-subset validation and existing direct-adjoint algebra tests execute in CI.


### 10.3. CI dependency correction for the W5 algebra unit suite

The first attempt to add `test_bqg_v10_w5_lorentzian_adjoint_closure.py` to the fast inventory workflow failed for an infrastructure/configuration reason: that workflow did not install NumPy, while the algebra test module imports NumPy at module load. The existing inventory, pruning-audit, and sparse-mask suites still passed in that run (5, 11, and 6 tests respectively); the newly added algebra suite did not execute.

The workflow has now been corrected to install `numpy` and `sympy` before running the W5 closure unit suite. This is a test-environment fix, not a change to the test assertions or the operator implementation. The rerun is required before claiming the new suite passes.
