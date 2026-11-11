# BQG research continuity checkpoint — 11 November 2026, 00:20 UTC

**Repository:** `Shtenco/synergy_binary_quantum_theory`  
**Active branch:** `research/lemaitre-bqg-causal-cone-20261110`  
**Baseline:** `7c4bd56c7d78bdba333246a7910f6161cac7d28a`  
**Branch HEAD at checkpoint preparation:** `2a32384a8cce3b746242cb67faa3ca1ba117655f`  
**Main branch:** not modified; no merge performed.

This file is the continuity ledger for the work consolidated on 11 November. It records implementation facts, scientific evidence boundaries, CI state and the next dependency-ordered tasks. It does not claim unexecuted tests passed.

## 1. What is now preserved in the repository

### 1.1. Earlier W5 Lorentzian line

- V09 raw W5 column: Euclidean (A), raw Lorentzian (B), 24 epsilon triples, scalar-projection and leakage diagnostics.
- V10 two-sided adjoint closure: separately ordered direct-adjoint construction, closure over forward and adjoint support, matrix-level adjoint residual.
- V10 acceptance hardening: PASS is conjunctive across nonzero forward/adjoint norms, adjoint residual, scalar projection, leakage, spin wall and LSINE's explicit `passed`.
- V11 common carrier: coherent node-local accumulation before one shared threshold, V09/V10 forward comparison, (B=B_{sym}+B_{anti}) decomposition, basis/provenance hashes.
- Joint-Gram hardening: basis identity checks, explicit coefficients and normalization source, rejection of empty/nonfinite inputs, structural right-kernel padding for rectangular SVD.
- Hermitian ordering ledger: sign convention for (K_{comm}), two-K-leg conversion and the distinction between relative classical coefficient and unknown absolute W5 normalization.
- Physical causal-cone audit: independent binary operational influence, causal order, retarded Green support, spectral reconstruction and universality requirements.

### 1.2. New production-linearity line

- `BQG_PRODUCTION_LINEARITY_AND_CLOSURE_AUDIT_20261111.md`: active call graph, threshold classification, L1/L2/L3 tests, closure-loss diagnostics and cutoff convergence protocol.
- `scripts/bqg_sparse_linear_accumulator.py`: strict reference accumulator with no state-amplitude pruning; fixed matrix sparsification is a separate post-assembly operation.
- `tests/test_bqg_production_linearity.py`: coherent subthreshold accumulation, complex additivity/homogeneity, permutation stability, fixed-mask linearity, nonfinite input rejection, and expected-failure characterization of legacy thresholded helpers.
- `.github/workflows/bqg-production-linearity-l1.yml`: dedicated scoped workflow, push-triggered by changes to the reference helper/tests/workflow and manually dispatchable.

The reference accumulator is **not yet wired into the V10 KKV production path**. Its presence is a test oracle and refactoring foundation, not a production fix.

## 2. Verified source-level finding

The active V10 KKV route has state-dependent thresholded accumulation in:
- `peter_weyl_covariant_composition_gate.py:48-54` (`COMP.add`);
- `peter_weyl_covariant_K_leg_gate.py:53-60` (`CK.add_state`);
- `peter_weyl_covariant_K_leg_gate.py:89-134` (branch amplitudes containing the input amplitude);
- `peter_weyl_covariant_K_composition_gate.py:43-51` (basis-column contributions multiplied by actual input coefficients and passed through thresholded accumulation);
- `peter_weyl_covariant_K_composition_gate.py:54-99` (K-leg composition);
- `peter_weyl_covariant_volume_leg_gate.py:98-145` (covariant branch projection).

A fixed mask on precomputed matrix elements is a different operation and can define a fixed approximate linear matrix. Thresholding a running state amplitude after each scaled basis-column contribution is generally nonlinear and order-dependent.

The V10 matrix adjoint residual therefore remains a useful matrix-level consistency diagnostic, but is not sufficient to prove that the actual thresholded production function is a linear operator on arbitrary superpositions.

Adjacent legacy helpers in `k5_peter_weyl_safe_hda_column.py` were not attributed to V10 unless a direct caller was confirmed. This prevents over-reporting the active defect surface.

## 3. L1 characterization semantics

The legacy helper tests are deliberately decorated with `unittest.expectedFailure`. This means the current known defect is recorded without making the whole characterization workflow red. If a legacy test unexpectedly passes, unittest reports an unexpected success and the run fails. After strict accumulation is integrated into the active production path, remove `expectedFailure` so the tests become hard acceptance gates.

The strict accumulator:
- stores all nonzero complex contributions by output key;
- uses `math.fsum` independently for real and imaginary parts;
- does not delete a state key because an intermediate partial sum is small;
- applies optional sparsification only to already assembled matrix columns;
- rejects nonfinite inputs and scaled/accumulated values.

No local or remote PASS is claimed until the corresponding workflow finishes and its logs are inspected.

## 4. Commit ledger

| Commit | Change |
|---|---|
| `a4fb1cfdfc5a5c99415a77ccb0eba933f2ed20fc` | V11 coherent-sum threshold fix |
| `44feba50f09dbe18e03237d93e9dcaf730c8db8a` | V10 conjunctive LSINE/spin-wall/adjoint gate |
| `4d328b9c934f97d055cd679e5c0348261e9f1523` | V10 acceptance contract documentation |
| `46f6d97168df2aeb9a2a218af7687f52abfc80d9` | Consolidated production-linearity and closure-loss audit |
| `248d8e78f918d3d26fca37cea98c3480eb8076cd` | Pinpoint active V10 pruning call sites |
| `e9d87b41668bc530675d5ab46f4a5b22be328a42` | Strict sparse accumulator, L1 tests and initial workflow |
| `2a32384a8cce3b746242cb67faa3ca1ba117655f` | Push-triggered L1 workflow |

## 5. CI snapshot at 2026-11-11 00:20 UTC

| Run | Meaning | Current status | Artifacts |
|---|---|---|---|
| [42929934053](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42929934053) | L1 accumulator audit at HEAD `2a32384` | queued | none |
| [42908538520](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908538520) | V10 cutoff Jmax2=10/12/14 | queued | none |
| [42908462865](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462865) | Earlier V10 cutoff Jmax2=10/12/14 | queued | none |
| [42908462954](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462954) | Adapter tests | queued | none |
| [42908462792](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908462792) | V11 common carrier | queued | none |
| [42908538513](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908538513) | Adapter tests on earlier head | cancelled | none |
| [42908538772](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42908538772) | V11 common carrier on earlier head | cancelled | none |

Queued means no conclusion yet; it is not a PASS and not a FAIL. The cancelled jobs are not failed scientific tests; they were cancelled runs. No numerical artifact or physical conclusion is inferred from any of these states.

## 6. Scientific evidence boundaries retained

1. V09/V10/V11 are structural finite-sector computations, not a completed physically normalized Lorentzian constraint.
2. V10's closure is threshold-relative; current `omitted_support` does not count all transitions removed by `amp_tol`. Scalar-projection rejected norm is a distinct loss channel.
3. The absolute W5 normalization remains open; do not fit it to a desired kernel.
4. A finite W5 nullspace is not automatically the full physical constraint kernel.
5. The reduced-model dispersion and front-velocity formulae are not the full BQG causal cone.
6. A physical causal cone requires independent agreement of operational binary influence, retarded propagation and spectral metric, with a clock and volume calibration.

## 7. Dependency-ordered next work

1. Inspect L1 run 42929934053: job status, each step, logs, expected failures and artifacts.
2. Confirm the strict reference tests pass and that the legacy characterization tests fail exactly as expected.
3. Add L2 tests for the actual active primitives, under the same zero-aware volume patch/cache conditions as V10.
4. Add L3 comparison between direct production action on arbitrary complex superpositions and matrix action assembled from basis columns.
5. Refactor active KKV path to use no state-amplitude pruning after arbitrary input coefficients are applied; preserve declared fixed matrix masks separately.
6. Add below-threshold forward/adjoint counts, L2 loss and max discarded amplitude to closure diagnostics; track internal stage losses separately.
7. Re-run closure at `amp_tol = 1e-8, 1e-10, 1e-12, 1e-14` and strict reference where feasible; compare spectra only on a common embedded basis.
8. Re-evaluate matrix adjoint consistency after production-linearity is established.
9. Continue physical normalization, physical-kernel and causal-cone work as separate gates.

**Current scientific priority:** prove that the finite-sector production KKV map is linear on its declared domain before interpreting its basis-column matrix as the operator whose adjoint, spectrum or kernel is physically meaningful.


## 8. L2 source audit extension — 11 November 2026, 00:32 UTC

A second source inspection was completed after the original checkpoint. The new document [BQG_L2_PRODUCTION_PRIMITIVE_LINEARITY_AUDIT_20261111.md](BQG_L2_PRODUCTION_PRIMITIVE_LINEARITY_AUDIT_20261111.md) records the exact map-by-map L2 test matrix and L3 reconstruction contract.

Additional source-confirmed boundaries:

- `peter_weyl_covariant_K_composition_gate.py:43-51`: `apply_HE_complete_state_custom` multiplies a unit-key result by an arbitrary input coefficient and passes it through thresholded `COMP.add`. This is a direct arbitrary-superposition nonlinearity risk.
- The same module's `apply_K_complete_custom` combines `V H_E` and `H_E V` components through further thresholded accumulation.
- `peter_weyl_covariant_K_leg_gate.py`: `add_state` thresholds branch sums; `apply_HE_complete_key` uses it while assembling oriented/adjoint branches; `apply_volume_complete_state` applies a final (10^{-11}) output filter; `apply_HE_complete_state` applies a final (10^{-10}) output filter.
- `peter_weyl_covariant_K_composition_gate.py` therefore contains distinct loss boundaries at primitive branch assembly, arbitrary coefficient composition, volume action, and commutator combination. They must be measured independently.
- V10's outer `amp_tol` pruning occurs after coherent node-map summation, but the removed forward/adjoint amplitudes are not currently counted. This external loss must not be merged with internal primitive losses or scalar-projection rejection.

The new L2 test plan prioritizes real production maps: `CK.apply_volume_complete_state`, `CK.apply_HE_complete_state`, `COMP.apply_HE_complete_state_custom`, `COMP.apply_K_complete_custom`, `COMP.C_K_component`, and `COMP.C_volume_component`. It requires additivity, complex homogeneity, permutation/cancellation tests, direct-superposition versus frozen basis-matrix reconstruction, and stage-specific loss ledgers.

**New commit:** `4f888c8658299af441097674fa08a6949a837be9` — add L2 production-primitive linearity audit.

**CI status remains unchanged:** runs `42930769907` and `42929934053` were observed queued, with no artifacts. The existing L1 workflow does not execute the new L2 test matrix. No L1/L2/L3 PASS is claimed.

**Next gate:** implement and execute the three highest-priority L2 tests before changing production semantics. Do not wire the reference accumulator into V10 until the exact insertion boundary and regression behavior are established.


## 9. Executable L2 gate added — 11 November 2026, 01:00 UTC

Added `tests/test_bqg_l2_production_primitives.py` and `.github/workflows/bqg-production-linearity-l2.yml`.

The L2 suite initializes a real Peter–Weyl basis state, converts it through the repository's covariant/complete adapters, applies the preregistered zero-aware volume patch and calls three actual production primitives: `CK.apply_volume_complete_state`, `CK.apply_HE_complete_state`, and `KCOMP.apply_HE_complete_state_custom`.

Each primitive has:
- a hard complex-homogeneity test at a safe scale above its output threshold;
- an expected-failure threshold-boundary additivity characterization, calibrated from the actual unit-state output so the test compares two individually subthreshold contributions with their coherent sum.

The workflow compiles the production modules, runs the suite on Python 3.11 and uploads the test log even if a test fails. The suite does not yet test the full `COMP.C_K_component` / `COMP.C_volume_component` maps. This is an executable L2 characterization, not evidence of full production linearity or physical validity.

**New commits:** `340f311f9eba27ff0e53c3a0255ca15f39f2d206` (L2 primitive tests), `060b50dca818e2552067e5f905ccf240a4da6b85` (L2 workflow), followed by documentation updates.

**CI status:** L2 run has only just been triggered by the scoped push. Do not claim PASS/FAIL until its job steps and uploaded log are inspected. The previous L1 runs remain separately tracked.

**Next:** inspect L2 safe-domain homogeneity and each expected-failure reason; then extend to the full covariant maps and L3 direct-superposition/matrix reconstruction.


## 10. L2 canonical-suite reconciliation and CI poll — 11 November 2026, 01:05 UTC

A duplicate L2 test file was briefly introduced during implementation, then removed after a repository recheck found the canonical suite already present. The duplicate was deleted in commit `37dec6c0f5f282c7d893f9b45b8c427b7030e121`; no production source or canonical test was changed by that cleanup. The branch now retains one authoritative L2 test file and one L2 workflow.

Canonical implementation:
- `tests/test_bqg_l2_production_primitives.py` (145 lines): genuine Peter–Weyl covariant/complete basis fixture; applies the zero-aware volume convention and clears affected caches.
- Hard safe-domain complex-homogeneity checks for `CK.apply_volume_complete_state`, `CK.apply_HE_complete_state`, and `KCOMP.apply_HE_complete_state_custom`.
- Three `unittest.expectedFailure` threshold-boundary additivity characterizations. Their scale is calibrated from the actual unit-state output so individual contributions lie below the relevant output threshold while their coherent sum crosses it.
- Empty fixture outputs are skipped by the characterization helper rather than silently counted as linearity evidence; safe-domain homogeneity tests also skip on empty reference output. The L2 report must therefore distinguish PASS, XFAIL, SKIP and ERROR and must inspect fixture reachability.

Workflow:
- `.github/workflows/bqg-production-linearity-l2.yml` compiles the active production modules and test, executes the canonical suite on Python 3.11, and uploads `l2-production-primitives.log` even when tests fail.
- L2 run `42935840345` (run #1, HEAD `060b50dca818e2552067e5f905ccf240a4da6b85`) — job `130245591568`: queued; no artifacts.
- L2 run `42935968925` (run #2, HEAD `f1a7c8d1941448b8453b6a895826f2a37658c9a5`) — job `130246013044`: queued; no artifacts.
- L1 run `42930769907`, job `130228642399`: queued; no artifacts.
- L1 run `42929934053`, job `130225807106`: queued; no artifacts.

The four jobs are still waiting in GitHub Actions. There is no test PASS/FAIL evidence yet. The latest branch HEAD after duplicate cleanup is `37dec6c0f5f282c7d893f9b45b8c427b7030e121`; it is exactly the branch tip and is 14 commits ahead of the prior checkpoint `eb72a1ab...`, with no commits behind. The two L2 runs were triggered by earlier commits on the same branch and are not proof about the final tip. `main` remains untouched.

Next gate: obtain the L2 logs, verify that the three fixture-reachability tests execute (not skip), verify all three safe-domain homogeneity tests pass, and verify the three expected failures fail for the intended threshold-boundary reason. If an expectedFailure unexpectedly succeeds, remove that marker and promote it to a hard acceptance test. Do not treat a green overall unittest result as proof of production linearity while expected failures remain.


## 11. L2 fixture hardening — 11 November 2026, 01:07 UTC

A test-quality weakness was identified in the canonical L2 suite: its generic map helpers used `skipTest` when the selected reference output was empty. A suite with an empty physical fixture could therefore show skips rather than a decisive fixture failure. To prevent a vacuous green result, three ordinary hard tests were added to require nonempty output from:
- `CK.apply_volume_complete_state`;
- `CK.apply_HE_complete_state`;
- `KCOMP.apply_HE_complete_state_custom`.

The existing safe-domain and threshold-boundary tests remain unchanged. Thus the suite now has three explicit fixture-reachability checks, three hard safe-domain homogeneity checks, and three expected-failure threshold-boundary characterizations. A skip in an individual characterization can no longer conceal that the same frozen input has no production response, because the independent fixture test for that map must fail.

**Commit:** `750978aade52455db8cf680ee364a4e36493f098` — require nonzero L2 production fixtures.

**New L2 run:** `42936480504`, HEAD `750978aade52455db8cf680ee364a4e36493f098`, job `130247689781`: queued, no artifacts. Earlier L2 runs `42935840345` and `42935968925` remain queued without artifacts. Both L1 runs `42929934053` and `42930769907` remain queued without artifacts.

**Current branch tip before this checkpoint update:** `750978aade52455db8cf680ee364a4e36493f098`. `main` remains untouched. No L1/L2 test result is yet available; do not report PASS, FAIL, or a physical conclusion.


## 12. Independent-basis L2 superposition check — 11 November 2026, 01:08 UTC

The L2 suite was extended beyond one-vector rescaling. A new test constructs two distinct genuine Peter–Weyl basis states (`basis_full_jhalf()[0]` and `basis_full_jhalf()[16]`), converts both through the repository's covariant/complete adapters, verifies that the keys differ and that each has a nonzero volume response, then compares:

`CK.apply_volume_complete_state(x+y, 0)`

against

`CK.apply_volume_complete_state(x, 0) + CK.apply_volume_complete_state(y, 0)`

at a common scale (10^6), chosen to keep the test away from the declared (10^{-11}) output cutoff. The relative L2 residual is required to be at most (5\times10^{-12}). This is a genuine independent-support additivity test, not merely a homogeneity test along one input ray.

Both independent basis responses are hard assertions, so an empty second vector cannot silently turn the test into a vacuous pass. This test currently covers the volume primitive only; independent-support tests for complete H_E and full covariant C(V)/C(K) remain future work.

**Commits:** `10f9e4964e5caacbd3f2a08209fd840a29ded40f` (add independent-basis test), `ec82224564516d9b902ac52de37b3e540bce4c22` (require both basis responses).

**New L2 run:** `42936669007`, HEAD `10f9e4964e5caacbd3f2a08209fd840a29ded40f`, job `130248321716`: queued, no artifacts. The immediately prior run `42936480504` also remains queued. Older L2 and L1 runs remain queued. No PASS/FAIL is available.


Follow-up hardening of the independent-support fixture was committed as `ec82224564516d9b902ac52de37b3e540bce4c22`: the test now explicitly asserts that the two converted complete-state keys are distinct and that **both** individual volume responses are nonempty before evaluating additivity. This prevents a zero-response second basis vector from making the comparison vacuous.

The corresponding L2 run `42936684504`, HEAD `ec82224564516d9b902ac52de37b3e540bce4c22`, job `130248370659`, is queued without artifacts. Run `42936669007` (the immediately preceding version) is also queued. The later documentation-only commits do not change the test or workflow contents.



## 13. Independent-support additivity for cached complete H_E — 11 November 2026, 01:10 UTC

Added a hard L2 test to `tests/test_bqg_l2_production_primitives.py` for `CK.apply_HE_complete_state` on two independent Peter–Weyl basis inputs, using basis indices 0 and 16. The test:
- converts both inputs through `CV.gauss_to_covariant` and `COMP.covariant_to_complete`;
- asserts the complete-state keys are distinct;
- requires each basis input to have a nonzero H_E response;
- compares direct production action on the combined state with the sum of the two individual responses;
- uses a common input scale of (10^6) and relative residual bound (5\times10^{-12}).

Source-level distinction: `CK.apply_HE_complete_state` applies cached unit-key columns, combines their arbitrary input coefficients, and then applies its final (10^{-10}) output filter. This test probes the cached-column extension on a safe-scale two-vector fixture. It does not test threshold-boundary behavior or the full (C(K)), (C(V)), or K-commutator maps. The custom composition map `KCOMP.apply_HE_complete_state_custom` remains a separate and more direct coefficient-dependent threshold risk because it sends each scaled column through `COMP.add`.

The new test, its L2 audit interpretation, and README status were committed on the research branch. The scoped push workflow should run on this test change; its run identifier and outcome must be captured separately. At 01:10 UTC, all previously tracked L1/L2 jobs were still queued and had no artifacts. No PASS/FAIL is claimed. `main` remains untouched; no merge is performed.



## 14. CI re-poll and repository integrity — 11 November 2026, 01:12 UTC

The latest source change and its documentation are now present on the research branch. Repository compare against checkpoint base `eb72a1ab1b6fac36c2f2ca83d42eec83062e63ef` reports **33 commits ahead, 0 behind**. The branch tip is `14e300b19e11bc8e25802781586fb0fac79b3c9f`; comparing that commit to the branch ref reports identical. The most recent four commits are:
- `4cc8638c89cd25481c0e99027e6ee50d7a4d30b0` — independent-support additivity test for cached complete H_E;
- `4e191ff8d8f220a4c9e2572b0c31701a10a20253` — L2 audit interpretation;
- `04ee7d94f5401dd4e61d525062fad3e975155df7` — continuity checkpoint update;
- `14e300b19e11bc8e25802781586fb0fac79b3c9f` — README status update.

### Latest available Actions poll

At 01:12 UTC, the available job and artifact endpoints still reported:

| Workflow family | Run IDs | Job state | Artifacts |
|---|---|---|---|
| L2 production primitives | `42935840345`, `42935968925`, `42936480504`, `42936669007`, `42936684504` | all queued | none |
| L1 accumulator | `42929934053`, `42930769907` | both queued | none |
| V10 full-node closure | `42908538520`, `42908462865` | all six matrix jobs queued | none |
| Adapter tests | `42908462954` | queued | none |
| V11 common carrier | `42908462792` | queued | none |
| Older adapter / V11 runs | `42908538513`, `42908538772` | completed, cancelled | none |

The new test commit `4cc8638` changes a path watched by the L2 push workflow. The push-triggered run identifier was not returned by the available commit-associated run listing, which is limited to pull-request-triggered runs; do not infer a run ID or a test outcome. Its job must be identified and inspected through the Actions run listing when available.

No L1/L2/V10/V11 PASS or FAIL is claimed from queued or cancelled jobs. No artifacts are available. The added independent H_E test has not been executed in this environment. `main` remains untouched and no merge was performed.
