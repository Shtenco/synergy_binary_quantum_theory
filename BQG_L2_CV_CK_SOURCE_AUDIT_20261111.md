# BQG L2 source audit: C(V), C(K) and CI gate — 11 November 2026, 01:33 UTC

**Repository:** `Shtenco/synergy_binary_quantum_theory`  
**Branch:** `research/lemaitre-bqg-causal-cone-20261110`  
**Verified branch HEAD:** `d8a4d1b47c59e122f13eb45d43017e186c42ca15`  
**Compare base:** `eb72a1ab1b6fac36c2f2ca83d42eec83062e63ef`  
**Compare:** 47 commits ahead, 0 behind.  
**Main:** untouched; no merge.

This note records a fresh source-level inspection and a CI re-poll. It deliberately does not claim that the L2 suite executed.

## 1. Canonical test and CI evidence

Canonical test:
- `tests/test_bqg_l2_production_primitives.py`
- Git blob SHA returned by GitHub at this inspection: `7579aa87dbb5f98d9905d9259e814034b1194c7c`
- Test coverage currently present: three hard fixture-reachability checks, three hard safe-domain complex-homogeneity checks, three expected-failure threshold-boundary characterizations, and three independent-support additivity checks (local volume, cached complete H_E, one full covariant C(V) component).

The following jobs were queried at 01:33 UTC and were all still queued, with no conclusion and no artifacts:
- L2: `42937467628` (job `130250982584`), `42936684504`, `42936669007`, `42936480504`, `42935968925`, `42935840345`.
- L1: `42930769907` (job `130228642399`), `42929934053` (job `130225807106`).

For run `42937467628`, the artifact endpoint returned an empty artifact list. A queued job is neither PASS nor FAIL. No test log was available, so fixture reachability, safe-domain homogeneity, XFAIL reasons, and independent-support residuals remain unverified by execution.

## 2. Exact C(V) composition structure

Source: `scripts/peter_weyl_covariant_composition_gate.py`.

The inspected implementation has these stages:

1. `direct_volume_covariant` converts the covariant state to the complete representation, applies `CK.apply_volume_complete_state`, and converts back.
2. `inverse_complete` applies the inverse holonomy branch and complete-charge projection.
3. `close_complete` applies the forward holonomy branch and covariant projection.
4. `C_volume_component` forms the hVh contribution by summing over the auxiliary index `k in range(2)`, applying the volume map between inverse and forward holonomy.
5. The direct term is added only when `i == j`; the hVh contribution is subtracted for every matrix component.

The accumulator `COMP.add` uses `TOL = 1e-11` and drops/removes running partial sums below the threshold. It is called inside inverse projection, closure, hVh accumulation, and final matrix-component assembly. Therefore the complete C(V) map still contains state-dependent threshold boundaries even though one safe-scale independent-support test is present.

**Scope of the current full C(V) test:** edge (0,1), matrix component (i,j)=(0,0), two independent Gauss Peter–Weyl basis inputs (indices 0 and 16), common scale 10^6, relative L2 residual ceiling 5e-12. It does not cover the other three matrix components. Since the direct term is structurally conditional on i=j, off-diagonal components must be tested as their own maps; they must not be inferred from the (0,0) result.

## 3. Exact C(K) composition structure and risk surface

Source: `scripts/peter_weyl_covariant_K_composition_gate.py`.

The current composition chain is:

- `apply_HE_complete_state_custom`: for every complete basis key, obtains a unit-key `CK.apply_HE_complete_key` result and sends that result, multiplied by the actual input coefficient `amp`, through thresholded `COMP.add`.
- `apply_K_complete_custom`: constructs `V H_E` and `H_E V) legs separately and combines them with `COMP.add) using opposite signs.
- `direct_K_covariant`: maps covariant input to complete form, applies the custom K map, then maps back.
- `C_K_component`: constructs the inverse-holonomy branch, applies custom K while source and target nodes are charged, closes the target, and combines the hKh term with the diagonal direct K term using thresholded `COMP.add`.

This is not merely a fixed matrix sparsification. The arbitrary input coefficient enters before a thresholded accumulation in `apply_HE_complete_state_custom`; the subsequent K commutator and C(K) composition introduce additional thresholded combination stages. Hence the same linearity question must be asked independently at each stage.

No C(K) additivity or L3 result is claimed by this source inspection.

## 4. Why the next test sequence is gated

Do not expand the matrix-component matrix or introduce C(K) tests before the existing full C(V) run is diagnosed. The correct sequence is:

1. Obtain run `42937467628` job log/artifact.
2. Confirm the three hard fixture tests execute rather than error; inspect all three safe-domain homogeneity outcomes.
3. Inspect each XFAIL and verify that the failure is specifically the intended threshold-boundary additivity assertion, not an unrelated failure.
4. Inspect the independent-support volume, cached H_E and full C(V) residuals.
5. If C(V) passes, extend the same frozen fixture to the remaining three matrix components, with hard nonzero checks where the exact map is expected to have nonzero support. If a component is identically zero for the chosen fixture, replace or supplement the fixture rather than treating a vacuous zero=zero comparison as evidence.
6. Only after C(V) is diagnosed, add a separate C(K) independent-support test for `COMP.C_K_component`, initially on a declared edge/component and with both individual responses required nonzero.
7. Then perform L3 reconstruction: freeze a finite input basis, compute each production basis column, form the matrix action on deterministic complex coefficient vectors, and compare it with direct production action on the same superpositions. Include phase changes, cancellation and near-threshold vectors.

Do not relax tolerances merely to obtain a green test. If a test fails, classify the failure by the actual threshold/caching/projection stage and preserve the failing residual and support diagnostics.

## 5. Scientific interpretation

The present state is a stronger and more explicit executable characterization plan, not a proof that the production KKV map is linear. A basis-column matrix can represent a linear extension even when the actual thresholded function on arbitrary superpositions is nonlinear. Consequently:

- a passing C(V) single-component test is local evidence for that finite fixture only;
- it does not establish all four components, C(K), the K commutator, or the full Lorentzian constraint;
- matrix-level adjoint residuals remain downstream diagnostics until production linearity is established on the declared domain;
- no physical normalization, physical kernel, HDA closure, or causal cone follows from this L2 test family.

The strict sparse accumulator remains a reference implementation and is not yet wired into the active V10 KKV path.

## 6. Record of this inspection

At 01:33 UTC the branch was unchanged from `d8a4d1b47c59e122f13eb45d43017e186c42ca15`. This note adds source-contract detail and an explicit sequencing gate only. No production operator, test, workflow, numerical threshold, or physical claim was changed. No local test execution was performed. `main` remains untouched.


## 7. L2 workflow dependency-closure correction — 11 November 2026, 01:38 UTC

A CI wiring defect was found by comparing the test's import graph with the workflow's `paths` and `py_compile` lists. The L2 test calls `COMP.C_volume_component`, but the workflow originally did not watch or explicitly compile `scripts/peter_weyl_covariant_composition_gate.py`. Thus a change to that production module could fail to trigger the L2 workflow, and its syntax was not covered by the dedicated compile step.

The workflow has now been corrected and the correction itself reviewed for duplicate entries. It explicitly watches and compiles the L2 import/dependency closure:
- `peter_weyl_covariant_composition_gate.py`;
- `k5_peter_weyl_safe_hda_column.py`;
- `charged_intertwiner_recoupling_gate.py`;
- `peter_weyl_lorentzian_K_block_gate.py`;
- `peter_weyl_covariant_K_projection_audit_gate.py`;
- `peter_weyl_volume_cross_rep_audit.py`;
- `peter_weyl_q_block_gauge_audit.py`;
alongside the already tracked L2 test, CK, KCOMP, CV and zero-aware volume migration modules.

This is a CI coverage correction, not a production-operator change and not a test result. The push should create a new workflow execution; its actual run identifier and status must be obtained and inspected. No PASS is inferred from the workflow file being updated. The earlier run `42937467628` remains separately identified as the test-bearing run for the C(V) gate and must still be inspected. Do not confuse a later workflow-configuration run with execution of a different test revision.


## 8. Re-audit of the actual CI queue and C(K) coefficient path — 11 November 2026, 01:57 UTC

### 8.1 Direct CI re-poll

The eight previously identified L1/L2 runs were queried again through both the job and artifact endpoints:

| Family | Run | Job | Job status | Conclusion | Artifacts |
|---|---:|---:|---|---|---|
| L2 | 42937467628 | 130250982584 | queued | none | none |
| L2 | 42936684504 | 130248370659 | queued | none | none |
| L2 | 42936669007 | 130248321716 | queued | none | none |
| L2 | 42936480504 | 130247689781 | queued | none | none |
| L2 | 42935968925 | 130246013044 | queued | none | none |
| L2 | 42935840345 | 130245591568 | queued | none | none |
| L1 | 42930769907 | 130228642399 | queued | none | none |
| L1 | 42929934053 | 130225807106 | queued | none | none |

This repeated observation establishes a queue/observability blockage for these tracked runs, not a numerical result. In particular, run 42937467628 still has no execution log or artifact. No PASS, FAIL, skip count, XFAIL cause, or numerical residual can be inferred. The connector still does not expose an all-push-run listing or workflow-dispatch operation, so the run created by the later workflow correction cannot be identified from the available interface. No run ID is guessed.

### 8.2 Additional C(K) source-level precision

The current `peter_weyl_covariant_K_composition_gate.py` source confirms three distinct coefficient-sensitive layers:

1. **Column construction:** `apply_HE_complete_state_custom` obtains a unit-key response from `CK.apply_HE_complete_key`, then passes that response and the actual arbitrary input coefficient `amp` to `COMP.add`. Since `COMP.add` thresholds the running destination after scaling, the retained support can depend on the magnitude and ordering of coefficients in a superposition.
2. **Commutator assembly:** `apply_K_complete_custom` forms `VH=V(HE(state))` and `HV=HE(V(state))` separately, then combines them with opposite signs using the same thresholded accumulator. Cancellation between these two terms is therefore also exposed to intermediate pruning.
3. **Outer covariant leg:** `C_K_component` sums the two auxiliary-index branches, closes the target endpoint, and combines the hKh term with the direct diagonal term through further thresholded `COMP.add` calls.

This distinguishes the custom state action from the cached-column `CK.apply_HE_complete_state`: the latter assembles cached unit-key columns by ordinary coefficient multiplication before its final output filter, whereas the custom map explicitly inserts each scaled column into the thresholded accumulator. Neither implementation should be presumed linear solely because its unit-key columns can be serialized as a matrix.

### 8.3 L2 acceptance contract sharpened

When Actions becomes observable, the report must retain separate outcomes for:
- syntax/compile closure;
- hard fixture reachability;
- hard safe-domain complex homogeneity;
- expected threshold-boundary failures, including their exact assertion messages;
- each independent-support additivity residual and output support size.

A red suite caused by an expected threshold-boundary characterization is not equivalent to a broken fixture; a green suite with XFAILs is not proof of linearity. Any unexpected success of an `expectedFailure` is a test-contract change requiring review, not something to ignore.

### 8.4 Scientific boundary

No production operator or test was changed in this re-audit. No local execution was performed. The active C(V) test still covers only one matrix component, and no C(K) additivity or L3 matrix-reconstruction result exists. The next mathematical gate remains the actual execution/diagnosis of the current C(V) fixture; the next software gate is to restore CI observability before expanding claims. `main` remains untouched.
