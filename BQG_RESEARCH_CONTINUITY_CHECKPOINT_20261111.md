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
- The same module's `apply_K_complete_custom` combines `V H_E) and `H_E V` components through further thresholded accumulation.
- `peter_weyl_covariant_K_leg_gate.py`: `add_state` thresholds branch sums; `apply_HE_complete_key` uses it while assembling oriented/adjoint branches; `apply_volume_complete_state` applies a final (10^{-11}) output filter; `apply_HE_complete_state` applies a final (10^{-10}) output filter.
- `peter_weyl_covariant_K_composition_gate.py` therefore contains distinct loss boundaries at primitive branch assembly, arbitrary coefficient composition, volume action, and commutator combination. They must be measured independently.
- V10's outer `amp_tol) pruning occurs after coherent node-map summation, but the removed forward/adjoint amplitudes are not currently counted. This external loss must not be merged with internal primitive losses or scalar-projection rejection.

The new L2 test plan prioritizes real production maps: `CK.apply_volume_complete_state`, `CK.apply_HE_complete_state`, `COMP.apply_HE_complete_state_custom`, `COMP.apply_K_complete_custom`, `COMP.C_K_component`, and `COMP.C_volume_component`. It requires additivity, complex homogeneity, permutation/cancellation tests, direct-superposition versus frozen basis-matrix reconstruction, and stage-specific loss ledgers.

**New commit:** `4f888c8658299af441097674fa08a6949a837be9` — add L2 production-primitive linearity audit.

**CI status remains unchanged:** runs `42930769907` and `42929934053` were observed queued, with no artifacts. The existing L1 workflow does not execute the new L2 test matrix. No L1/L2/L3 PASS is claimed.

**Next gate:** implement and execute the three highest-priority L2 tests before changing production semantics. Do not wire the reference accumulator into V10 until the exact insertion boundary and regression behavior are established.
