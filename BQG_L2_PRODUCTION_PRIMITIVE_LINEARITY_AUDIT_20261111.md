# BQG L2 production-primitive linearity audit — 11 November 2026

**Repository:** `Shtenco/synergy_binary_quantum_theory`  
**Branch:** `research/lemaitre-bqg-causal-cone-20261110`  
**Inspected reference HEAD:** `eb72a1ab1b6fac36c2f2ca83d42eec83062e63ef`  
**Purpose:** extend the L1 source audit into an exact map-by-map L2 test plan, using the active V10 KKV production call graph. This is a source-level audit and preregistration document; it does not claim that L2 tests have run or passed.

## 1. New source-confirmed findings

The active implementation contains more than the two already documented dictionary accumulators. Thresholding occurs at several distinct mathematical layers. They must be audited separately rather than grouped under one generic "pruning" label.

### 1.1. Composition accumulator: arbitrary input coefficients are thresholded

In `scripts/peter_weyl_covariant_K_composition_gate.py`:

- `apply_HE_complete_state_custom` obtains a unit-column-like `part` from `CK.apply_HE_complete_key(key,...)`, then calls `COMP.add(out, part, amp)` (lines 43–51).
- `amp` is the coefficient of an arbitrary input superposition. `COMP.add` compares the running sum `dst[key] + amp*part[key]` with `TOL=1e-11`, retaining or deleting the key.
- Consequently, the operation applied to an arbitrary state is not guaranteed to equal the linear extension of its basis columns. This is a direct L2/L3 target, not merely a hypothetical risk.

The same module's `apply_K_complete_custom` uses this path for both `V H_E` and `H_E V`, then combines the two maps through further `COMP.add` calls. The component-wise commutator assembly therefore has another threshold boundary.

### 1.2. K-leg branch accumulator and projection

In `scripts/peter_weyl_covariant_K_leg_gate.py`:

- `add_state` thresholds the running complex sum at `1e-11`.
- `apply_HE_complete_key` invokes `add_state` while accumulating orientation, adjoint and primitive-sequence branches.
- `project_branch_complete_charges` applies coefficient cutoffs while resolving local charge sectors and also filters the projected output by `tol`.
- `apply_volume_complete_state` filters its output amplitudes at `1e-11`.
- `apply_HE_complete_state` accumulates cached basis-column outputs, then filters the final state at `1e-10`.

Interpretation must be staged: a cutoff used to define a fixed basis-column approximation is not, by itself, proof that the assembled matrix action is nonlinear. A cutoff applied after arbitrary input coefficients have been mixed, or to a state-dependent output, is nonlinear in general. Each boundary therefore needs its own controlled test and loss counter.

### 1.3. V10 closure has an additional external threshold

In `scripts/bqg_v10_w5_lorentzian_adjoint_closure.py`, after node maps have been coherently accumulated into `f_total` and `ad_total`, the closure executor retains only amplitudes satisfying `abs(z) > amp_tol`. This is a declared support approximation and is distinct from the internal thresholds above.

The current closure diagnostics count scalar-projection rejections and support keys omitted because of `max_states`/`max_depth`, but do not record the forward/adjoint amplitudes removed by `amp_tol`. A complete loss ledger must measure this external loss and internal losses independently.

## 2. L2 test matrix (required before production refactor)

Run against the real implementations, not synthetic stand-ins. Freeze the source node, target node, finite spin wall, zero-aware volume convention, cache state and all thresholds in each test record.

| Target | Required checks | Failure mode isolated |
|---|---|---|
| `COMP.add` | coherent subthreshold sum; cancellation/re-entry; all permutations | order-dependent running-sum deletion |
| `CK.add_state` | same as above, complex phases | branch accumulator nonlinearity |
| `CK.apply_HE_complete_key` | branch-order permutations; coefficient scaling at the branch boundary; compare with no-prune reference | primitive branch loss |
| `CK.apply_volume_complete_state` | homogeneity/additivity around (10^{-11}); output loss norm | state-output pruning after fixed volume action |
| `CK.apply_HE_complete_state` | additivity/homogeneity at (10^{-10}); compare direct superposition to cached basis columns | post-assembly state pruning |
| `COMP.apply_K_complete_custom` | additivity/homogeneity and separate VH/HV loss accounting | thresholded commutator assembly |
| `COMP.C_K_component` | direct action on (x+y) vs sum of actions; coefficient scaling; edge-index permutations | complete covariant K-leg map |
| `COMP.C_volume_component` | same checks with the exact V10 zero-aware environment | complete covariant volume-leg map |

For each map F, report both absolute and relative residuals:

```text
r_add_abs = ||F(x+y) - F(x) - F(y)||_2
r_hom_abs = ||F(alpha*x) - alpha*F(x)||_2

r_add_rel = r_add_abs /
  max(||F(x+y)||_2, ||F(x)||_2 + ||F(y)||_2, epsilon)

r_hom_rel = r_hom_abs /
  max(||F(alpha*x)||_2, |alpha|*||F(x)||_2, epsilon)
```

Use deterministic complex vectors with (a) amplitudes well above all thresholds, (b) amplitudes straddling each threshold, (c) coherent sums of individually subthreshold terms, and (d) cancellation followed by re-entry. Include several non-real (alpha), not just positive real scalars. Preregister numerical tolerances from precision and conditioning before inspecting the results.

## 3. L3 whole-map reconstruction

On one frozen finite input basis (e_j), define columns (M_j=F(e_j)) using the exact same cache and zero-aware environment. For deterministic complex coefficient vectors (c), independently evaluate the production function (F(sum_j c_j e_j)) and compare it to (Mc=sum_j c_jM_j).

Report:
- input/output basis hashes and dimensions;
- all fixed matrix masks and their thresholds;
- direct production output support and matrix-reconstructed support;
- absolute and relative reconstruction residual;
- additivity and homogeneity residuals;
- internal stage losses and final closure `amp_tol` losses separately;
- deterministic seed and coefficient-vector hash.

A matrix built from basis columns is not an adequate substitute for this test: the direct production call on a nontrivial superposition is the independent side of the comparison.

## 4. Loss accounting contract

At every threshold boundary, before deleting an entry, accumulate:

- discarded term count (not only unique keys);
- discarded complex L2 norm squared, with a clear convention: per discarded contribution or per coherently combined key;
- maximum discarded amplitude;
- threshold value and stage identifier.

For the closure-level coherent maps, count below-threshold forward and adjoint keys after coherent node summation, with L2 norm and maximum amplitude. Do not combine these with scalar-projection leakage, basis-reconstruction leakage, state/depth caps or internal branch losses. They represent different failure channels.

Avoid double counting: the ledger should identify whether a loss is measured before or after coherent combination. For each stage, document whether its counter is contribution-wise or key-wise.

## 5. Strict-linear architecture decision

The safest near-term design is not to replace every helper blindly. Preserve the current finite-sector construction while separating two modes:

1. **Column construction:** compute each basis column with explicit representation/spin walls and separately reported projection leakage. Any fixed coefficient mask must be frozen and hashed.
2. **Operator application:** assemble the frozen sparse matrix once, then apply it to arbitrary coefficient vectors with the strict accumulator and no state-amplitude pruning.

The direct reversed-order adjoint should then be evaluated on the same frozen domain and compared to the adjoint of the same frozen matrix. If the existing primitive pipeline is retained for independent direct-action comparison, its nonlinear thresholds must remain visible and must not be conflated with the strict matrix action.

Do not wire the L1 reference accumulator into the live V10 path until L2/L3 establish the exact insertion boundary and the regression suite is executable. In particular, do not silently change thresholds or projection semantics while trying to make linearity tests pass.

## 6. CI and evidence boundary

At the time of this audit, L1 runs `42930769907` and `42929934053` were still queued and had no artifacts. The existing L1 workflow runs reference-accumulator tests and imports/compiles selected production modules; it does not execute the L2 map-by-map matrix above. No L1 PASS, L2 PASS, production-linearity PASS, spectral convergence or physical conclusion is claimed here.

**Next gate:** implement L2 tests first for `CK.apply_volume_complete_state`, `CK.apply_HE_complete_state`, and `COMP.apply_HE_complete_state_custom`, because these isolate the state-level thresholds that sit directly on the active KKV path. Then extend to `C_K_component`, `C_volume_component`, and finally the L3 whole-map reconstruction.


## 7. L2 executable characterization added — 11 November 2026, 01:00 UTC

Two executable artifacts were added after the source audit:

- `tests/test_bqg_l2_production_primitives.py`: constructs a real finite Peter–Weyl basis state, converts it through the existing covariant/complete representation adapters, applies the same zero-aware volume convention as the research composition path, and calls the actual production functions:
  - `CK.apply_volume_complete_state`;
  - `CK.apply_HE_complete_state`;
  - `KCOMP.apply_HE_complete_state_custom`.
- `.github/workflows/bqg-production-linearity-l2.yml`: compiles the relevant modules, runs the L2 suite on Python 3.11, and preserves the execution log as an artifact.

The tests have two roles:

1. Hard safe-domain complex-homogeneity checks at a scale well above the output thresholds.
2. `unittest.expectedFailure` threshold-boundary characterizations. For each map, the test calibrates an input coefficient from the actual unit-state output so that two individually subthreshold contributions are compared with their coherent sum. These are intentionally expected to expose the currently documented output/accumulator threshold behavior.

The fixture is a real finite-sector key, but the tests are still finite-sector numerical characterization, not a proof of full operator linearity, convergence, or physical correctness. The L2 workflow has just been added; no run result is claimed until its job and artifact are inspected. The tests do not yet exercise `COMP.C_K_component` or `COMP.C_volume_component`; those remain subsequent integration-level tests.

**Next gate:** inspect the first L2 workflow execution. If safe-domain homogeneity fails, investigate cache/patch contamination or a genuine fixed-domain mismatch before relaxing tolerances. Confirm that each expected failure occurs for the intended threshold reason rather than a fixture/import exception.


## 7. Executable L2 characterization now present

The canonical executable suite is `tests/test_bqg_l2_production_primitives.py`; its runner is `.github/workflows/bqg-production-linearity-l2.yml`. It uses a real finite Peter–Weyl basis vector, maps it through the covariant/complete representation adapters, applies the same zero-aware volume convention used by the V10 composition experiment, and clears dependent caches.

Three production maps are exercised:
- `CK.apply_volume_complete_state`;
- `CK.apply_HE_complete_state`;
- `KCOMP.apply_HE_complete_state_custom`.

For each map, a safe-domain complex-homogeneity check is a hard assertion. A separate expected-failure test calibrates the input scale from the actual unit-column output so that each of two separate contributions falls below the output cutoff while their coherent sum crosses it. These are explicit characterizations of the current state-dependent thresholds, not acceptable operator behavior.

The suite also checks the fixture is a single valid complete covariant basis key. Empty map outputs cause the safe-domain test to skip and are not to be interpreted as a pass; the CI log must be inspected for skips and fixture checks. The workflow preserves the complete unittest log as an artifact with `if: always()`.

At the 01:05 UTC poll, L2 runs `42935840345` and `42935968925` were both queued and had no artifacts. The two L1 runs `42929934053` and `42930769907` were also queued without artifacts. No execution outcome is claimed. The L2 workflow is a first primitive gate only: it does not yet exercise `COMP.apply_K_complete_custom`, `COMP.C_K_component`, `COMP.C_volume_component`, or the L3 arbitrary-superposition matrix reconstruction.

One temporary duplicate L2 test file was removed in commit `37dec6c0f5f282c7d893f9b45b8c427b7030e121`; the canonical suite and workflow above remain the single source of truth.


## 8. Fixture reachability is a hard gate

The canonical L2 test suite was hardened in commit `750978aade52455db8cf680ee364a4e36493f098`. It now has separate nonempty-output assertions for the frozen Peter–Weyl input under all three tested maps: volume, cached complete (H_E), and custom complete (H_E) composition. This prevents the helper's empty-output `skipTest` behavior from masking a fixture that never reaches the operator.

The current suite therefore separates three distinct outcomes:
- fixture reachability: hard PASS/FAIL;
- safe-domain complex homogeneity: hard PASS/FAIL (unless the map helper itself reports empty, which must be read alongside the hard fixture test);
- threshold-boundary nonlinearity: expected failure while the current state-amplitude cutoffs remain.

The newly triggered run `42936480504` (HEAD `750978aade52455db8cf680ee364a4e36493f098`, job `130247689781`) was queued at 01:07 UTC with no artifacts. Earlier L2 runs `42935840345`, `42935968925`, and L1 runs `42929934053`, `42930769907` were also queued. No test verdict is available.


## 9. Independent-support additivity for the volume primitive

The canonical test file now also checks additivity on two distinct Peter–Weyl basis keys, rather than relying only on rescaling one input vector. It uses `basis_full_jhalf()[0]` and `basis_full_jhalf()[16]`, verifies distinct complete covariant keys and nonzero individual volume responses, then compares the direct production action on their sum with the sum of separate actions.

The common input scale is (10^6) and the relative residual bound is (5\times10^{-12}), intended to keep the comparison safely above the (10^{-11}) output-amplitude cutoff. This extends L2 from one-dimensional homogeneity/threshold characterization to independent-support additivity for `CK.apply_volume_complete_state`. It does not yet establish additivity of complete (H_E), the K commutator, or the full covariant legs.

The scoped L2 workflow was triggered as run `42936669007` (job `130248321716`), queued with no artifacts at 01:08 UTC. No result is claimed.
