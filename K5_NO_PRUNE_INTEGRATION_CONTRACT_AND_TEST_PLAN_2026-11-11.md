# K5 no-prune integration contract and numerical acceptance plan — 2026-11-11

**Repository:** `Shtenco/synergy_binary_quantum_theory`  
**Research branch:** `research/consolidated-bqg-audit-20261111`  
**Inspected producer:** `scripts/k5_peter_weyl_safe_hda_column.py`  
**Inspected producer blob:** `e7e1e0705c4b833dcee167a2036ef5f8929a4388`  
**Purpose:** specify a safe, testable route from the standalone coherent-pruning primitive to the actual K5 Peter–Weyl HH producer. This is an implementation contract and source audit, not a claim that the producer has already been refactored or that a no-prune HH run has passed.

## 1. Executive finding

The K5 producer has a chain of amplitude cutoffs at mathematically different locations. Removing only the visible dictionary threshold is insufficient. In particular, local intertwiner projection, primitive path output, cached T output, H/HH composition, and final commutator construction all have separate loss points.

The producer therefore has three distinct tasks:

1. **Preserve the historical result:** retain a clearly named legacy path with its exact old thresholds and regression constants.
2. **Create a strict reference path:** thread an explicit no-prune policy through every amplitude-producing and caching layer; coherently aggregate all contributions by canonical output key before any optional threshold.
3. **Establish evidence:** compare strict sparse action with a tiny dense oracle, test adjoints bilinearly, inspect the raw volume operator before symmetrization, and only then run larger HH/shell computations.

No default should silently change during this migration. No no-prune call may return a legacy-thresholded cache entry.

## 2. Exact cutoff map in the inspected source

The following are source-level observations from blob `e7e1e0705c4b833dcee167a2036ef5f8929a4388`; line references refer to that exact blob.

| Location | Current operation | Why it matters |
|---|---|---|
| `add_dict`, lines 191–195 | Running sum is retained only if `abs(v)>1e-13`; otherwise an existing key is deleted | Partial cancellation can delete a key and later recreate it; result depends on contribution order |
| local channel extraction, lines 249–256 | Local intertwiner coefficient is ignored unless `abs(c)>1e-12` | A physically nonzero local channel can disappear before global output assembly |
| primitive path assembly, lines 258–264 | Product amplitude is filtered at `1e-12`, then final output map at `1e-11` | Contributions can be removed before they are coherently combined with other paths |
| T cache, lines 273–282 | Cached T output is filtered at `1e-10` | Cached values are already approximate; downstream code cannot recover removed support |
| state/HH path, lines 284–293 and 326–331 | `prune_state` defaults to `1e-8`; HH composition applies it | The two-hit state is threshold-dependent |
| top-level H and commutator, lines 333–342 | H0/H1 and final commutator are pruned at `1e-8` | Reported support and norms are not complete unpruned operator action |

The separate graph-shell wrapper was previously source-audited as adding another `1e-10` threshold. That wrapper must be included if it is used in the final shell pipeline.

**Important:** the cutoffs above are not one global tolerance. A sum of their numerical values, a count of dropped keys, or a single final norm is not a valid propagated error estimate. Each loss must be measured at its own stage and propagated through downstream operator norms.

## 3. Required policy architecture

Use an immutable policy object with at least:

- `mode`: `legacy`, `no_prune`, or `aggregate_then_prune`;
- `thresholds`: named stage-to-threshold mapping, not one ambiguous global scalar;
- `policy_version`: stable identifier included in run metadata;
- `ledger`: stage-specific discarded-vector records.

The standalone helper currently supports only `no_prune` and `aggregate_then_prune`. Keep legacy behavior in the producer adapter rather than pretending that post-aggregation pruning reproduces historical intermediate thresholding.

### 3.1. Mode semantics

**legacy**
- Reproduce the existing implementation exactly, including each intermediate cutoff and order of operations.
- Preserve the historical regression values.
- Label all resulting output as legacy-thresholded, not pruning-certified.

**no_prune**
- No amplitude threshold at any stage.
- Exact mathematical zeros may be omitted; no nonzero value may be removed by a tolerance.
- Coherently aggregate all terms sharing a canonical output key.
- If resource limits are exceeded, terminate explicitly with a resource-limit status. Never silently fall back to legacy.

**aggregate_then_prune**
- First complete coherent aggregation for the named stage.
- Apply only the explicitly declared stage threshold after aggregation.
- Record the exact discarded vector, its L1/L2 norms, maximum amplitude, key count and digest.
- Never use this mode to claim equivalence with the old legacy path.

## 4. Threading requirements and cache safety

The policy must be an explicit argument through the full dependency chain:

1. `apply_primitive_oriented`
2. `apply_T_oriented_basis`
3. `T_cached`
4. `apply_T_cached_state`
5. `apply_H_cached_state`
6. `prune_state` / state composition
7. `compose_on_sparse`
8. top-level H0/H1/commutator assembly
9. graph-shell column generation and leakage construction.

Every cache whose output can depend on a policy or threshold must include the complete policy identity in its cache key. In particular, a strict request must never reuse a legacy T result. A safe migration can either:
- add a frozen, hashable policy key to every relevant `lru_cache` signature; or
- keep separate cache namespaces/functions for legacy and strict modes.

The second option is easier to audit initially. Cache identity must include all numerical conventions that affect the value: spin wall, orientation convention, policy version, stage thresholds and relevant volume convention.

## 5. Coherent aggregation rules

The existing helper `scripts/bqg_sparse_pruning.py` is the reference for the aggregation primitive. It is deliberately separate from the K5 producer today.

For each output key (k), form
[
a_k=\sum_{r:\,k_r=k} a_r
]
before deciding whether to keep or discard it. Use deterministic ordering and stable real/imaginary summation. Do not delete a partial sum because it temporarily falls below a threshold.

Two specific source locations need careful treatment:

- `add_dict`: the legacy running-sum deletion must remain untouched in legacy mode; strict mode should collect contributions and aggregate at the end of the stage.
- `apply_primitive_oriented`: local channel coefficients and products must not be tolerance-filtered in no-prune mode. The resulting full canonical output keys must be coherently combined before any stage-level pruning.

Do not confuse a small local coefficient with an exact zero. Exact zeros can be dropped; a numerical tolerance is a policy decision and must be ledgered.

## 6. Numerical acceptance tests

### Gate A — primitive aggregation

Use deterministic synthetic contributions and actual K5-produced branches.

- Permute the same branch list and require the same canonical sparse map within a predeclared floating-point bound.
- Constructive accumulation: several individually subthreshold contributions to one key must survive when their coherent sum exceeds the threshold.
- Cancellation: large opposite contributions plus a small residual must retain the residual if it is nonzero.
- Complex phase: repeat with non-real amplitudes and cancellation in both real and imaginary parts.
- Verify no-prune retains every finite nonzero canonical amplitude.

### Gate B — tiny dense oracle

Choose the smallest nontrivial spin wall and a reduced set of links/legs for which the complete operator can be explicitly materialized as a dense matrix. Compute:
- each basis column using the strict sparse path;
- the dense matrix action on deterministic complex vectors;
- the strict sparse action directly on those same vectors.

Require agreement of the two paths in both complex amplitudes and L2 norm. Record basis ordering, input/output dimensions, source hash, policy hash and exact spin wall. A synthetic dictionary oracle alone does not satisfy this gate; at least one oracle must use the actual K5 branch/recoupling routines.

### Gate C — adjointness

The reverse-sequence construction is structurally plausible but is not independently certified by the source inspection. On a frozen finite domain, choose independent deterministic complex vectors (x,y) and test
[
\langle x,Ty\rangle=\langle T^\dagger x,y\rangle.
]
Test each primitive path, each T term, the assembled node operator, and the Hermitian completion separately. Include boundary states at the spin wall; the test domain must make clear whether the adjoint is the finite-wall adjoint or the restriction of a larger-wall operator.

Do not infer adjointness merely because the final matrix was symmetrized.

### Gate D — volume integrity

In `volume123_matrix`, the current code constructs (Q) and immediately replaces it by ((Q+Q^\dagger)/2) before diagonalization. This prevents the reported result from revealing the raw assembly defect.

The corrected diagnostic must:
1. retain `Q_raw`;
2. measure (\|Q_{raw}-Q_{raw}^\dagger\|/\max(\|Q_{raw}\|,\epsilon));
3. record the Hermitian and anti-Hermitian parts separately;
4. reject or explicitly justify any non-negligible anti-Hermitian component;
5. only then form the volume operator from the justified Hermitian matrix.

Also test volume positivity, spectral reconstruction residual, and covariance under the orientation/epsilon transforms used by the producer.

### Gate E — threshold ladder and error propagation

After strict mode passes the small oracle, run a named threshold ladder such as (10^{-8},10^{-10},10^{-12},10^{-14}), plus no-prune where feasible. Compare on a common embedded basis:
- support changes;
- output L2 difference;
- commutator norm and fixed/outside-shell fractions;
- eigenvalue and eigenspace changes for the relevant Gram/master matrices.

For stage errors (\delta_s), report the bound
[
\|x_n-\widetilde x_n\|_2
\leq\sum_s
\|A_{n-1}\cdots A_{s+1}\|_2\,\|\delta_s\|_2.
]
If downstream operator norm products are not known or bounded, report the stage losses without claiming a global error certificate.

## 7. Graph-shell acceptance

The no-prune K5 column is a prerequisite for the coherent leakage-Gram construction, not a replacement for it. For shell projection (P_S), define
[
C_v=P_SH_vP_S,\qquad L_v=(I-P_S)H_vP_S.
]
With complete output columns and the correct orthogonal/Gram-weighted inner product,
[
P_SH_v^\dagger H_vP_S=C_v^\dagger C_v+L_v^\dagger L_v.
]
Build the full matrix (G_v^{leak}=L_v^\dagger L_v), matching outside-shell canonical keys across different input columns. A scalar outside norm per column is not enough because it discards off-diagonal interference.

The compressed-only master and the pullback-corrected master must be reported separately. Their kernels are not interchangeable.

## 8. Cross-branch continuity: V10 production-linearity audit

At the time of this update, `main` has two commits not contained in this research branch:
- `8d55ff41aaad648fdde92986e98fa9c5a8f8db54`: adds `BQG_L2_PRODUCTION_PRIMITIVE_LINEARITY_AUDIT_20261111.md`;
- `254f19fdab41d32b920bf454f987aa7c0fdba34d`: adds `BQG_RESEARCH_CONTINUITY_CHECKPOINT_20261111.md` and records the L2 audit extension.

Those documents concern the active V10 KKV production call graph in different scripts. Their source-confirmed findings must not be automatically attributed to `k5_peter_weyl_safe_hda_column.py`: the continuity checkpoint explicitly distinguishes adjacent legacy helpers unless a direct caller is confirmed. The useful shared methodology is the L1/L2/L3 hierarchy:
- L1: accumulator and helper-level characterization;
- L2: actual production primitive additivity/homogeneity;
- L3: direct arbitrary-superposition action versus frozen basis-matrix reconstruction.

The K5 producer requires its own corresponding gates above. This branch is intentionally not merged with `main` in this update; its divergence is recorded rather than hidden. No history rewrite or merge is performed here.

## 9. Current status and acceptance language

**Completed**
- Source inspection of the exact K5 producer blob and its nested thresholds.
- Independent deterministic coherent-pruning helper and eight committed tests.
- Design of stage-specific ledgers and a telescoping error bound.
- Identification of the pre-symmetrization volume diagnostic gap.
- Exact leakage-Gram identity and the requirement to preserve off-diagonal outside-shell interference.

**Not completed**
- Policy threading into the K5 producer.
- Direct execution of the committed pruning unit-test file in a local checkout.
- Actual K5 tiny dense-oracle test.
- Bilinear adjoint tests on K5 production primitives.
- Raw volume anti-Hermitian measurement.
- No-prune HH run, threshold ladder, or pruning-certified leakage Gram.

**Release rule:** do not describe the current K5 HH output as unpruned, pruning-certified, adjoint-certified, or physically HDA-closed. It remains a legacy-thresholded finite regression until the corresponding gates above pass with reproducible run manifests.


## 10. Producer integration update — 2026-11-11

The K5 producer now has an explicit CLI mode switch:

- --mode legacy (default): retains the historical threshold path and regression comparison;
- --mode no_prune: bypasses local intertwiner, primitive-output, cached T, state, HH-composition, top-level H/commutator, and V5-input amplitude cutoffs. It coherently aggregates primitive path contributions, T-sequence maps, state superpositions, node-H contributions, HH compositions, the commutator, and the V5 superposition using the shared deterministic sparse aggregator.

The T_cached cache signature includes mode, so legacy and no-prune values occupy distinct cache keys. A MemoryError in no-prune mode emits an explicit resource_limit JSON status and exits nonzero; it does not fall back to legacy. The volume routine now records the raw pre-symmetrization relative anti-Hermitian defect for each encountered spin triple before retaining the prior Hermitian symmetrization used by the historical calculation.

Added scripts/test_k5_sparse_policy_wiring.py with six focused checks: tiny-amplitude retention in no-prune accumulation, legacy cutoff preservation, no-prune state retention, legacy state pruning, policy-bearing cache signature, and raw volume defect recording. CI installs NumPy and SymPy and runs both pruning and K5 wiring test scripts on push and pull-request path changes.

Not yet certified: the full no-prune K5 HH run has not been executed; the new focused test scripts have not yet returned a confirmed CI PASS; no actual K5 dense-matrix oracle or production bilinear adjoint test has been completed. The V5 input coefficient selection is now threshold-free in no-prune mode, but the no-prune result still needs a real run manifest and independent comparison. The raw volume defect is exposed diagnostically, not yet gated against a declared tolerance or independently derived exact value.


## 11. Exact implementation snapshot and source-level threshold review

Latest K5 producer blob after integration: a87cc00c3be07b4a66c249c16ec840cb8cb5629d (17378 characters). The source-level review confirms that each remaining absolute-amplitude cutoff is guarded by the legacy branch or occurs after an early return from the no-prune branch:
- add_dict: the 1e-13 update cutoff is only applied when mode is not no_prune;
- primitive local channel 1e-12 and output 1e-12 / final 1e-11 filters are bypassed by no_prune;
- T cached-state 1e-10 filtering is after the no-prune aggregate-and-return branch;
- prune_state and HH/top-level 1e-8 cuts are bypassed by no_prune;
- the V5 coefficient 1e-12 input cutoff is used only in legacy mode.

This is a source-path inspection, not yet a dynamic proof that every execution route is complete. The committed tests are six K5 wiring tests and twelve standalone pruning-helper tests. Their latest blobs are 476255880c123289b7f35d13b4bc81d03f1a4f50 and a1c2cbbb772b814355a520194e4bb32b19aa2383; the workflow blob is ffd43126b2d6003a7e551a5f51e2d5d59bde60de. No test execution result is available in the repository status query, so no green CI claim is made.

No-prune can consume substantially more memory. The CLI catches MemoryError and emits resource_limit without switching modes. It does not yet catch every possible OS-level kill or external timeout; such termination may leave no JSON artifact. A future run wrapper should add wall-time and memory-budget supervision outside the Python process.


## 12. Direct policy-path regression expansion

The K5 wiring suite now contains nine focused tests (up from six). Three new tests exercise behavior rather than only function signatures:
1. a mocked primitive branch with a local intertwiner coefficient of 1e-13 is retained by no_prune and rejected by legacy's 1e-12 local cutoff;
2. two HH input-state contributions are coherently combined before thresholding in no_prune, while the legacy 1e-8 state cutoff remains intact;
3. the same T_cached call arguments with distinct policy modes return distinct values, demonstrating that the mode participates in the cache key.

The tests are intentionally small and deterministic. They test the producer's policy wiring with controlled primitives, not the physical K5 operator's full numerical output. Combined with the 12 helper tests, the committed suite now contains 21 tests. Latest K5 test blob: b83733713aca6bcfa03fdf9212ec7b7f97613e5c. No CI PASS is claimed; combined commit status remains empty. A local clone/test attempt was blocked because this execution environment could not resolve github.com, so no local execution claim is made either.


## 13. Acceptance-gate correction: legacy regression is not a no-prune oracle

A final control-flow review found that the first integration version applied the historical regression targets to both modes. That is logically invalid: the no-prune output is expected to differ from the legacy-thresholded output, so a nonzero process exit based solely on legacy target mismatch would misclassify a completed strict calculation.

The producer now uses apply_acceptance_gate(result, mode):
- legacy: compares the three frozen historical regression targets and sets passed_regression accordingly;
- no_prune: sets calculation_completed=true, leaves passed_regression and regression_errors as null, and labels the scope unpruned_calculation_completed_not_physically_certified;
- unsupported modes are rejected.

The CLI returns zero for a completed no-prune calculation, but that return code means computation completed only; it is not a numerical or physical PASS. Legacy mode still exits according to the historical regression gate. Two additional tests cover this separation, bringing K5 wiring tests to eleven and the combined committed suite to twenty-three.

Latest producer blob: 7db5bae2bb22535cc99f9d86f5dda0f024db7fbc. Latest K5 test blob: 0d43625f3d82f079f3f70c8b547fb4ac88039b8b. The producer, tests and docs have been fetched back from GitHub. No CI execution result or local test execution is available; no test PASS is claimed.
