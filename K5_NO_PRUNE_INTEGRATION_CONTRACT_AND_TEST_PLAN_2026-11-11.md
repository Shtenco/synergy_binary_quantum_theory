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


## 14. Independent synthetic dense-oracle and adjoint harness

Added scripts/test_bqg_sparse_dense_oracle.py and wired it into the sparse-pruning workflow. The test builds a deterministic 5-dimensional complex sparse operator with duplicate matrix-element contributions, a tiny nonzero amplitude, and an exact coherent cancellation. It independently assembles the dense matrix without calling the sparse aggregation helper, then checks:
1. sparse assembly versus dense matrix assembly;
2. sparse matrix-vector action versus dense matrix-vector multiplication;
3. the complex Hilbert-space adjoint identity <x,Ay> = <A†x,y>, using the independently assembled dense matrix;
4. retention of a tiny nonzero coefficient and removal of an exactly cancelled coefficient.

Scope boundary: this is a synthetic algebraic oracle for sparse aggregation/action and complex adjoint conventions. It is not a dense oracle for the production K5 Hamiltonian, the Peter-Weyl volume tensor, or the physical K5 column. The production-operator oracle and adjoint checks remain open.

The workflow now runs three scripts: 12 sparse-pruning helper tests, 11 K5 policy-wiring tests, and 4 synthetic dense-oracle/adjoint tests, for 27 committed test methods total. No test execution or CI PASS is claimed until logs are available.


## 15. Production K5 raw-volume operator exposure and small-sector tests

The K5 producer now exposes volume123_raw_operator(s1,s2,s3), which returns the unsymmetrized triple-product Q before the existing Hermitian completion. volume123_matrix now consumes this shared raw operator, records the pre-symmetrization anti-Hermitian defect, and then retains the historical completion and sqrt(abs(Q)) construction. This makes the raw object directly testable rather than inferring its properties from the symmetrized output.

Added scripts/test_k5_volume_raw_operator.py with three test methods over doubled-spin triples (1,1,1), (1,1,2), (1,2,2), (2,2,2):
1. raw Q is Hermitian before symmetrization to a relative defect <= 5e-14;
2. the derived volume matrix is Hermitian and positive semidefinite within stated floating-point tolerances;
3. the recorded raw defect and raw norm agree with independently recomputed diagnostics from Q.

These are tests of the actual production K5 spin/volume routines in small finite sectors, not a full HH computation and not a proof of continuum volume properties. The workflow now includes this test file. Committed test-method inventory: 12 sparse helper + 11 K5 wiring + 4 synthetic dense oracle/adjoint + 3 production raw-volume tests = 30. No test run or CI PASS is claimed in this update.


## 16. Production link-hit adjoint controls

Added scripts/test_k5_hit_adjoint.py. It checks the actual production Clebsch-Gordan link-hit matrices for the j=1/2 -> j=1 transition in all four binary endpoint channels. Reversing edge orientation and swapping endpoint channel labels must return the conjugate-transpose matrices and the same real normalization. Additional structural tests verify that path reversal and sequence adjunction are involutions.

Scope: this is a local link-hit adjoint control. It does not yet establish that the complete primitive, T, or node-H maps are adjoints on the full coupled K5 state space. Those require independent assembled matrix-element tests, including coherent sums and all recoupling channels.

The CI workflow now runs five test scripts: 12 sparse-helper + 11 K5 policy-wiring + 4 synthetic dense-oracle/adjoint + 3 raw-volume + 3 link-hit-adjoint = 33 committed test methods. No test run or CI PASS is claimed.


Workflow wiring review: all five test scripts are now present in both push and pull_request path filters, and all five have explicit execution steps. This closes a CI-configuration omission where the newest test files were initially present only in push filters. The final workflow blob SHA is 7caec171b996b3c662b4445681764ff24f1d2a1f. The workflow's actual run result is still not visible; configuration correctness is not a CI PASS.


## 23. First production K5 primitive dense-adjoint oracle — 2026-11-11

Added `scripts/test_k5_primitive_dense_adjoint.py`. Unlike the synthetic sparse/dense harness and the local CG hit-matrix test, this test calls the actual production `apply_primitive_oriented` routine and independently materializes two finite-sector matrices:

- forward primitive on the closed triangle path `(0,1,2,0)`;
- reversed primitive generated by the production adjoint-sequence convention.

The input basis contains all five-node recoupling-channel combinations with all ten links at doubled spin 1. The output basis explicitly spans all eight assignments of doubled spins 0 or 2 on the three triangle edges, with the other seven links fixed at doubled spin 1 and all allowed node recoupling channels included. The closed loop changes two incident links at each affected node, preserving the parity needed for gauge-invariant intertwiner sectors. The finite spin wall is `Jmax2=2`; thus the test is explicitly a finite-wall adjoint comparison, not an unrestricted-Hilbert-space theorem.

For each of the four binary endpoint-channel pairs, the test compares the independently materialized reverse matrix with the conjugate transpose of the forward matrix. It fails closed if the production action produces a key outside the declared output basis. This is the first committed test in this line that probes a complete coupled K5 primitive rather than only its individual CG link-hit factors.

**Evidence boundary:** the test and CI wiring are committed, but neither the new test nor the workflow has yet been executed or returned a CI result. A passing result would certify only this one closed-loop primitive, this finite spin wall, this basis convention and these channel sectors. It would not certify every primitive path, each T term, the assembled node-H operator, full HDA closure, or physical continuum behavior.

The workflow now includes this sixth test script. The committed method inventory is 12 sparse-helper + 11 K5 policy-wiring + 4 synthetic dense-oracle + 3 production raw-volume + 3 local link-hit adjoint + 1 production primitive dense-adjoint = **34 test methods**. No PASS is claimed.


## 24. Orthonormal-basis prerequisite added to the production primitive oracle — 2026-11-11

The finite-sector dense-adjoint test now independently checks the local recoupling basis before comparing operator matrices. New method:
`test_local_recoupling_bases_are_orthonormal_on_all_test_sectors`.

For the all-doubled-spin-1 input sector and each of the eight output assignments on the three closed-triangle edges, it constructs every local oriented intertwiner tensor at all five vertices and verifies the full local Gram matrix against the identity. This guards a crucial hidden assumption: a matrix equality with a conjugate transpose is the correct adjoint test only when the chosen domain and codomain bases use the declared orthonormal inner product.

The added check does not itself certify the global tensor-product basis or physical graph-changing inner product; it verifies the local recoupling factors used to form that basis. The production primitive matrix comparison remains a separate test and is not weakened or symmetrized.

The production primitive test file now contains two test methods. Total committed inventory: **35 test methods** across six scripts. The workflow already invokes the file, so no workflow command change is needed. Neither the expanded test nor CI has been executed/confirmed; no PASS is claimed.


## 25. Derived-volume spectral covariance gate — 2026-11-11

Added a direct test of the production spectral map (V=\sqrt{|Q|}), separate from raw-(Q) Hermiticity and alternating-tensor covariance. For mixed spins `(1,2,2)), all six tensor-factor permutations are represented by explicit permutation unitaries. The test requires
\[
V_\pi=U_\pi VU_\pi^\dagger.
\]
The parity sign in the raw (Q) transformation disappears under the absolute-value spectral map. This catches a mismatch in the derived volume reconstruction that a raw-Q-only test would not detect.

**Inventory update:** `scripts/test_k5_volume_raw_operator.py` now has five test methods; the six workflow scripts contain **37 committed test methods** by source inventory. The existing workflow already invokes the volume test file. This gate is committed but not yet execution-verified. It checks one mixed-spin finite sector and does not certify all spin triples, epsilon-orientation conjugations, the four-valent intertwiner restriction, or the full HH calculation.


## 26. Volume square-root spectral reconstruction gate — 2026-11-11

The volume suite now checks the defining spectral identity V²=|Q|, not only Hermiticity and nonnegative eigenvalues of V. For each of the four existing finite spin triples, the reference |Q| is assembled from a fresh eigendecomposition of the Hermitian part of the raw production Q, then compared with the square of the production volume123_matrix result.

The volume test file now has six methods; the six workflow scripts contain **38 committed test methods** by source inventory. This new gate is committed but not execution-verified. It tests finite spectral reconstruction and does not establish the physical correctness of the chosen triple-product operator, its orientation convention, or full K5 HH closure.


## 27. Complete T-monomial production adjoint gate

Added a matrix-element adjoint check for one complete five-hit monomial from the actual `T_sequences` generator. It includes the local volume insertion and three path factors, then compares every nonzero forward output coefficient with the corresponding return coefficient from the reversed production sequence.

Acceptance identity:
\[
\langle y|A|x\rangle=\overline{\langle x|A^\dagger|y\rangle}.
\]

Scope is explicitly finite: one monomial, one all-spin-1/all-(K=0) input ket, and doubled-spin wall 2. This does not certify the coherent 32-monomial T sum, all orientation triples, or the assembled node Hamiltonian. The test is committed and included through the already-existing primitive test workflow step; no execution or CI PASS is claimed.


## 28. Orientation enumerator audit and complete achievement register update — 2026-11-11

The production orientation generator `oriented_specs(v)` is now covered by an independent combinatorics test at all five K5 vertices. The test derives the expected ordered list without calling the generator's internal construction: for each omitted-neighbor index (r=0,1,2,3), it preserves the inherited order of the remaining triple, emits its three cyclic rotations, and assigns common sign ((-1)^r). It checks exact sequence equality, cardinality 12, and uniqueness.

New test method:
`test_all_twelve_oriented_specs_have_exact_combinatorics`.

New audit artifact:
`K5_ORIENTED_SPEC_COMBINATORICS_AUDIT_2026-11-11.md`.

### Current test inventory

The six workflow scripts contain 40 methods by source inspection:
- sparse pruning helper: 12;
- K5 pruning/policy wiring: 11;
- synthetic dense oracle: 4;
- production raw-volume operator: 6;
- production link-hit adjoint: 3;
- production primitive/orientation: 4.

This count is not an execution result. GitHub Actions run #27 (`42944352617`) was queued at the last check; run #26 (`42942768669`) was also queued. No local run and no CI PASS are claimed.

### Remaining hard acceptance gates

1. Full coherent 32-monomial T matrix-element adjointness, including all coefficient signs and interference.
2. Assembled (H_v=(T_v+T_v^\dagger)/2) Hermiticity with all 12 signed orientation specifications.
3. Explicit spin-wall leakage/clipping accounting and stability as (J_{\max}) grows.
4. Full HH composition in no-prune mode, with resource accounting.
5. Full leakage Gram (sum_vL_v^\dagger L_v), including off-diagonal interference and the correct graph-changing Hilbert inner product.

The legacy path remains preserved. No-prune completion is not a physical acceptance result. Do not merge or rebase this research branch onto `main` without explicit approval.


## 29. Full coherent T-sum adjoint acceptance gate — 2026-11-11

Added `test_full_32_term_T_sum_matches_reverse_matrix_element`. It evaluates all 32 production T monomials for a fixed finite-sector input, accumulates their signed forward contributions independently, then evaluates the corresponding reverse matrix element by applying all 32 reversed production sequences. The test checks one deterministic, largest-magnitude nonzero matrix element:
[
\langle y|T|x\rangle=\overline{\langle x|T^\dagger|y\rangle}.
]

This closes an important gap between a single-monomial check and the coherent T sum. It does not yet establish matrix-wide adjointness: only one selected element, one input ket, one orientation triple and one spin wall are covered.

Source inventory is now 41 methods across six workflow scripts. This is not a test execution result. CI status must be verified from the actual Actions job; no PASS is implied by commit presence.


## 30. Clarification of full T-sum test scope and CI state

The 32-monomial test now checks one selected matrix element of the complete coherent T sum. It is not yet a matrix-wide T/T† adjoint test. The next gate must expand to multiple input columns and a declared finite domain/codomain basis, with explicit handling of basis states clipped by the spin wall.

At the latest check, Actions runs #26 (`42942768669`), #27 (`42944352617`) and #28 (`42944740425`) were all queued; no test PASS is claimed. Run #28 is attached to code commit `ba4dfc48db4f1a2ea5caea8ef7a3160fc63bdd00`. Later commits only update audit documentation.


## 31. Local oriented-volume adjointness gate

The production volume suite now tests self-adjointness of `apply_volume_tensor_oriented` itself, not merely Hermiticity of the un-oriented spectral matrix. For deterministic complex tensors at all five vertices and all four declared spin triples, it checks the Frobenius identity
[
\langle A,V_vB\rangle_F=\langle V_vA,B\rangle_F.
]
This isolates the epsilon-conjugation convention used to orient the local volume action.

The volume test file now has seven methods. Total workflow source inventory is 42 methods. This remains finite local tensor-space coverage and is not a full primitive/H adjoint certificate. Test commit: `50c807d03c8616dce639487afc83704ec19c5b36`; execution and CI PASS are not claimed.


## 32. Matrix-wide T/T† and assembled H acceptance protocol — 2026-11-11 02:18 UTC

### Current source and CI checkpoint

The current six workflow-invoked scripts contain 42 test methods by source inspection (12 + 11 + 4 + 7 + 3 + 5). The latest observed branch head before this documentation update was `3e4b45e944ecbe92e731b2d0e7cd497d9665afd2`; comparison against `main` was diverged, 109 ahead and 2 behind. The branch remains isolated; no merge/rebase/main write was performed.

Actions runs #26–#29 (`42942768669`, `42944352617`, `42944740425`, `42944880600`) were checked through their jobs. All four jobs remained queued, with no conclusion or steps. No local test execution and no CI PASS are confirmed. The latest code test commit is `50c807d03c8616dce639487afc83704ec19c5b36`.

### Achievements to preserve as separate evidence classes

- deterministic coherent sparse aggregation and no-prune execution contract;
- raw volume defect capture before historical Hermitian completion;
- raw alternating-tensor covariance under all six leg permutations;
- derived `sqrt(|Q|)` covariance under all six permutations;
- independent spectral reconstruction `V²=|Q|` on four finite spin triples;
- local oriented-volume tensor self-adjointness at all five vertices and four triples;
- local link-hit adjoint tests and local recoupling-basis orthonormality;
- production closed-loop primitive dense-adjoint matrix oracle;
- exact 12-spec orientation combinatorics at all five vertices;
- one complete five-hit T monomial adjoint matrix-element check;
- full coherent 32-term T/T† check for one selected matrix element.

These do not yet establish a matrix-wide T adjoint, assembled node-H Hermiticity, unrestricted operator adjointness, full HH, HDA closure, or continuum physics.

### Required matrix-wide T/T† protocol

The current 32-term gate uses one input ket and one selected output matrix element. Its next extension must declare domain and codomain bases explicitly and compare rectangular matrices, not only a chosen element. The finite test should use multiple orthonormal input columns, coherent no-prune aggregation of all 32 signed monomials, and an independently evaluated reverse action from all 32 adjoint sequences.

For every column, record the full set of output keys. For the reverse action, record outputs outside the declared input domain. Do not silently project away or drop escaped keys: report forward and reverse boundary leakage separately, and require zero leakage only when the chosen subspaces are intended to be closed. If the subspaces are intentionally rectangular/non-invariant, state that the equality is only the projected finite-subspace identity and do not promote it to unrestricted adjointness.

Keep the test small enough for the existing five-minute workflow timeout. Before increasing the number of columns, estimate the number of production primitive applications (32 per forward column and 32 per reverse column per tested output ket) and use a deterministic, mathematically justified basis selection. Do not weaken tolerances to accommodate a mismatch.

### Assembled node-H gate

After the constituent T gate, independently materialize the signed 12-spec sum for `H_v=(T_v+T_v†)/2`. Compare its matrix with its conjugate transpose on the declared finite sector. Preserve the per-specification sign and ordering; do not replace the production assembly with an already-symmetrized matrix. Report boundary leakage and raw anti-Hermitian defect before any diagnostic symmetrization.

### Subsequent sequence

1. matrix-wide coherent T/T† with boundary accounting;
2. assembled 12-spec node-H finite-sector Hermiticity;
3. spin-wall ladder and clipped-transition ledger;
4. full no-prune HH composition with resource-limit manifest;
5. complete leakage Gram `Σ_v L_v†L_v), including off-diagonal interference and the correct graph-changing inner product.

No-prune completion remains a computation status, not a physical acceptance result. Preserve the legacy path and never use fallback to legacy to conceal a resource limit.


## 33. First 2×2 coherent production T/T† block gate

Added `test_full_32_term_T_two_by_two_matrix_block_matches_reverse` in commit `3bb5e42bc946791ae96ecc6f57d041a4c1fd7554`; the follow-up commit `92ab5397876aabfbe509b179d8efaa56ec64ad59` prints its leakage ledger into CI logs. The test uses two distinct input recoupling basis kets, all 32 signed T monomials, coherent no-prune accumulation, two deterministic output rows, and independently accumulated reverse columns. It compares the resulting 2×2 block to the conjugate transpose of the forward block.

The gate also reports forward amplitudes outside the selected output rows and reverse amplitudes outside the selected input kets (support count, L2 norm, max absolute value). These are projected-subspace leakage metrics, **not** finite-spin-wall clipping metrics. It deliberately does not require these values to vanish: the selected 2D subspaces are not claimed invariant.

The primitive/orientation script now has 6 methods and the six workflow scripts have **43 methods by source inspection**. The existing workflow already invokes this script. The source was fetched back from GitHub and inspected. A local test execution could not be started because the container could not resolve github.com; Actions run discovery for the push commit returned no run, so execution/CI PASS is unconfirmed.

**Interpretation:** this is the first multi-input/multi-output projected block test for the full coherent T sum. It is stronger than the earlier single-element test, but still not a full matrix-wide finite-sector adjoint certificate. Next gates: observe runtime and numerical output; expand to complete declared finite domain/codomain matrices with wall clipping accounting; then test the assembled signed 12-spec H_v.
