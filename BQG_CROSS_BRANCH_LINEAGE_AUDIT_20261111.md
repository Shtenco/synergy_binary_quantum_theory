# Cross-branch BQG audit — inventory and deep source tranche — 2026-11-11

## Status and scope

This is a reproducible cross-branch inventory and a deep source audit of the two currently active research frontiers. It is **not** a claim that every line of every historical branch has already been semantically reviewed. The complete branch-head manifest is stored in `BQG_BRANCH_HEAD_INVENTORY_20261111.json`.

The review distinguishes:
1. source facts confirmed from the current branch snapshots;
2. tests committed but not executed;
3. numerical/finite-sector results;
4. mathematical and physical claims that remain open.

No production formula, threshold, operator, or physical claim was changed by this audit record. `main` was not modified.

## 1. Branch inventory and lineage

GitHub branch inventory snapshot:
- branch refs: **168**;
- unique branch-head commits: **157**;
- main HEAD: `254f19fdab41d32b920bf454f987aa7c0fdba34d`;
- causal-cone branch HEAD at snapshot: `78a0d58085bbc69ed310d348348ce6695b1dff70`;
- consolidated research branch HEAD immediately before writing the inventory: `5350eadf30bc3865392db0f5fcf1ab2ff7e8fe2a`.

The inventory itself was committed to the consolidated research branch as:
- `521595504fab369475ca3e3bbb8b4e74ed8b4094`;
- file blob: `cb3e489f3e8a896bca2510bda4e751684ef35c55`.

The inventory JSON records every branch name, its HEAD SHA, duplicate-HEAD groups, and the two principal research refs. Its consolidated-branch SHA is the pre-inventory snapshot by design; the inventory commit cannot include its own resulting HEAD.

Current compare against main after the inventory commit:
- `research/consolidated-bqg-audit-20261111`: **106 ahead / 2 behind**, diverged;
- `research/lemaitre-bqg-causal-cone-20261110`: **141 ahead / 4 behind**, diverged.

These are independent research lines. Do not merge, rebase, or cherry-pick them as one undifferentiated history. The main branch remains untouched.

Tree-size snapshot:
| Ref | Tree entries | Python files | Test files under tests/ | Workflows |
|---|---:|---:|---:|---:|
| main | 724 | 210 | 15 | 35 |
| consolidated research | 739 | 218 | 15 | 36 |
| causal-cone research | 742 | 218 | 20 | 39 |

Counts are source-tree inventory measures, not counts of unique mathematical modules or executed tests.

## 2. Branch-wide Python audit is now stale

The previous `PYTHON_BRANCH_WIDE_AUDIT_2026-11-10.md` reports:
- 157 branch snapshots and 146 unique branch-head commits;
- 18,674 Python-file appearances across snapshots;
- 674 unique Python paths;
- 757 path/content versions;
- 754 unique Python-source blobs;
- 752 AST parses successful and 2 historical syntax-error blobs.

The current inventory has 168 branch refs and 157 unique HEADs. Therefore the previous scan is valuable historical evidence but is **not** a complete scan of the current branch population.

The two syntax-error blobs previously identified were:
- `scripts/pl_lorentzian_group_worker_v2.py:49` on the V23/V24/V25 and independent-proof lineage;
- `tests/test_v21_c1_checkpoint_repair.py:455` on `fix/v21-c1-fail-closed-20261108`.

The previous report states that these blobs were absent from the then-current main snapshot. That does not make the historical branch snapshots valid runnable checkpoints.

The existing scanner is static: it does not execute the 752 parseable source versions. The prior audit also explicitly excluded Python embedded in Markdown/YAML and Python payloads inside tracked archives/base64 bundles. Those embedded sources need a separate inventory and provenance treatment.

The main-branch workflow `.github/workflows/python-all-branches-audit.yml` supports manual dispatch and a weekly schedule, and fetches all remote branch heads. No new full scan was executed in this audit session. A refreshed run must record its exact branch-head map and artifact digest before the new 157-head snapshot can be marked AST-audited.

## 3. C(V)/C(K): confirmed source-level linearity boundary

Audited branch: `research/lemaitre-bqg-causal-cone-20261110`, HEAD `78a0d58085bbc69ed310d348348ce6695b1dff70`.

The active source files were fetched and inspected:
- `scripts/peter_weyl_covariant_composition_gate.py` — blob `1ef229b1dd4b117c890fdce5f7b0025199c5f534`;
- `scripts/peter_weyl_covariant_K_leg_gate.py` — blob `60b9b3ec878cd41383a0d0ac7862db852c9d9552`;
- `scripts/peter_weyl_covariant_K_composition_gate.py` — blob `4f837adb14e5051df85044196e7bd14a9b670d94`;
- `scripts/peter_weyl_covariant_volume_leg_gate.py` — blob `7d1fb7fe364e8e05431dcc109a954365c3274b80`.

### Confirmed risk

`COMP.add` and `CK.add_state` threshold running complex partial sums. The K-composition path also scales basis-column outputs by the actual arbitrary input coefficient before passing them through `COMP.add`. The VH/HV commutator terms and outer hKh/direct composition have further thresholded accumulation boundaries.

Thus a fixed matrix mask and a state-dependent threshold are not interchangeable. The latter can make the production map non-additive, non-homogeneous, and order-dependent near the cutoff. A matrix assembled from unit basis columns does not prove that the production state-to-state function equals that matrix on arbitrary complex superpositions.

The strict accumulator `scripts/bqg_sparse_linear_accumulator.py` is a reference/test oracle. The continuity checkpoint explicitly states that it is **not yet wired into the active V10 KKV production path**.

### Actual L2 test coverage

`tests/test_bqg_l2_production_primitives.py` (blob `7579aa87dbb5f98d9905d9259e814034b1194c7c`) contains:
- three hard fixture-reachability checks;
- three hard safe-domain complex-homogeneity checks;
- three expected-failure threshold-boundary characterizations;
- three independent-support additivity tests.

The full covariant C(V) integration test currently covers only edge (0,1), matrix component (i,j)=(0,0), and two independent Peter–Weyl basis inputs (indices 0 and 16), at common scale 10^6. It does not cover the other three matrix components. There is no C(K) additivity result and no L3 direct-superposition-versus-frozen-matrix reconstruction result.

The expected-failure threshold tests are diagnostic characterization, not accepted linearity. Their exact failure reasons still require actual logs.

### CI state

The L2 workflow dependency-closure correction is present:
- workflow blob `7a8398a40f7efce21587958782f0e81843a9cb83`;
- it watches and compiles the C(V) composition module and the listed Peter–Weyl/charge/volume/projection dependencies.

However, workflow configuration is not execution evidence. At the latest inspected Actions inventory:
- L2 runs #1–#10, including latest run `42940595199`, remained queued;
- L1 runs `42929934053` and `42930769907` remained queued;
- the latest L2 job had no conclusion and no artifact.

No L1/L2 PASS or FAIL, numerical residual, or test log is claimed.

## 4. K5 no-prune and graph-shell audit

Audited branch: `research/consolidated-bqg-audit-20261111`, inventory commit `521595504fab369475ca3e3bbb8b4e74ed8b4094`.

### What has materially improved

The current `scripts/k5_peter_weyl_safe_hda_column.py` source (blob `209a19b466ba28d37fba4d31f9dda0157ec44152`) has an explicit `legacy|no_prune` mode. The no-prune route:
- retains nonzero local recoupling coefficients and primitive outputs;
- coherently aggregates primitive contributions;
- uses mode-aware T-cache identity;
- aggregates T, H, HH and commutator maps at the main composition boundaries;
- avoids silently falling back to legacy mode after `MemoryError`;
- returns `calculation_completed` with `passed_regression=null` and an explicit not-physically-certified scope.

This is a real producer-integration milestone. A zero exit code in no-prune mode means computation completed, not that numerical or physical acceptance passed.

### Confirmed graph-shell defects remain in the current source

`scripts/bqg_k5_graph_sector_master_shell.py` (blob `fa26bedfd33c2be83d8da9c67b7555e49051c038`) still:
1. calls `K5.prune_state(..., 1e-10)` in `apply_H_to_key`, without selecting the producer's no-prune mode;
2. records only scalar per-column outside-shell norm, not the coherent off-diagonal leakage Gram;
3. constructs `M=(M+M†)/2` in `master_matrix`, then measures `norm(M-M†)` in `run`. That Hermiticity diagnostic is tautological after symmetrization and cannot reveal the raw defect;
4. classifies zero modes with `scale=max(max(abs(eigenvalues)),1.0)`, which inserts an absolute unit-scale floor and breaks positive-rescaling invariance for small-norm matrices;
5. accepts a weak structural condition (post-symmetrization Hermiticity, loose minimum-eigenvalue bound, and at least one graph group), which is not a physical-kernel, gap, or convergence certificate.

The exact pullback identity, under a complete action and the correct orthogonal/Gram-weighted Hilbert-space decomposition, is:
`P H† H P = C† C + L† L`.
Omitting `L†L` can create spurious zero modes in the compressed-only master. The identity itself does not require H to be self-adjoint; self-adjointness is a separate physical/operator obligation.

### Leakage-Gram helper is not production-integrated

`scripts/bqg_sparse_leakage_gram.py` (blob `f37e42aa8b05961f1a7f8623465eefa5c38034dd`) correctly encodes the coherent off-diagonal Gram construction relative to supplied complete columns and includes a scale-relative PSD/nullspace diagnostic. Its own input contract requires complete, coherently aggregated columns in the correct orthonormal basis.

But:
- the helper is not wired into the K5 shell producer;
- no dedicated leakage-Gram unittest file is present in the current consolidated tree;
- the helper's embedded self-test uses bare `assert` statements;
- the current `bqg-sparse-pruning.yml` workflow does not execute or path-watch this helper.

Therefore the helper is a useful algebraic kernel, not an accepted production leakage correction.

### Production-level finite tests and their limits

The consolidated branch contains 40 test methods by source inspection across six workflow-invoked scripts. The current workflow's latest inspected runs #26 (`42942768669`) and #27 (`42944352617`) were queued with no artifacts. No test execution or CI PASS is claimed.

The production primitive test file (blob `639120c270eeb88bc730751215e44024d17ba852`) now includes:
- a finite-wall closed-triangle primitive matrix adjoint test;
- local recoupling-basis orthonormality checks;
- 12 oriented-spec combinatorics checks at all five vertices;
- a five-hit matrix-element adjoint test for only `T_sequences(...)[0]`.

The last test is one monomial, one all-spin-1/all-K=0 input, and spin wall 2. It does not establish adjointness of the full coherent 32-monomial T sum or assembled node Hamiltonian.

The oriented-spec test checks that the generator matches a separately written expected list, but both production and test use the same sign rule `(-1)**omitted`. It is an implementation regression check for the declared convention, not an independent derivation of the physical orientation sign.

The raw-volume tests cover four doubled-spin triples: (1,1,1), (1,1,2), (1,2,2), (2,2,2). They check raw Hermiticity, derived-volume positivity/Hermiticity, permutation covariance in one mixed-spin case, and `V²=|Q|`. The production routine still records the raw anti-Hermitian defect and then symmetrizes Q; the standalone HH result does not make this defect a hard acceptance gate over every encountered spin triple.

## 5. Existing finite-sector scientific evidence: preserve its scope

The consolidated audit records a strong finite-sector algebra diagnostic:
- orientation signs restore S5 covariance across the ten node pairs in the stated finite model;
- the projected SO(5) Lie skeleton is exact;
- the full operator bracket has graph leakage `sqrt(296/297) ≈ 0.998315`.

This is evidence of a fixed-sector closure obstruction in that finite diagnostic, not a successful physical HDA closure.

The K5 reachability enumerator reports 120 Gauss-admissible assignments and an 816-dimensional upper bound after one H, and 4193 assignments and a 24364-dimensional upper bound after an HH pair, reaching j=5/2. These are representation-selection counts; they do not prove that every listed state has nonzero amplitude. The j=5/2 wall is scoped to the stated two-hit histories, not arbitrary higher shell depth.

## 6. Dependency-ordered next gates

1. Refresh the all-branch static Python inventory against all 168 refs / 157 unique HEADs, and include embedded code/archive payloads as a separate inventory.
2. Execute and diagnose the queued L1/L2 runs before expanding C(V) claims; then test the other three C(V) matrix components with non-vacuous fixtures.
3. Add independent-support and complex-phase tests for `C_K_component`, then L3 direct action versus a frozen basis-column matrix on cancellation and near-threshold vectors.
4. Thread no-prune mode into graph-shell generation and preserve exact source/policy/cache provenance.
5. Measure raw H_v and master Hermiticity before symmetrization; replace the unit-floor spectral threshold with a scale-relative/backward-error rule and scale sweeps.
6. Add a dedicated unittest suite and CI workflow for the leakage-Gram helper, then wire the helper to complete production columns.
7. Verify the full 32-term T adjoint, all 12 signed orientation terms at operator level, assembled H_v Hermiticity, finite-wall leakage, and only then full HH/leakage-Gram convergence.
8. Keep the physical projector, physical inner product, HDA structure function, Lorentzian normalization and causal cone as open until their independent derivation and convergence gates are met.

## 7. Evidence discipline

- A source inspection is not a test run.
- A queued workflow is neither PASS nor FAIL.
- A finite-dimensional identity is not automatically a continuum or physical theorem.
- A no-prune calculation completing is not a numerical acceptance result.
- A compressed-shell nullspace is not the physical constraint kernel.
- No branch should be merged into main merely because its documentation or CI workflow exists.

**Main remains untouched. No scientific claim has been promoted by this audit.**


## 8. Newer successful branch-wide static scan — run #73

A later successful run was found and its artifact was downloaded and inspected:
- workflow run: [42823223090](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42823223090), run #73;
- run HEAD: `1f9f09c4320abace3939120d8eb8c9af2f992516`;
- artifact ID: `13801636701`;
- artifact SHA-256: `7fae3985481e5ecea86b5d582e3467db3fa931ea9b5ee2cec81039ea549cfd65`;
- job and all audit/upload steps: completed successfully.

This run used the expanded scanner revision that scans standalone Python, explicitly marked Python fences/heredocs, notebooks, and supported Python sources inside tracked archives/base64 archives. It supersedes the earlier 2026-11-10 report counts for that run's 167-branch snapshot.

### Run #73 inventory

- 167 branch refs / 156 unique branch-head commits;
- 21,038 Python-file appearances across branch snapshots;
- 713 unique Python paths, 837 path/content versions, 834 unique Python-source blobs;
- 52,798 tree entries across the 156 unique snapshots;
- 1,383 parsed source units, including standalone, inline, and embedded archive sources;
- 2,327 assertion statements;
- 14 explicit random-seed occurrences.

Embedded/inline coverage:
- 22,895 branch-snapshot inline-source file occurrences;
- 942 unique inline-source blobs and 216 extracted Python blocks;
- 215 archive occurrences, 10 unique archive blobs scanned;
- 983 archive members, 335 embedded Python sources, 1,442,660 embedded-source bytes;
- zero unsafe archive members and zero archive-processing errors.

### Run #73 syntax findings

The run found exactly two source syntax errors:
1. Inline Python heredoc `.github/workflows/lorentzian-pruned-first-column.yml#heredoc-PY-2`, line 86, blob `a0869d6eb4474206015c0d4efaff8e01ed0aff23`, present on branch heads `4cb8ae20b72768d556973361820edb1e3137ae05` and `ab24171a4bb327d190551bdf68957bd10a819c1d`.
2. `scripts/pl_lorentzian_group_worker_v2.py`, line 49, blob `6eb7ba71ddfec1051eb4874399b224768c680372`, present on branch heads `524240db085485ac31678e82778f585ffd38d199`, `b58872dcb23404fd509421e14b975a98d2214a3e`, and `fddd49427c636e92dfe79b96034862aeba629927`.

The previously reported `tests/test_v21_c1_checkpoint_repair.py:455` syntax error is not present in this later scan; the `fix/v21-c1-fail-closed-20261108` branch had advanced to a corrected head before run #73. The historical report's two-error list must not be copied forward unchanged.

### Main snapshot in run #73

The run scanned main HEAD `23669209c83fe0aa2741612e5b2b119562ee6d7d`:
- 210 Python files;
- zero standalone syntax errors;
- zero inline Python syntax errors;
- `main_ast_status=PASS`.

Main-only static findings included 250 assertion nodes, 8 broad exception handlers, 7 process-execution call sites and 2 network call sites. The run reported no main-branch `pickle.load` or unguarded archive-extraction finding. These are static pattern counts, not proof that every call site is unsafe or safe.

Across all branch snapshots the scanner reported 100 broad-exception handlers, 55 process-execution call sites, 32 pickle-deserialization call sites, 17 network call sites, 12 archive-extraction call sites and one bare `except`. The bare handler is at `aggregate_existing/bqg_s4sign_aggregate_capacity.py:61` on the `research/depth6-s4sign-proof` snapshot. It catches `BaseException` while parsing a filename and should be narrowed to explicitly expected parse errors. Broad exception counts require semantic review; worker-boundary serialization is not automatically a defect.

### Why another refresh is still required

Run #73 predates the current branch census:
- run #73: 167 refs / 156 unique heads;
- current inventory: 168 refs / 157 unique heads;
- the additional branch is `research/consolidated-bqg-audit-20261111`, created after run #73;
- current main HEAD is `254f19fdab41d32b920bf454f987aa7c0fdba34d`, later than run #73's main HEAD.

Thus run #73 is the latest confirmed complete static baseline, including inline/archive coverage, but it is not a scan of the current 168-ref snapshot. The next branch-wide run must be tied to the current branch-head manifest and current main SHA. Do not claim current 168/157 coverage until that artifact exists.
