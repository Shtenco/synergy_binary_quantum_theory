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

The scanner flagged two parse errors, but manual source inspection classifies them differently:
1. `scripts/pl_lorentzian_group_worker_v2.py`, line 49, blob `6eb7ba71ddfec1051eb4874399b224768c680372`, is a **confirmed Python IndentationError**. The metadata assignments after `d=json.loads(...)` are indented one level too far. The blob is present on branch heads `524240db085485ac31678e82778f585ffd38d199`, `b58872dcb23404fd509421e14b975a98d2214a3e`, and `fddd49427c636e92dfe79b96034862aeba629927`.
2. Inline Python heredoc `.github/workflows/lorentzian-pruned-first-column.yml#heredoc-PY-2`, line 86, blob `a0869d6eb4474206015c0d4efaff8e01ed0aff23`, is a **template-dependent parser false positive**, not a confirmed workflow syntax failure. The snippet contains assignments such as `idx=${{ matrix.idx }}`; GitHub Actions substitutes these expressions before the shell executes the heredoc. The raw pre-expansion block is not valid Python, but the rendered block can be valid. The same blob occurs on branch heads `4cb8ae20b72768d556973361820edb1e3137ae05` and `ab24171a4bb327d190551bdf68957bd10a819c1d`.

The scanner's `iter_inline_python_blocks` currently parses the raw heredoc without first accounting for GitHub Actions `${{ ... }}` interpolation. Its tests cover Markdown fences and plain Python heredocs, but not templated heredocs. The inline finding should therefore be classified as template-dependent until a template-aware parser or a rendered-workflow test is added.

The previously reported `tests/test_v21_c1_checkpoint_repair.py:455` syntax error is not present in run #73; the `fix/v21-c1-fail-closed-20261108` branch had advanced to a corrected head. Thus run #73 contains one confirmed Python syntax defect plus one scanner false positive, not two confirmed source defects.

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


## 9. Audit-scanner precision issue found during manual syntax triage

The manual check of the run #73 inline syntax finding found a limitation in the scanner itself:
- `iter_inline_python_blocks` extracts raw Python heredoc content from YAML/shell files;
- it does not normalize GitHub Actions template expressions such as `${{ matrix.idx }}`;
- the AST parser then reports a syntax error on an unrendered template expression.

The exact flagged workflow assigns `idx`, `b`, and `c` from `${{ matrix.* }}` values. GitHub Actions substitutes those values before the shell launches Python. This is therefore a false positive for runtime syntax, although it is correct that the unrendered template text is not standalone Python.

The existing inline-audit tests cover ordinary Python fences/heredocs, indentation and multiple fences, but not GitHub Actions interpolation. A future scanner improvement should:
1. distinguish ordinary Python blocks from GitHub-templated blocks;
2. substitute syntactically safe typed placeholders or run a rendered-workflow fixture before AST parsing;
3. retain the original source line map;
4. report unresolved template syntax separately from confirmed Python syntax errors.

The other run #73 finding, `pl_lorentzian_group_worker_v2.py:49`, remains a genuine syntax error. No historical source was modified as part of this audit.


## 10. Independent-control matrix transcription defect — Habr audit branch

Audited branch: `audit/habr-emergent-graph-independent-control-20261111`, commit `553c2936a2c53352dfd12e4645255124422584a5`.

The commit message says it preserves exact rational half-integer entries, but the actual diff changes six numeric coefficients in `audit_integer_kernel_full_matrices.py` incorrectly. The parent version (blob `76b2fe7565e8884d71bb2f7ca0a6c64dfaf6d76c`) had the original values as floats; the current version (blob `71625cd4a0e6c0daf0cdc9603cac8192c82bda8a`) converts them to different rationals:

| Matrix row | Previous value | Current expression | Current value | Correct exact conversion |
|---|---:|---|---:|---:|
| Ksi_0 | 1/2 | `sp.Rational(0,2)` | 0 | `sp.Rational(1,2)` |
| lambda_B0 | 1/2 | `sp.Rational(0,2)` | 0 | `sp.Rational(1,2)` |
| Ksi_minus | 1/2 | `sp.Rational(0,2)` | 0 | `sp.Rational(1,2)` |
| phi_meson | 1/2 | `sp.Rational(0,2)` | 0 | `sp.Rational(1,2)` |
| rho_meson | 5/2 | `sp.Rational(4,2)` | 2 | `sp.Rational(5,2)` |
| K_star_meson (70-row matrix only) | -5/2 | `sp.Rational(-4,2)` | -2 | `sp.Rational(-5,2)` |

The 69-row matrix therefore has five corrupted coefficients; the 70-row matrix has six. Since the reproducer constructs `A=(2*M.T).applyfunc(int)`, each error changes an integer entry of A by exactly one. This is not a formatting-only difference: the input integer matrices for the Smith decomposition and both LLL calculations have changed.

The committed `FINDINGS_019.md` and `STATUS.md` still report the earlier exact SNF/LLL results, including saturated-kernel indices `2^179` and `2^180`, and explicitly state the K-star row as `(2,0,-5/2,0,6,-1)`. Those results were produced against the parent matrix, not the current script. They are not verified for the current blob.

The script's internal assertions verify rank, kernel multiplication, unimodularity and lattice identities **for whatever matrix is currently transcribed**. They do not verify that the matrix entries match the frozen external source revisions. There is no row-by-row source digest or independent matrix fixture. Consequently, a transcription error can pass the script's internal algebra checks while invalidating the external-source claim.

No Actions run was found for this branch in the inspected run inventory. The branch's `core-regression.yml` is restricted to pushes and pull requests targeting `main`, so it does not automatically validate this branch's script.

**Disposition:** quarantine the current `audit_integer_kernel_full_matrices.py` result. Restore the six exact values above in a separate corrective change, rerun the full 69×6 and 70×6 SNF/LLL calculation, compare every matrix row against the frozen source blobs, and regenerate `FINDINGS_019.md` and its JSON artifacts. Do not silently overwrite the prior result; preserve both source versions and record the corrected output as a new finding.


## 11. Exact recomputation of the Habr matrix variants

To quantify the transcription damage, the 69×6 and 70×6 matrices were reconstructed from the parent-version values and from the current committed values, using exact SymPy rationals. The same calculation was applied to `A=2M^T), Smith normal decomposition, saturated integer kernels, row-wise primitive kernels and LLL (`delta=3/4`). SymPy 1.14.0 was used; this was a local exact computation, not a GitHub Actions run.

| Matrix/version | Rank | SNF diagonal of 2Mᵀ | Primitive-row index in saturated kernel | Power of 2 | Max coefficient: saturated LLL / primitive LLL |
|---|---:|---|---:|---:|---:|
| M69, parent values | 6 | (1,2,2,2,2,2) | 2^179 | 179 | 2 / 16 |
| M69, current committed values | 6 | (1,2,2,2,2,2) | 2^179 | 179 | 2 / 16 |
| M70, parent values | 6 | (1,2,2,2,2,2) | 2^180 | 180 | 2 / 16 |
| M70, current committed values | 6 | (1,2,2,2,2,2) | 2^181 | 181 | 2 / 16 |

The unchanged SNF diagonal and unchanged maximum LLL coefficient do **not** rescue the transcription:
- For M69, the saturated kernels are different despite the same rank, SNF diagonal, index and max-coefficient summary. Five of 63 basis vectors from each version fail the other version's exact kernel equation.
- For M70, the primitive-row index doubles from 2^180 to 2^181. Six of 64 basis vectors from each version fail the other version's exact kernel equation.
- Thus the summary statistics in `FINDINGS_019.md` happen to remain the same for M69, but the actual kernel basis/subspace is not the same; the M70 index is explicitly stale.

This independently confirms that the current script and the frozen result document are inconsistent. The exact matrix entries, not just rank and a few summary invariants, must be treated as part of the certificate.

A further fail-open issue: `audit_integer_kernel_full_matrices.py` uses Python `assert` for its rank, exact-kernel, unimodularity and basis checks, then prints literal `PASS` labels. Running it with `python -O` removes those assertions while leaving the PASS output. It has no `__debug__` refusal or explicit runtime `require` helper. This script must not be used as a certification gate under optimized Python until that is corrected.

**Revised disposition:** the parent-version computation reproduces the recorded 2^179/2^180 indices. The current committed M70 script reproduces 2^181, and both current matrices use different kernels. Correct the rational entries, add exact row-fixture checks, replace optimization-sensitive assertions, rerun both complete matrices, and only then publish a new result revision.


## 12. Tetrahedral charged-volume branch: CI failure is a near-zero relative-error false positive

Audited branch: `audit/pl48-volume-failclosed-20261110`, HEAD `b26ae55f64e3f9646c91735f005a35b30b4bb34b`. The current backend blob is `c50e1e8597c21a1624c864da1abc5ea56c8b1f97`.

### Actual CI evidence

The latest inspected `tetrahedral-charged-volume-audit` run was `42848173310` at source HEAD `7edee78d2e06a2486fc6f614570f7c8268bb4473`. Its backend blob is identical to the current branch backend. It failed in `test_tetrahedral_volume_spectral.py` before the uniqueness, S4 covariance, and full charged-volume audit steps could execute.

The traceback reports:
`raw tetrahedral Q leaks out of its total-(J,M) sector: absolute=5.551115123125783e-17, relative=1.0`.

A second run, `42848284454`, at the current branch HEAD, also failed in its preflight `Tetrahedral charged-volume audit` step:
`absolute=4.9065389333867974e-18, relative=1.0`.

The absolute residuals are at or below ordinary double-precision roundoff. The backend computes:
`relative = err / max(norm(QB), 1e-30)`
and rejects whenever this relative value exceeds `1e-10`. When the mathematically expected output `QB` is zero or nearly zero, a roundoff residual of order 10^-17 becomes a relative defect of 1.0. The code does not apply an absolute tolerance to this near-zero case, and the exception omits the spin tuple, J2 and M2 that produced it.

This is a **numerical gate defect**, not evidence of a 10^-17 absolute operator leakage. It is not a PASS either: the test is currently red and the remaining gates are skipped.

### Required correction

Use a mixed backward-error condition, for example:
`err <= atol + rtol * norm(QB)`,
where `atol` is derived from machine epsilon, the operator/block scale and the norm of the normalized basis input. Record `norm(QB)`, `err`, `atol`, `rtol`, `spins_local`, `J2`, and `M2` in every failure. Do not simply raise the tolerance until the test turns green. Add a deliberate exact-zero-output fixture and a nonzero-output leakage fixture so the two regimes are separately tested.

The local file `verification_results/TETRAHEDRAL_CHARGED_VOLUME_AUDIT_LOCAL.json` says `passed=true`, but it is explicitly a prior local prototype and does not override the failing current CI runs. No CI artifact was produced by the failing runs.

### Uniqueness theorem: correct result, numerical rank implementation

The uniqueness gate constructs an integer 96×4 constraint matrix but computes its rank using `numpy.linalg.svd` with a floating threshold of `1e-12`. The markdown describes the rank result as exact; the implementation is numerical.

An independent local SymPy exact-integer recomputation gives:
- exact constraint rank: 3;
- exact nullspace generator: `(-1,1,-1,1)`, projectively `(1,-1,1,-1)`;
- zero failures of the representation homomorphism over all 24×24 permutation pairs;
- exact zero defect for all 24 alternating-character identities.

Thus the one-dimensional S4 pseudoscalar result is independently supported, but the executable rank gate should use exact rational/integer rank and nullspace rather than SVD.

### Additional source-level precision issues

- `TETRA_Q_DIAGNOSTICS["raw_frobenius_norm"]` stores `max(norm(Q),1.0)`, not the raw Frobenius norm. Preserve the actual norm separately from the tolerance scale.
- M-independence and `V^4-Q^2` diagnostics also use unit-floored scales. Report absolute and genuinely relative residuals separately so scale dependence is visible.
- The current workflow's preflight failure prevents the later S4 covariance and Euclidean-invariance checks from being rerun in CI on the current code. The local JSON values remain historical until the workflow is green.

**Disposition:** retain the tetrahedral completion as a strong finite-sector candidate and retain the exact uniqueness result, but keep the charged-volume production gate RED until the mixed-tolerance correction, exact-rank gate, diagnostic context and fresh artifacts are committed and executed.


## 13. Lorentzian ordered-triple linearity workflow: repeated vacuous zero-output failure

Audited branch: `audit/pl48-volume-failclosed-20261110`.

The latest inspected `lorentzian-superposition-linearity` run was `42839123248` at HEAD `36b8cd90466df79d6bab18c36722f7379bbb3b06`. Its script blob `2067bc9e7e7c9822ea349dbe55cf66151910ff63` is unchanged at the current branch HEAD.

The artifact `PETER_WEYL_LORENTZIAN_SUPERPOSITION_LINEARITY.json` reports:
- `passed=false`;
- `superposition_support=0`, `column_sum_support=0`;
- `superposition_norm=0`, `column_sum_norm=0`;
- `relative_linearity_error=0.0`;
- `superposition_output_nonzero=false`;
- `C_K` cache: zero hits, zero misses, zero entries;
- `HE_reduced` cache: zero hits, zero misses, zero entries;
- C(V) and inverse caches were exercised, but the K-composition path was not reached.

The gate correctly fails because the nonzero-output requirement is false. The reported zero relative error is vacuous: both sides of the linearity comparison are the zero map on this chosen fixture. Since `ordered_triple_covariant_from_gauss` skips the K legs when the first C(V) output is empty, the run provides no evidence about C(K) or the full Lorentzian ordered triple.

The current fixture is hard-coded to source node 0, ordered edges (2,3,4), and the first two entries of `PW.basis_full_jhalf()`. Those choices produce no output at `Jmax=5/2`. The correct next step is a deterministic fixture search over reachable basis keys and valid ordered triples, then freeze a nonzero fixture with:
- nonzero C(V) support and norm;
- nonzero first C(K) support and norm;
- nonzero second C(K) or final scalar Gauss output;
- a complex superposition whose output remains nonzero;
- direct-vs-column reconstruction checked with a strict accumulator and fixed support.

The current helper `add` prunes every accumulated value below `TOL=10^{-10}`; this is another state-dependent numerical boundary in a test whose purpose is to certify linearity. The test should use a no-prune/stable reference accumulator, or a fixed matrix mask, and separately characterize the production threshold.

The workflow `.github/workflows/lorentzian-superposition-linearity.yml` also has an incomplete path filter. The test imports and calls the K5 producer, C(V), ordered-triple, C(K), and recoupling modules, but the workflow watches only the superposition script, logical-projection script, Gauss-action script and its own YAML. Changes in important operator dependencies can therefore fail to trigger this test.

**Disposition:** this is not evidence that the physical operator is nonlinear. It is a correctly detected non-vacuity/fixture failure, and it currently does not exercise C(K). Keep the Lorentzian superposition claim unproven until a nonzero end-to-end fixture and dependency-closed CI run pass.


## 14. Collective FEM spectral-dimension precursor: useful independent evidence, not a killer pass

The `audit/pl48-volume-failclosed-20261110` branch contains an intrinsic tetrahedral FEM spectral-dimension calculation on the canonical barycentric PL-S3 carrier. The method is distinct from graph-hop fitting:
- local FEM stiffness is assembled from each tetrahedron's intrinsic Gram matrix;
- a lumped mass matrix defines the symmetric Laplacian `M^{-1/2} K M^{-1/2}`;
- L0–L2 are diagonalized directly;
- L3 (40,256 vertices / 221,184 tetrahedra) uses 24 deterministic Rademacher/Hutchinson probes;
- the L3 time window is predicted from the L2 peak using h² scaling and scanned on a preregistered fixed interval.

The checked-in local artifact reports:
`D_peak = [1.9085210169, 2.3877002596, 2.7553265235, 2.8529008753 ± 0.0217976734]`.
Total intrinsic PL volume is conserved to relative spread `1.9984e-15`. The L3 central value remains `0.1470991247` away from 3, so the declared static `0.10` threshold is **not met**.

The code correctly marks this as `PRECURSOR_ONLY`, says not to populate `D_space_metric`, and does not promote a static background observable to the dynamic E/S/R collective killer. Its `passed=true` refers only to internal checks (volume conservation, monotone peak sequence, L3>L2 and finite uncertainty), not to the GR target comparison.

The artifact is local, not a CI artifact on this branch. The workflow's push trigger is restricted to `research/bcqg-core-candidate-v1`; on this branch it would require manual `workflow_dispatch`. No matching Actions run was found for `audit/pl48-volume-failclosed-20261110`.

**Disposition:** retain as a useful independent static FEM precursor and a nontrivial monotone refinement trend; do not call it a three-dimensionality proof, a killer pass, or a dynamical collective result. More probes and independent seed checks can improve the uncertainty estimate, but cannot remove the current 0.1471 central-value miss by themselves.


## 15. Q2 reversible history dilation: lower bound is valid, uniqueness claim is too strong

Audited ref: `research/q2-unified-complex-structure`, HEAD `83f45dfc19cdef884e5290f111bf46dc0e10d1e7`.

The information-preservation argument is valid under its stated premise. If four orthogonal active states are mapped isometrically to one common rank-one no-link state tensored with an environment, then the four environment states must be mutually orthogonal:
`<e_k|e_l>=delta_kl`.
Therefore the environment/history channel has dimension at least four. Under the specific active-plus-transition direct-sum model, the minimal carrier dimension is at least 4+4=8.

However, the gate does not independently construct or test the isometry that forces this Gram matrix:
- it sets `channel_gram = I4.copy()`;
- it then computes the rank of that preselected identity matrix;
- it separately constructs one particular 8×8 permutation `U=[[0,S4],[I4,0]]`.

Thus `reversible_no_link_channel_requires_rank4` is true by construction in the code. The abstract inner-product proof in the markdown supplies the real lower-bound argument; the executable check is only a witness for one realization.

More importantly, the source repeatedly calls the carry cycle **unique**. The lower-bound theorem does not determine the return map from the four history channels to the active states. That is an additional dynamical choice. For example, with active/transition block decomposition, the chosen model fixes active `k -> transition k -> active k+1`; the rank argument alone only requires four orthogonal history labels. Other reversible unitary completions and other channel-return permutations are possible unless the return rule, phase convention and equivalence under history-basis relabelling are specified.

**Correct statement:** an eight-dimensional reversible dilation is minimal in the declared 4-active + rank-one-bottleneck model, and the repository gives an exact order-eight carry-cycle realization. Uniqueness of that dynamics is not established by the channel-rank argument. State whether uniqueness is intended literally, up to unitary equivalence, or only after fixing a particular carry map.

## 16. Q2 winding-character gate has an optional-input fail-open mode

In `scripts/q2_winding_phase_character_gate.py` (blob `0744b95e2d151942519e30a6217fd92fffcf6c3c`), `--history-json` is optional. If omitted:
- `history_input_ok` remains initialized to `True`;
- `history_evidence` remains `None`;
- the script can still return `passed=true` after checking the abstract group-character algebra.

The abstract theorem `Omega(w)=u^w` is correct conditional on an additive integer winding group and multiplicative unit-norm weights. But a gate advertised as validating the **history-to-winding** pipeline should not pass without its upstream history evidence. The phase-refinement workflow does pass `--history-json`, so its intended CI path is not affected; direct CLI use can nevertheless produce a misleading green status.

Required hardening:
1. make `--history-json` required for an end-to-end gate; or split the tool into clearly named `character_algebra_only` and `history_to_character` modes;
2. in end-to-end mode validate the upstream schema, exact status, artifact provenance and required winding fields;
3. include the upstream artifact SHA in the output;
4. add a test asserting missing history input fails closed.

Do not describe a standalone abstract character identity as independent evidence that the repository's physical histories actually generate integer winding.

## 17. Q2 unified-J CI dependency closure is incomplete

The wrapper `scripts/q2_unified_complex_structure_compatibility_gate.py` imports and calls:
- `q2_history_fourier_real_complex_structure_gate.py`;
- `phase_invariant_quadratic_weight_gate.py`;
- `realification_quantum_dynamics_gate.py).

But `.github/workflows/q2-unified-complex-structure.yml` path filters watch only:
- `Q2_UNIFIED_COMPLEX_STRUCTURE_COMPATIBILITY.md`;
- the wrapper script;
- the workflow YAML itself.

The workflow executes the three imported gates, but changes to those imported dependencies do not trigger it automatically. This recreates the same dependency-closure defect previously found in the Peter–Weyl L2 workflow.

There is a related separation: `q2_orientation_complex_conjugation_gate.py` is not invoked by the unified-J workflow, and the phase-refinement workflow does not path-watch or execute that orientation/conjugation gate. Therefore the claimed cross-layer package currently has individually testable pieces but no single dependency-closed CI job that verifies the complete chain from q2 label reflection through orientation pseudoscalar, history Fourier block, winding character, quadratic invariant and realified dynamics.

Required correction:
- add every imported source and its own fixtures/data to path filters;
- run the orientation/conjugation gate in the unified package;
- add a final aggregator that records each upstream artifact SHA and refuses missing evidence;
- keep physical time, Born rule, Hamiltonian and matter interpretation explicitly open.

**Disposition:** the q2 real-complex-structure chain is a mathematically coherent and valuable finite-dimensional representation programme. The exact algebraic identities are not invalidated by these audit findings. What is not yet established is uniqueness of the physical history dynamics or a single independently reproducible CI certificate covering every layer.


## 18. W5 scalar-covariance blocker: pruning ablation is specified, but not executed on the audited branch

Audited ref: `audit/depth6-stageb-ledger-reconcile-20261110`, HEAD `56fa40e4efae47fd600487200f11cf0b0bb80555`.

The W5 pilot's reported forbidden-sector weight is `5.390512613173637e-20`, with allowed-sector weight `6.641081550043974`. The independent zero-tolerance diagnostic reported forbidden-sector weight `4.016710228941321e-30` and total H-sine norm squared `6.641081686796968`. This is a large relative difference in the tiny forbidden tail, but the allowed-sector result also differs by about `1.36752995e-7` (about `2.06e-8` relative). Therefore the issue cannot be dismissed solely as an irrelevant forbidden tail.

The code identifies at least three numerical boundaries that must be separated:
- `project_branch_complete_charges` production tolerance `1e-11`;
- `add_state` production accumulator tolerance `1e-11`;
- an internal basis coefficient cutoff `1e-13` which remains active even when the caller passes `tol=0`.

The ablation design is thoughtful: it compares projection and accumulator tolerances independently, flat versus per-oriented-spec accumulation, reversed spec order, and a `math.fsum` real/imaginary reference on the same generated branches. The newer one-pass script `scripts/w5_scalar_tolerance_ablation_onepass.py` (blob `73e1235999b75a7665ed7d374da28c4e6eddd4fc`) generates each branch once and evaluates both projection tolerances on that same branch, which is preferable to repeated expensive branch regeneration.

However, its workflow `.github/workflows/w5-scalar-tolerance-ablation.yml` is configured with `push.branches: [main]`, and the inspected Actions inventory contains no run named W5 scalar tolerance ablation. The audited branch is 3 commits ahead and 5 behind main, so its version of the workflow/script has not been demonstrated by a branch-specific CI execution. The workflow invokes the one-pass script, not the older multi-pass script with the similar name.

**Disposition:** W5 remains a genuine scalar-covariance blocker. The ablation is a well-designed diagnostic, not a result. Run it against the exact audited branch HEAD (manual dispatch or a temporary branch-aware trigger), preserve the frozen production guard, and inspect all eight full-map comparisons before attributing the residual to numerical noise or changing any physical operator. No five-seed W5 covariance, five-node closure, physical Gram, or joint physical kernel is certified by the evidence reviewed here.


## 19. Active-habitat master: recovered successful artifact contains a strong finite leakage result

The v1 workflow did run successfully:
- run `33304603146`, branch `audit/active-habitat-master-v1`, HEAD `43f8d00ea237e2123ec07d88c1aa2d3d91f683e2`;
- job `active-habitat-master`: success;
- artifact `active-habitat-master-evidence`, ID `9730364258`, downloaded and inspected.

The artifact reports:
- active one-hit dimension: 4 (seed plus three independent one-hit states);
- orthonormality defect: `6.6651e-16`;
- node leakage/full-action operator-norm ratios: node 0 `0.7894673138`, node 1 `0.7894673138`, node 2 `0.7745482093`;
- full master eigenvalues: `[7.0797521, 11.2566282, 11.8042842, 11.8719658]`;
- projected-shortcut eigenvalues: `[1.9892391, 2.3751064, 2.7154066, 7.0797521]`;
- leakage Gram eigenvalues: `[0, 8.9665323, 9.1565593, 9.7300345]`;
- full-versus-projected spectral distortion: `0.6783828012`;
- extra projected zero modes at relative thresholds 1e-8, 1e-10, 1e-12: `[0,0,0]`;
- master/leakage decomposition relative error: `4.7827e-16`.

This is a substantive finite result: the one-hit active habitat is very far from invariant under the three-node Euclidean graph-changing constraints. The projected shortcut substantially changes the master spectrum, but in this particular four-dimensional diagnostic it does **not** create extra apparent zero modes at the three tested relative thresholds. The result is therefore stronger and more precise than a generic statement that projection creates spurious physical states: it demonstrates large positive leakage and a 67.8% spectral distortion, without observing an extra zero mode in this specific sample.

The artifact's `passed=true` is correctly an implementation/decomposition pass; the scientific status is `LEAKAGE_PRESENT_PROJECTED_MASTER_NOT_EQUIVALENT`. It is not a physical-kernel result.

### Numerical/operator boundary

The source implementation applies `PW.prune_state(..., 1e-8)` after every `H_v) action and its sparse helper routines also discard amplitudes below `1e-13`. Thus the reported matrices are built from a thresholded finite action. The identity `M_full-M_proj=M_leak` is checked to high precision for those generated columns, but that identity does not independently prove that the thresholded action equals the unpruned linear Hamiltonian on arbitrary superpositions. Add a no-prune reference and a complex-superposition linearity test before promoting the result to an operator-domain statement.

### V2 preregistration status

`ACTIVE_HABITAT_MASTER_GAP_FLOW_V2_PREREGISTRATION.md` freezes a deeper P0/P1/P2 nested-cone experiment before reading v1's result. That is good experimental discipline. But in the audited v2-preregistered branch:
- no `active_habitat_master_gap_flow_v2_gate.py` exists;
- no depth-two implementation exists;
- no branch-specific Actions run exists;
- the inherited workflow push trigger is restricted to `audit/active-habitat-master-v1`, not the v2-preregistered branch.

Therefore v2 is a genuine preregistration, not an executed experiment. Its next step is implementation against the frozen P2 definition, without changing its 1e-11 dependence rule or the five preregistered zero-mode thresholds after seeing the v1 spectrum.


## 20. Five-step representation-RG ancestry: one-step certificates are real, global composition remains untested

The five consecutive one-step channel-flow matrices are backed by a successful historical workflow run `42797583866` and persisted JSON shards for `s2=3,4,5,6,7`. Their row-stochasticity and paired-eigenvalue checks are valid finite representation-RG evidence. The maximum-product ancestry path is also correctly labeled as a local scalar-overlap candidate.

The direct two-step composition script `scripts/bqg_two_step_channel_composition.py` (blob `22f7357810329f2df44ab815ea3c87b465553828`) is a meaningful next experiment. It recomputes both one-step matrices, checks them against the saved references, propagates the actual intermediate rank-2 carrier, tests invariance under a U(2) basis rotation, and compares:
- direct composed overlap `G_direct`;
- scalar matrix product `F_first @ F_second`;
- maximum-product intermediate-channel path;
- intermediate off-diagonal block norms and within-channel anisotropy.

No run of `bqg-two-step-channel-composition` was found in the inspected Actions inventory. Its workflow push trigger is restricted to `main`; the audited branch `audit/depth6-stageb-ledger-reconcile-20261110` has 3 commits ahead and 5 behind main. The workflow supports manual dispatch, but no manual branch run or result artifact was found.

This means the five-step maximum-product path remains a **local ancestry statistic**, not yet a demonstrated global microscopic refinement lineage. The missing comparison is not cosmetic: scalar multiplication `F_j F_{j+1}` is equivalent to an intermediate channel measure-and-prepare map and can discard coherence. The direct `G_direct` computation is precisely the experiment needed to test whether that information matters.

**Disposition:** preserve all five one-step matrices and their historical successful run. Do not promote their product to a global refinement fidelity or physical transition probability. Execute both source pairs `s2=3` and `s2=4` at the audited branch HEAD, retain all artifacts and compare direct composition defects before making any multi-scale physical interpretation.


## 21. Cross-branch analytical bridge: exact finite identities and conditional continuum chain

This section joins the independently audited Q2, active-habitat, K5-volume, K5-HDA, representation-RG and binary-dimension branches. It is an analytical dependency map, not a declaration that the complete physical theory is closed.

### 21.1 Binary selector and the causal-volume fixed point

The frozen local selector is
[
q+2=2^q,qquad qinmathbb N, qge1.
]
Its unique integer solution is (q=2): (q=1) fails, (q=2) satisfies it, and for (qge2) the increment (2^q-(q+2)) is strictly increasing after its zero at 2.

The frozen route rewrite sets (B=2^q=4). With (E_0=1), active-edge growth and graph-distance growth are
[
E_g=(2B)^g=8^g,qquad L_g=2^g.
]
The exact vertex count is
[
N_g=2+Bsum_{r=0}^{g-1}(2B)^r
=2+4rac{8^g-1}{7}
=rac{4,8^g+10}{7}.
]
Hence
[
d_g=rac{log(N_g/N_{g-1})}{log 2}
=3+log_2!left(1-rac{35}{16,8^{g-1}+40}ight).
]
For (gge2), the correction is negative and increases monotonically to zero, so (d_g\nearrow3). This is an exact asymptotic exponent for the frozen route rewrite. It does not prove that every dynamically generated graph or quantum ensemble obeys that rewrite.

### 21.2 Independent rotation/HDA dimension estimator

For a q-level local Hilbert space, the traceless Hermitian observable space has real dimension (D_q=q^2-1), and (dim PSU(q)=D_q). If the physical interpretation requires the available local frame action to cover all of (SO(D_q)), dimension counting requires
[
D_q=rac{D_q(D_q-1)}2.
]
For (D_q>0), this gives (D_q=3), hence (q=2). For q=2, the adjoint map (SU(2)/\mathbb Z_2\to SO(3)) is indeed onto. For q>2, dimension equality fails.

This is a conditional isotropy selector, not a derivation that the microscopic dynamics realizes full isotropy. Under standard D-dimensional HDA/ADM assumptions, the inverse DeWitt trace coefficient is (c_D=1/(D-1)), so
[
D_{HDA}=1+c_{eff}^{-1}.
]
Agreement with the binary rotation selector requires (c_{eff}\to1/2). This is an independent cross-check only if (c_{eff}) is extracted from the microscopic kinetic form without fitting it to D=3.

### 21.3 Exact K5 local volume algebra

On three distinct link factors, let
[
Q_{123}=\epsilon_{ijk}J_1^i\otimes J_2^j\otimes J_3^k.
]
With Hermitian (J_r^i) and distinct tensor factors commuting,
[
Q_{123}^\dagger
=\epsilon_{ijk}J_1^i\otimes J_2^j\otimes J_3^k
=Q_{123}.
]
For a leg permutation \(\pi\), with permutation unitary (U_\pi),
[
Q_{\pi(1)\pi(2)\pi(3)}
=\operatorname{sgn}(\pi)U_\pi Q_{123}U_\pi^\dagger.
]
Define (V_{123}=\sqrt{|Q_{123}|}) by spectral calculus. Then
[
V_{123}^\dagger=V_{123}\succeq0,\qquad
V_{123}^2=|Q_{123}|,\qquad
V_{\pi}=U_\pi V_{123}U_\pi^\dagger.
]
The last identity follows because (|-A|=|A|) for Hermitian A and functional calculus is unitary covariant. Thus orientation-odd Q and orientation-even derived volume are compatible exactly. The production code's pre-eigensolver symmetrization is only defensible as roundoff cleanup after the raw anti-Hermitian defect is independently shown negligible; the finite tests now record this defect before symmetrization.

At a K5 vertex, the oriented combinatorics has four choices of omitted neighbour and three cyclic orders, hence 12 distinct signed specifications, with sign ((-1)^r) for omitted-neighbour position r. This proves the finite orientation-list combinatorics; it does not prove adjointness of the complete graph-changing T sum.

### 21.4 Exact adjoint law for the graph-changing monomials

For a primitive sequence (T=A_1A_2\cdots A_n) on a common finite inner-product space,
[
T^\dagger=A_n^\dagger\cdots A_2^\dagger A_1^\dagger.
]
The K5 production reverse sequence must therefore reverse path order, exchange input/output magnetic labels, and retain the self-adjoint volume insertions. For the complete signed 32-term sum (T=\sum_{\alpha=1}^{32}s_\alpha T_\alpha),
[
T^\dagger=\sum_{\alpha=1}^{32}s_\alpha T_\alpha^\dagger
]
because the signs are real. This identity is exact; a single matrix element or selected 2×2 block does not establish it on the entire finite sector. The present source has a selected 2×2 coherent block, projection-leakage ledger and a spin-wall 2/4/6 experiment. The fixed-input ceiling argument gives maximum doubled spin 3 for that specific T input, so walls 4 and 6 should agree for that input. It is not a universal wall theorem for arbitrary inputs or HH.

### 21.5 Full-versus-projected master: exact leakage identity

Let P be an orthogonal projector, (H_v) node constraints, and (L_v=(I-P)H_vP). Then
[
M_{full}=\sum_v(PH_v^\dagger H_vP),\qquad
M_{proj}=\sum_v(PH_vP)^\dagger(PH_vP)
]
satisfy the exact identity
[
M_{full}-M_{proj}
=\sum_vPH_v^\dagger(I-P)H_vP
=\sum_vL_v^\dagger L_v\succeq0.
]
This is the correct operator-level reason that premature compression can change the master spectrum. The v1 active-habitat artifact gives a finite example with large leakage and 67.8% spectral distortion; it does not establish the unrestricted physical kernel. The identity should be evaluated with unpruned coherent columns and off-diagonal interference retained.

### 21.6 K5 fixed-sector SO(5) skeleton is a projection, not full HDA closure

The oriented vector-5 (J_{max}=1/2) diagnostic reports for all ten node pairs:
[
\|[H_i,H_j]P_{full}\|=3906.6220702801,
\qquad
\epsilon_{graph}=\sqrt{37/69},
\]
and fixed-sector norm fractions (32/69) and graph-changing fraction (37/69). The projected generators (Q_{ij}=-iP_{full}[H_i,H_j]P_{full}) span a ten-dimensional Frobenius-orthogonal space. Their projected overlapping brackets reproduce the \(so(5)\) structure constants with normalization 32, but only (1/297) of the full bracket norm lies in that projected span; the orthogonal fraction is (296/297).

Therefore the SO(5) skeleton is a real exact/finite projected algebraic pattern, but it is not the full graph-changing constraint algebra. A physical HDA claim requires the full commutator to match the diffeomorphism/structure-function target plus Gauss terms, including graph-changing channels, in a regulator-safe common domain.

### 21.7 The conditional binary-HDA chain and its unproved hypotheses

If (i) q=2 is the full local spatial frame carrier, (ii) the microscopic dynamics restores full SO(3) isotropy, (iii) a local two-derivative metric Hamiltonian emerges, and (iv) the standard first-class HDA is recovered, then
[
D_{rot}=3,quad c_{eff}=1/2,quad A_{eff}B_{eff}=1,quad z=1,
quad N_{TT}=2.
]
The TT equation is (ddot h_{TT}=A_{eff}B_{eff}\nabla^2h_{TT}), so (A_{eff}B_{eff}=1) gives a unit relativistic tensor cone in the chosen lapse/time normalization. In D spatial dimensions, first-class spatial diffeomorphism plus Hamiltonian constraints leave ((D+1)(D-2)/2) configuration degrees of freedom; at D=3 this is two.

The chain is conditional because the repo has not yet established the microscopic hypotheses (i)-(iv) together. In particular, the finite K5 HDA diagnostic has substantial graph leakage and the regulator-safe full HH/commutator comparison remains outstanding.

### 21.8 What must be proved before calling this a complete analytic theory

The remaining proof obligations are logically distinct:
1. **Microscopic operator definition:** one frozen no-prune action, common dense domain, inner product and exact adjoint convention for every primitive.
2. **Full constraint algebra:** complete graph-changing ([H(N),H(M)]), not only fixed-sector projections, equals the discrete diffeomorphism generator with the derived structure function plus Gauss terms.
3. **Regulator control:** a proven per-history reachability bound or a controlled increasing-wall limit; finite wall equality for one input is not enough.
4. **Physical projector:** the full leakage Gram and kernel must be computed without replacing (M_{full}) by a compressed shortcut.
5. **RG composition:** direct microscopic refinement overlaps (G) must be compared to scalar (F_1F_2); the latter is an intermediate measure-and-prepare channel unless a factorization theorem is proved.
6. **Lorentzian completion:** KKV adjoints must be reversed in operator order; partial auxiliary trace cyclicity is unavailable for operator-valued entries. Physical ordering and normalization must be independently fixed.
7. **Continuum/IR theorem:** establish convergence of embedded operators, projectors, structure functions and observables, not only convergence of one finite exponent or one spectral gap.
8. **Empirical map:** derive dimensionful predictions with fixed constants and uncertainty bounds, then test against independent observations.

**Status:** exact finite identities have been connected into one dependency graph. The full analytical/physical closure is not yet proved. No theorem above is promoted beyond its stated hypotheses.


## 22. Q2 orientation/history/metric/scalar source: exact identities and carrier obstructions

Audited source branches:
- `research/q2-orientation-intrinsic-metric-nogo`;
- `research/q2-history-directed-laplacian-factorization`;
- `research/q2-oriented-flux-history-observable`;
- `research/q2-relational-metric-source`;
- `research/first-q2-scalar-effective-action`.

### 22.1 Orientation is a genuine gauge-scalar flux observable, but not a linear intrinsic-metric coordinate

For four spin-1/2 face fluxes, the oriented triple product
[
Q_{or}=\epsilon_{abc}J_1^aJ_2^bJ_3^c
]
restricted to the two-dimensional four-valent singlet carrier obeys the exact identity
[
Q_{or}^{logical}=\frac{\sqrt3}{4}Y_L,
\qquad
Y_L=\frac4{\sqrt3}Q_{or}.
]
It commutes with all total SU(2) generators and changes sign under odd face permutations. Thus (Y_L) is not merely a label in this frozen carrier.

At the same time, the exact logical-to-intrinsic-metric Jacobian has columns
[
B_{XYZ}=(\mathrm{vec}M_X,0,\mathrm{vec}M_Z),
qquad
B_{XYZ}^TB_{XYZ}=\operatorname{diag}(9/2,0,9/2).
]
Consequently (\partial g/\partial Y=0), and for any source covariance (\Sigma),
[
C_{metric}=B_{XYZ}\Sigma B_{XYZ}^T
]
has rank at most two. Even (Sigma=I_3) cannot restore the missing orientation direction. There is no contradiction: the flux pseudoscalar is orientation-sensitive, while the intrinsic metric is orientation-even. A nonzero orientation-odd Lorentzian response must be sought in flux/triad, connection, extrinsic-curvature or history-current observables, not forced into a linear parity-even metric source.

### 22.2 Complete directed history difference exactly factors the undirected C8 Laplacian

For the eight-state reversible history shift (U_8) and (P_\pm=(I\pm Y_L)/2), define
[
W=P_+\otimes U_8+P_-\otimes U_8^{-1}.
]
Then (W) is unitary and
[
W+W^\dagger=I_{geom}\otimes(U_8+U_8^\dagger).
]
For the full directed difference (\Delta_W=W-I),
[
\Delta_W^\dagger\Delta_W
=2I-W-W^\dagger
=I_{geom}\otimes L_h,
quad
L_h=2I-U_8-U_8^\dagger.
]
Since W is normal, (\Delta_W\Delta_W^\dagger) gives the same operator. On (U_8|\theta\rangle=e^{i\theta}|\theta\rangle),
[
\Delta(\theta)=(\cos\theta-1)I-\sin\theta J,
\quad J^2=-I,
]
and
[
\Delta(\theta)^\dagger\Delta(\theta)
=4\sin^2(\theta/2)I=2-2\cos\theta.
]
The odd Hermitian current alone has square (\sin^2\theta=L_h-L_h^2/4) and therefore has the extra finite-lattice zero at (\theta=\pi). The complete forward difference has only the trivial C8 zero. This is an exact finite graph identity, not yet a physical Dirac operator, physical time evolution or fermion construction.

### 22.3 The relational projector/source construction is exact as a positive control

For the C8 control step (R=J), define
[
V|\psi\rangle=\frac1{\sqrt8}\sum_{t=0}^7|t\rangle\otimes R^t|\psi\rangle.
]
Then (V^\dagger V=I) and (VV^\dagger=P_{rel}). For any geometry operator O,
[
O_{rel}=\sum_t|t\rangle\langle t|\otimes R^tOR^{-t},
\quad O_{rel}V=VO,
\quad V^\dagger O_{rel}V=O.
]
For the two shape sources (K=j_XX+j_ZZ), with (K^2=r^2I), (r^2=j_X^2+j_Z^2),
[
Z(j)=\tfrac12\operatorname{Tr}(e^K)=\cosh r,
\qquad W(j)=\log\cosh r.
]
The mean-field magnitude is (s=\tanh r). The exact radial Legendre transform is
[
\Gamma_{shape}(s)=s\operatorname{artanh}s+\tfrac12\log(1-s^2),
\quad 0\le s<1,
]
with
[
\Gamma'(s)=\operatorname{artanh}s,
\quad
\Gamma''(s)=\frac1{1-s^2}>0.
]
Its small-s expansion is (s^2/2+s^4/12+s^6/30+s^8/56+O(s^{10})), and the limit at (s\to1^-) is \(\log2\). This is an exact local nonlinear 1PI functional of the registered finite control.

However (R=J) is the declared positive-control step, not the actual graph-changing gravitational history amplitude. Therefore this \(\Gamma\) is not a physical cosmological scalar action or graviton effective action.

### 22.4 Exact scalar-carrier obstruction

The local q=2 metric tangents (M_X,M_Z) are tracefree with respect to the regular background (g_0):
[
\operatorname{Tr}(g_0^{-1}M_X)
=\operatorname{Tr}(g_0^{-1}M_Z)=0.
]
The conformal direction is orthogonal to this rank-two tangent. The present source has no independent lapse response and no nonfactorizing inter-block history amplitude. If its only justified multi-block extension is (Z_N=\prod_b Z_b), then
[
\partial^2 W_N/(\partial j_{A,b}\partial j_{B,c})=0
\quad(b\ne c).
]
Thus no momentum-dependent scalar Poisson kernel or cosmological potentials follow from the local control. The next structurally motivated carrier is the collective radial/volume channel: fixed (j=1/2) freezes the absolute four-valent volume, while the equal-spin (j=1) sector is the first where it is non-scalar. The symmetric-block (2\to j=1) result is conditional on symmetric blocking; the microscopic dynamics has not yet been shown to select it.

**Disposition:** the q=2 chain now has a precise split: orientation flux and directed history algebra are exact finite observables; the intrinsic metric deliberately forgets orientation; the source Legendre transform is exact but belongs to a positive control; physical Lorentzian history, conformal/lapse response and connected inter-block scalar dynamics remain open.


## 23. Whole-repository branch census and scope correction

The refreshed GitHub branch census at this checkpoint contains **168 branch refs and 157 unique HEAD commits**. The family inventory includes 95 `research/*` branches, 22 `audit/*` branches, 10 `exp/*` branches, 3 `docs/*` branches, and the remaining implementation, recovery, release, CI and maintenance branches. Duplicate branch names sharing the same HEAD are retained in the branch inventory but count as one source snapshot for content analysis.

Important scope correction: the previously committed branch-wide Python AST audit (`PYTHON_BRANCH_WIDE_AUDIT_2026-11-10.md`) covered an earlier snapshot of **157 branches / 146 unique HEAD commits**, not the refreshed 168 / 157 census. It reported 18,674 Python-file occurrences across branch snapshots, 674 unique Python paths, 757 path/content versions, 754 unique source blobs, 752 parsed sources and 2 historical syntax errors. Those are valid historical audit results for that earlier census, but they must not be presented as a complete AST scan of the refreshed 168/157 snapshot.

Therefore the current full-repository work has two separate scopes:
- **Branch census:** refreshed to 168 refs / 157 unique source snapshots.
- **Deep analytical cross-branch review:** currently covers the linked Q2 orientation/history/source chain, active-habitat leakage identity, K5 volume and graph-changing adjoint obligations, finite HDA obstruction, RG composition caveat, binary dimension selector and conditional HDA uniqueness chain.

The full all-branch content/AST refresh and exhaustive symbolic theorem extraction over all 157 distinct snapshots remain unfinished. The next repository-wide pass must first scan the 11 additional unique HEAD snapshots relative to the previous 146-head AST baseline, deduplicate blobs, then classify every distinct theorem/specification and map exact branch ancestry. It must not silently infer that a branch is covered merely because its name appears in the census.

### Analytical proof ledger format required for every branch family

For each candidate theorem, record:
1. exact source definition and branch/commit/blob;
2. mathematical objects, domains, inner products and gauge quotients;
3. assumptions and quantifiers;
4. line-by-line symbolic derivation;
5. equivalence transformations and conditions for each;
6. counterexamples/no-go results from sibling branches;
7. executable finite gate and whether it was actually run;
8. cutoff/pruning/roundoff error statement;
9. refinement/continuum theorem or explicit open obligation;
10. final status: exact finite identity, conditional theorem, finite numerical evidence, conjecture, falsified ansatz, or untested.

A result is promoted to the common theory only when its hypotheses are compatible with upstream and downstream branches. A formula with an unproved bridge is retained as a conditional node, not merged into a stronger theorem.


## 24. Five-block microscopic P-carrier -> signed spatial TT precursor: a real new bridge, with a precise boundary

New branch deeply inspected: `physical-tt-kernel-leading-ir`, HEAD `e3118e1a710b83d23622af57d9c41ed66d7da2c4`.

Primary source files:
- `BQG_MULTIBLOCK_TT_WILSON_PRODUCTION_STATUS.md`, blob `c640d361ecc63c6b596e0c9f9d8be190d8ba3df7`;
- `scripts/bqg_fiveblock_p_carrier_producer.py`, blob `45eae813787e535a40cd9ff1a49dc30e8159e025`;
- `.github/workflows/bqg-signed-five-block-bridge.yml`, blob `f8a29c24992b0b457265595051b4f71ee2361573`;
- `scripts/bqg_multiblock_tt_wilson_extractor.py`, blob `b74ddcc41a525aa561660fb7a62d55835c503032`;
- `PHYSICAL_TT_KERNEL_LEADING_IR_DERIVATION.md`, blob `00f9e8b1a046fd795a6d07d2f70d263c5db69e3b`.

### 24.1. Frozen microscopic object and carrier

The signed operator is frozen as
[
G_{frozen}=-\frac23 H_E^{sine}-\frac{32}{9}S,
\qquad
S=-\frac i2(L_{raw}-L_{raw}^{\dagger}).
]
No (R_{op}), GR ratio, TT target or fitted Wilson coefficient enters the P-carrier producer.

The carrier is constructed on the full L1 16-cell habitat (384 fine tetrahedra, 768 dual links). Five coarse parent tetrahedra are selected: central parent 0 and its four face-neighbours 1,2,4,8. Each parent supplies 24 strict-interior q4 source columns, combined through the frozen S4 parity/coset rule into six edge tangents:
[
P\in\mathbb C^{960\times30},\qquad
P^\dagger P\simeq I_{30}.
]
The branch report records rank 30, relative Gram defect (4.5324665183683945\times10^{-17}), common unnormalised six-edge norm square (0.6982713481539777), cross-parent maximum overlap 0, and zero cross-parent basis-label collisions. Every normalized column has support 32. The serialized NPZ/direct-Gram relative defect is (1.7786683001560746\times10^{-16}).

The geometrical reason for no recoupling at this stage is precise: each strict-interior P column changes only links internal to its own parent; all parent-boundary links and exterior labels remain at background. Thus distinct parent supports are disjoint. Full shared-face SU(2) recoupling belongs after (H_E^{sine}) or (S) acts, because those images can touch/cross shared faces. Applying shared-face recoupling prematurely to P would change the carrier being defined.

### 24.2. Exactness qualification discovered by source inspection

The producer is finite and explicit, but the word “exact” must be qualified at the amplitude level. Its `add()` routine discards amplitudes with (|z|\le 10^{-10}) (`TOL=1e-10`), and the q4 source dictionary is filtered again by that threshold. Therefore the recorded P is the output of a fixed-threshold sparse representation, not a proof that every mathematically nonzero amplitude in the untruncated microscopic state has been retained.

A strict follow-up must run a cutoff ladder (\epsilon=10^{-8},10^{-10},10^{-12},10^{-14}) (and, where feasible, a no-prune reference), then compare:
[
\|P_\epsilon-P_{\epsilon'}\|,\quad
\sigma_{min}(P_\epsilon),\quad
\|P_\epsilon^\dagger P_\epsilon-I\|,
\quad
\text{support and boundary-label invariants}.
]
The source also symmetrizes each Gram before reporting its Hermiticity defect; consequently the reported zero Hermiticity defect is not an independent raw-adjoint diagnostic. Future gate should record (|K_{raw}-K_{raw}^\dagger|/\|K_{raw}\|) before symmetrization, for both local six-edge and global 30-column Grams. These are audit hardening requirements, not evidence that the current rank/Gram result is false.

### 24.3. The actual next operator calculation is still missing

After constructing P, the next required calculation is:
[
P\xrightarrow{\ H_E^{sine},\ S\ }\{H_E^{sine}P,SP\}
\xrightarrow{\text{complete shared-face recoupling}}\mathcal H_{E/S}
\xrightarrow{\text{remove }\operatorname{span}P,\;\text{whiten residual}}Q_Q,
\qquad W=[P,Q_Q].
]
Then compute (E_W=W^\dagger H_E^{sine}W), (S_W=W^\dagger SW), and the leakage ratios
[
\eta_O=\frac{\|(I-WW^\dagger)OW\|}{\|OW\|},
\quad O\in\{H_E^{sine},S\}.
]
Only after closure of this target-independent amplitude basis may one assemble
[
G_W=-\frac23E_W-\frac{32}{9}S_W
=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix}
]
and form the orthonormal Schur complement
[
C_{eff}^{orth}=A-BD^+B^\dagger.
]
A zero/gapless Q mode coupled to P must be promoted, not hidden by a fitted (i\eta). Then return to labelled q coordinates, subtract the normalized-state background term, transport by the measured q-to-metric map (M_{hq}), build the five-block spatial symbol and TT-project it. The only permitted output at that stage is (\mathbf c_{micro,spatial}^{BQG}), never (\mathbf c_{BQG}^{IR}).

The workflow file contains a pull-request path trigger and manual dispatch, and has explicit guards for the 61 face patterns, frozen upstream SHA, P rank, Gram and disjoint-support invariants. In this pass I inspected the workflow definition but did not independently retrieve a matching Actions run/artifact; do not silently turn a workflow specification into an executed CI result.

### 24.4. The physical TT leading form is conditional, not yet BQG output

The separate ADM/HDA derivation gives, inside its registered local two-derivative family,
[
c=\frac12,\qquad AB=1,
]
and for two TT polarizations
[
S_{TT}^{(2)}=\frac1{4A}\sum_{A=1}^{2}
\int d^4x\,[\dot h_A^2-(\nabla h_A)^2].
]
Hence the *conditional ADM/HDA reference kernel* is
[
K_{E,TT}^{ADM}=Z_T(\omega_E^2+k^2)I_2+O(\partial^4),
\quad
K_{TT}^{ADM}=Z_T[-(\omega+i0)^2+k^2]I_2+O(\partial^4).
]
This yields (m_{TT}=0), two polarizations and (c_T=1) under those premises. It does not compute the BQG physical projector, history measure, connected (\Gamma_{phys}), physical residue/Newton normalization, or physical Wilson vector. The source explicitly keeps `physical_TT_kernel_frozen = false`. Also (z_{constraint}\ne\omega_{physical}) absent an independent clock/history map.

## 25. Physical TT and scalar channels: complete finite symmetry spaces, then physicalization

### 25.1. Six-dimensional quartic TT quotient is the general finite-regulator target

The corrected onsite/momentum distinction is essential. At zero momentum (or a fully stabilised direction), the S4 onsite six-edge kernel has the three-orbit form (aI+bA+cO). At generic nonzero momentum only covariance is required:
[
C_6(\omega,g\mathbf k)=U_gC_6(\omega,\mathbf k)U_g^{-1};
]
one must not impose fixed-direction invariance and collapse the answer to three scalar functions.

For the parity-even quartic metric response, the unrestricted representation count is 13 singlets before TT reduction. Imposing (\operatorname{tr}h=0) and (h_{ij}k_j=0) leaves exactly
[
\dim\mathcal W^{(4)}_{TT,S_4}=6.
]
The six Reynolds invariants (W_1,\ldots,W_6) and exact six-observable extraction matrix are recorded in `S4_TT_QUARTIC_COMPLETE_BASIS.md`. Three high-symmetry directions (100),(110),(111) have rank five only; the preregistered generic direction (120) closes the system. The exact extraction determinant is (1/699840000\ne0). Thus the production microscopic result must first be a full six-vector, not a two-number ansatz.

### 25.2. The six -> one -> zero hierarchy is a conditional theorem

The physical-pole branch contains a useful nested theorem:
[
\mathcal W^{(4)}_{TT,S_4}:6
\supset
\mathcal W^{(4)}_{TT,SO(3)}:1
\supset
\mathcal W^{(4)}_{TT,\text{local Lorentz metric massless pole}}:0.
]
First, spatial SO(3) invariance forces a parity-even real symmetric kernel on the two TT polarizations to be scalar, leaving only (eta_s a_*^2k^4I_{TT}). Second, if the physical vacuum is locally Lorentz invariant, the only gapless spin-2 field is the metric graviton, and the quadratic 1PI action has a local analytic derivative expansion, then
[
K_{TT}(s)=sI_{TT}F(a_*^2s),\qquad
s=-\omega^2+c^2k^2,\quad F(0)=Z_T>0.
]
The original massless branch remains (s=0) at every finite local derivative order. Therefore a nonzero physical massless-pole k4 shift requires a real failure of at least one premise (preferred foliation/tensor order, extra gapless mixing, Lorentz breaking, nonlocality, or non-Poincare background). A regulator orientation alone is not such a physical evasion.

This theorem is conditional EFT structure, not a result that BQG already satisfies the assumptions. It changes the scientific decision tree: extract all six finite-regulator coefficients blindly; test their SO(3) collapse; then test whether the isotropic term survives on the Lorentz-invariant massless pole. The bare lattice values are UV controls, not automatically IR observables.

### 25.3. Scalar channel: local exact control is not a propagating physical scalar

The q=2 relational source has an exact finite local Legendre transform, but its carrier is only a rank-two tracefree shape slice. Its conformal/volume direction and independent lapse response are absent; a product extension factorizes across blocks and has exactly zero connected cross-block Hessian. Hence this source alone cannot yield (\Phi,\Psi,\rho_{hist}(a)), a scalar (k)-dependent kernel, or a dark-sector prediction.

The correct physical scalar object is the connected relational history kernel
[
C_{AB}(\Delta t;b,c)
=\langle\delta\mathcal O_A(t,b)
\delta\mathcal O_B(t+\Delta t,c)\rangle_{phys,c},
]
followed by a justified relational-time transform and inversion on physical support:
[
\Gamma_{scalar}^{(2)}(\omega,k)=C_{phys}^{-1}(\omega,k).
]
The six S3 transfer functions of the shared-face pair can constrain the tensor form of the nearest-neighbour response, but their physical-history weighting and values remain uncomputed. The three variables (	au) (master heat depth), (z) (constraint resolvent) and (\Delta t_{rel}) (physical relational separation) must remain distinct; neither (	au) nor (z) is physical frequency.

### 25.4. Immediate proof obligations

1. Reproduce the five-block P carrier under a cutoff ladder and raw-adjoint diagnostics.
2. Execute the target-independent E/S image producer with full shared-face recoupling and no target pruning; publish wall-closure and leakage bounds.
3. Independently reconstruct the Schur complement from the full compressed block matrix and verify the gapless-Q promotion rule.
4. Freeze the six-Wilson extractor and verify its rational matrix/inverse independently (symbolic exact arithmetic), then compare direct directional observables against reconstruction.
5. Derive the actual physical projector/rigging map and relational clock; demonstrate the connected metric/scalar history kernel and its regulator/refinement limit.
6. Test the constraint algebra through the same derivative order as the physical TT/scalar kernel; only then identify the physical pole and observables.

Classification: exact finite representation/quotient identities; conditional EFT theorem; finite numerical carrier evidence at fixed threshold; unexecuted downstream physicalization. No complete analytic proof of the full BQG theory is claimed.


## 26. Independent symbolic verification and checkpoint after this tranche

The six-observable TT extraction matrix printed in section 25 was reconstructed independently in exact rational arithmetic (SymPy, not floating-point linear algebra). The independent calculation returns
[
\det A=\frac1{699840000}
]
and reproduces the exact rational inverse printed in `S4_TT_QUARTIC_COMPLETE_BASIS.md`. Thus the six-direction reconstruction is algebraically nonsingular; this verifies the extraction map, not the physical existence or measured values of the six Wilson coefficients.

Checkpoint:
- consolidated research branch after this report update: `7c00e711a3730150d85d6e9bd22f71787ac3988a`;
- inventory refresh commit: `9d5512336439a5d4fb54982dd79d850efed9bb96`;
- main HEAD remains `254f19fdab41d32b920bf454f987aa7c0fdba34d`;
- current compare main -> consolidated research: **146 ahead / 2 behind**, diverged;
- main was not modified or merged.

The branch manifest remains 168 refs / 157 unique HEADs. This tranche deep-inspected selected physical TT, scalar-history, HDA and operator-micro-dynamics branches; it is not a completed semantic audit or fresh AST scan of all 157 unique heads. The next pass must continue by deduplicating the remaining unique commit trees and blobs, then extract each candidate theorem together with its executable source, test fixture, cutoff policy, workflow dependency closure and actual run/artifact status.


## 27. CI evidence recovery and independent local recomputation (2026-11-11 audit pass)

### 27.1. Actual Actions status on the frozen physical TT branch

The inspected branch `physical-tt-kernel-leading-ir` is at `e3118e1a710b83d23622af57d9c41ed66d7da2c4`. GitHub Actions associates three completed successful runs with that exact commit:

| Workflow | Run | Result |
|---|---:|---|
| `core-regression` | `34542820459` | success; 4/4 jobs |
| `BQG frozen PL Hermitian S adapter` | `34542820360` | success; 1/1 job |
| `bqg-signed-five-block-bridge` | `34542820317` | success; 2/2 jobs |

The `core-regression` jobs were `canonical-core`, `higher-shell`, `l1-q4-metric`, and `directional-regge`. Their completed steps include canonical truth-ledger validation, Python compilation, q=2 dimension/manifold/qubit gluing, graph-link representation, Plebanski/Urbantke reconstruction, 32D Peter-Weyl normalization, two-/three-node HDA controls, Lorentzian support bounds, held-out Regge L6, TT propagator/vacuum, exact six-dimensional S4 TT quotient, Feshbach identity, nearest-block S3, six-Wilson predictor and scale-translator algebra. All listed steps completed successfully.

The five-block run completed both `fiveblock-p-carrier` and `bridge-selftests`. The former reran all 61 complete six-link face-recoupling patterns, constructed the 30-column carrier, and passed all hard guards. The latter passed both signed-operator assembler and Schur/metric/TT extractor selftests. The adapter preflight separately passed frozen upstream SHA, direct-Hermitian identity and 24-term bookkeeping.

These are successful historical runs attached to the exact frozen branch commit (artifact timestamps 2026-09-10), not a claim that a fresh workflow was dispatched on 2026-11-11. The consolidated audit branch itself has no associated PR-triggered workflow run.

### 27.2. Downloaded artifact: independent five-block recomputation

Downloaded artifact ID `10178025209`, name `bqg-fiveblock-p-carrier`, from run `34542820317`. Its SHA-256 for the serialized NPZ is
[
	exttt{f6aa57478dafa952fe8b7cddd26eaa2ea41b4d74bbb9f4a50f50d1363c79a637}.
]
The NPZ contains `basis_spins2: (960,768)`, `basis_K2: (960,384)`, and `P_real/P_imag: (960,30)`. I reconstructed the complex P matrix directly from the downloaded arrays, without calling the producer, and recomputed:
[
\operatorname{rank}(P^\dagger P)=30,
\qquad
\frac{\|P^\dagger P-I\|_F}{\|P^\dagger P\|_F}
=4.4269925586208663\times10^{-16}.
]
The raw reconstructed Gram anti-Hermitian relative norm is (1.0802605865622026\times10^{-34}); the eigenvalues round to ([1,1]), and maximum column norm defect is (3.3306690738754696\times10^{-16}).

A stronger structural property is visible in the serialized artifact: every one of the 960 global basis rows has exactly one nonzero entry across all 30 columns, and every column has exactly 32 nonzero entries. Thus the 30 *stored* columns have pairwise disjoint support, not merely disjoint support between different parent blocks. The smallest stored nonzero amplitude magnitude is (0.02723397555824756), far above the producer's (10^{-10}) retention threshold. This does **not** prove that amplitudes discarded before serialization were exactly zero; it proves that all retained serialized amplitudes are well separated from the threshold.

The artifact's own JSON records `P_Gram_I30_relative_defect=1.812986607347358e-16`, `serialized_P_Gram_relative_defect=2.771854347690155e-16`, rank 30, zero cross-block overlap/collisions and all boundary/exterior checks true. The independent recomputation differs slightly at machine-roundoff scale, as expected from sparse-vs-dense accumulation order.

### 27.3. Full six-link recoupling artifact

The same artifact includes the prerequisite `COLLECTIVE_FACE_COMPLETE_RECOUPLING.json`:
- all 61 ordered q4 face patterns tested;
- recoupling dimensions 16 through 144;
- total doubled-spin support ({0,2,4,6,8});
- worst left/right unitarity defect (4.637311052790768\times10^{-15});
- historical maximal-(J_3) Dicke subblock defect (1.1102230246251565\times10^{-16});
- all five hard checks true.

The current Markdown production-status page reports nearby but non-identical values for the P Gram and worst recoupling defect ((4.5324665\times10^{-17}) and (5.1388141\times10^{-15}), respectively). These are all machine-precision-scale results, but they are not the values in the downloaded run-34542820317 JSON. Treat the run-specific JSON/NPZ and hashes as canonical for that run; do not silently combine values from different production snapshots. The status page should be annotated with run IDs and artifact hashes in a future provenance cleanup.

### 27.4. Independent local reproduction of finite physicalization controls

The `physicalization-truth.yml` workflow is present on the inspected branch, but its triggers are restricted to pushes on `main` and pull requests targeting `main`; it does not run merely because this research branch was updated. I therefore reproduced the mathematical controls independently from the inspected definitions rather than reporting an unavailable branch CI run.

**Finite master-constraint projector.** For stacked constraints (mathcal C:mathcal H\to\mathbb C^m\otimesmathcal H) and a positive-definite matrix (G>0), the exact identity is
[
M_G=\mathcal C^\dagger(G\otimes I)\mathcal C,
\qquad
\langle\psi|M_G|\psi\rangle
=\|(G^{1/2}\otimes I)\mathcal C\psi\|^2.
]
Since (G^{1/2}\otimes I) is invertible,
[
\boxed{\ker M_G=\ker\mathcal C=\bigcap_a\ker C_a.}
]
This is an exact finite-dimensional theorem. Strict positive definiteness matters: for merely semidefinite (G), the kernel can enlarge to include nonzero (mathcal C\psi\) in (ker(G\otimes I)).

I reproduced the registered seed-1729, (8)-dimensional, four-constraint numerical gate with four independently sampled positive (G) matrices. The common kernel rank was 2 in all four cases; maximum constraint-on-kernel residual (1.3519683390710954\times10^{-15}); maximum pairwise zero-projector difference (3.4902265593583572\times10^{-15}); heat-kernel convergence identity maximum error (3.3306690738754696\times10^{-16}). This validates the finite algebra and numerical implementation, not the full gravity constraint family or the refinement/rigging limit.

**Finite C8 relational history.** Reconstructed the exact SymPy matrices for (S_{clock}), (J), (G=S\otimes J), and
[
P_{rel}=\frac18\sum_{t=0}^{7}G^t.
]
All nine core identities passed exactly: (J^2=-I), (J^8=I), (G^8=I), normalized history, (G)-invariance, Hermitian/idempotent combined projector, generation of history from the seed, clock-only projector idempotence, and annihilation of the nontrivial history by clock-only averaging. This remains a positive-control clock, not a derived gravitational clock.

**Finite q=2 metric-source control.** Independently reconstructed (V), relationalized (X,Z,Y), and the frozen shape-to-metric Jacobian. All 18 checks passed exactly: (V^\dagger V=I), (VV^\dagger=P_{rel}), relational observables commute with (G), intertwine and compress to their original operators, and
[
B^TB=\frac92 I_2,\quad
\operatorname{rank}C_{metric}=2,\quad
C C^+ C=C,\quad C^+ C C^+=C^+.
]
The Moore-Penrose inverse is valid only on the two-dimensional shape tangent. This does not add a conformal scalar or a connected inter-block physical kernel.

**Schur algebra stress test.** For a random positive (18\times18) Hermitian block matrix with a six-dimensional P sector, direct minimization over Q reproduced (A-BD^{-1}B^\dagger) for 20 probes with maximum relative discrepancy (7.50\times10^{-16}). A decoupled Q zero mode is compatible with the Moore-Penrose Schur complement; a coupled Q zero mode is not. This is an independent algebra regression, not a BQG amplitude calculation.

### 27.5. Independent higher-shell artifact verification

Downloaded higher-shell artifact ID `10178516490`. From its stored 32x32 complex matrices, I independently reconstructed:
[
K=K^\dagger>0,\qquad
\Lambda=\Lambda^\dagger>0,
\]
with (operatorname{rank}K=32), (lambda(K)\in[4.3060759870,13.3527813527]), and
[
\lambda(\Lambda)\in[10.6357598783,15.0599276660].
]
The stored block-Lanczos matrices satisfy
[
\frac{\|B_1^2-K\|}{\|K\|}=2.94\times10^{-15},
\qquad
\frac{\|B_2^2-\Lambda\|}{\|\Lambda\|}=2.97\times10^{-15}.
]
Thus the reported non-scalar higher-shell (Lambda) is independently reproduced from the archived matrices; its relative distance from a scalar identity is (0.0944046183). Scope remains the finite Euclidean (H_{E0}+H_{E1}) calculation at (J_{max}=5/2).

### 27.6. What the recovered CI does *not* prove — important negative controls

The recovered 30-file core evidence package contains both positive and negative results. The following distinctions are mandatory:

1. **HDA:** the old factorized Hamiltonian is exactly ruled out; the route-normal principal-symbol construction passes, but `full_quantum_HDA_closed=false`. The remaining task is to couple the route factor to the full Peter-Weyl Lorentzian geometry Hamiltonian on the same cylindrical domain and evaluate the off-shell joint commutator.
2. **Regge/EH:** the finite-lattice Regge/EH (c_2) ratios are (0.4539164) at (L=5) and (0.4670130) at (L=6), not the target (0.5). The (c_3/c_2) relative discrepancies are (18.10\%) and (12.83\%). The artifact explicitly says full Einstein IR is not proved; this is not a closure result.
3. **Qubit-to-Einstein:** the oracle geometry reconstruction passes to high precision, but `microscopic_rule_derived=false`; the face-qubit dynamics has not been shown to generate that input geometry.
4. **TT vacuum vs observer smoothing:** the reduced free TT Gaussian vacuum has (P(k)\sim k^{-1}) and block-RMS exponent 1. The observer-smoothing defect has a different exponent near 2.0017. The artifact explicitly rejects identifying the smoothing defect with the free TT vacuum spectrum without a changed interacting/history/RG universality class.
5. **Physical scale:** the translator's (lambda_{R,eff}=1) is an input placeholder, not a value derived by the reduced transfer. Its (eta_2) is a bare reduced-transfer scalar projection; the physical (eta_2), physical six-vector and common scale remain open. The printed SI/LVK values are therefore conditional diagnostics, not frozen BQG predictions.
6. **Joint regulator limit:** the exponent (1.0064429) and zero truncation error above (J_{max}=5/2) (Euclidean) / (13/2) (Lorentzian support bound) apply to the declared finite fixed-input HH family. A uniform theorem for unbounded collective input spin under refinement is still open.
7. **Registered theory status:** the archived integration artifact records 34 structural core gates and `structural_candidate_closed=true`, while `physical_projector_history_closed=false`, `physical_tt_kernel_frozen=false`, `ir_six_wilson_vector_frozen=false`, `common_physical_scale_calibrated=false`, `physicalization_complete=false`, and experimental confirmation false.
8. **L1 metric split:** the (E/T_2) splitting (0.08359565) is a finite Euclidean tangent-Gram precursor, not the dynamic effective-action Hessian or physical quartic anisotropy.

### 27.7. Immediate next calculations

1. Make the five-block producer threshold-stable by re-running it from the frozen upstream source at a cutoff ladder and, if feasible, a no-prune sparse reference. The downloaded NPZ alone cannot reconstruct discarded amplitudes. The local container could not clone GitHub (DNS/network unavailable), so no claim of a fresh source-level rerun is made.
2. Add raw pre-symmetrization Hermiticity measurements to the P producer and rerun its workflow.
3. Build and archive the actual (H_E^{sine}P) and (SP) amplitude images, complete their shared-face recoupling, whiten the target-independent Q support and measure the wall/leakage. This is still the decisive missing microscopic calculation.
4. Preserve the no-go/negative-control results in the integration ledger: full quantum HDA, finite Regge/EH discrepancy, oracle-vs-generated geometry, TT vacuum mismatch and uncalibrated scale.
5. Run the physicalization truth validator on a PR targeting main (or reproduce the full ledger locally) while keeping every physical gate open until the actual graph-changing rigging/history and connected metric response exist.

Classification: successful historical CI for the frozen physical branch; independently recomputed archived finite P and higher-shell matrices; independently reproduced finite master/projector/relational/Schur algebra; several exact finite identities; substantial structural and HDA progress; full physical BQG projector/history, connected physical TT/scalar kernels, joint quantum HDA, IR six-vector, common physical scale and experimental validation remain open.
