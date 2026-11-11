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
