# Branch-wide Python audit — 2026-11-10

## Scope and reproducibility

Repository: `Shtenco/synergy_binary_quantum_theory`  
Audit head on `main`: `86dd015d5aefd6181a7301f85c53383e2312f86c`  
Audit workflow: [python-all-branches-audit](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42805938831)  
Machine-readable artifact: [PYTHON_ALL_BRANCHES_AUDIT.json](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42805938831/artifacts/13757833131)

The audit enumerated every branch head returned by GitHub, deduplicated identical commit heads, fetched every distinct Python-source blob reachable from those snapshots, and parsed the sources with Python's AST. The scan did not execute repository code.

| Measure | Result |
|---|---:|
| Branches | 157 |
| Unique branch-head commits | 146 |
| Python-file occurrences across branch snapshots | 18,674 |
| Unique Python paths | 674 |
| Unique path/content versions | 757 |
| Unique source blobs fetched | 754 |
| Sources parsed successfully | 752 |
| Syntax errors | 2 |
| Total tree entries across unique branch heads | 46,182 |

At the original audit head, `main` contained 728 tree entries and 198 `.py` files. A later main snapshot at commit `99dd3c2a31610b0ab9b9c56cd1b3937cddda070c` contained 739 tree entries and 207 `.py` files. The latest checked main tree at commit `b965eb8992c517ed944d5c33e5d45094beb451e1` contains 744 tree entries, 210 `.py` files, and eight tracked ZIP/archive files; the expanded branch-wide refresh for the newer head is tracked separately. The larger repository-wide count comes from the many historical and experimental branches. The 18,674 figure counts appearances across branch snapshots, not distinct modules; after deduplication there are 674 distinct paths and 757 path/content versions. There are 754 distinct blob contents, of which 752 parsed successfully and two failed syntax parsing.

## Findings

### 1. Syntax errors in a historical research branch

Two distinct Python blobs fail AST parsing. Both are absent from the current `main` snapshot (the refreshed audit parsed all 198 main Python files with zero syntax errors). The first malformed blob is shared by four branch heads — `research/bqg-v23-exact-finite-certificates-20261108` (`5b6d53bfc826193f5008633a753e76cb2458862a`), `research/bqg-v24-finite-index-constraint-gates-20261108` (`524240db085485ac31678e82778f585ffd38d199`), `research/bqg-v25-local-HDA-charge-20261108` (`b58872dcb23404fd509421e14b975a98d2214a3e`), and `audit/bqg-independent-proof-20261108` (`fddd49427c636e92dfe79b96034862aeba629927`). The second malformed blob occurs on `fix/v21-c1-fail-closed-20261108` (`2626dcbed6bf6d91179810792bb7db5343b08477`):

- `scripts/pl_lorentzian_group_worker_v2.py:49`: the metadata assignments following `d=json.loads(...)` are indented one level too far, causing `IndentationError: unexpected indent`.
- `tests/test_v21_c1_checkpoint_repair.py:455`: a literal escaped newline was inserted inside a source line, causing `SyntaxError: unexpected character after line continuation character`.

The later V26 branch has a different, parseable worker version; this does not retroactively make the V23 branch executable. Do not treat the V23 branch as a valid runnable checkpoint without repairing or explicitly freezing it as archival.

### 2. Artifact ingestion and pickle trust boundary — refreshed status

The initial branch-wide audit correctly identified unsafe archive extraction,
unrestricted pickle deserialization, unvalidated artifact directory names, and
missing digest checks in the then-current aggregator path. **That description
is superseded by the hardening now present on main.**

Current protections:

- `proof_tools/safe_archive.py` rejects absolute and parent-traversal paths,
  backslashes, duplicate destinations, links and device members. Tar extraction
  requires Python's `filter="data"`; ZIP extraction is performed member by
  member and rejects links/special entries.
- `proof_tools/safe_pickle.py` uses a restricted unpickler. It permits a
  bounded set of inert built-ins and, when explicitly requested, NumPy numeric
  arrays/scalars. It rejects unlisted globals and validates that arrays are
  numeric and non-object. Both current-main depth-6 streaming aggregators now
  call `safe_load_numpy_pickle`, not unrestricted `pickle.load`.
- `scripts/download_workflow_artifacts_paginated.py` validates artifact names,
  streams ZIP downloads to disk, checks the GitHub-reported size and SHA-256
  digest, and uses safe ZIP extraction. Its redirect handler refuses non-HTTPS
  redirects and strips the GitHub bearer token on cross-origin redirects.
- `scripts/extract_base64_tar_safe.py` applies compressed and uncompressed
  archive-size limits before safe extraction of the tracked proof bundles.

These are substantive mitigations, but not a proof that pickle is an ideal
long-term interchange format. Restricted unpickling still has a parser and
resource-exhaustion attack surface: a hostile numeric payload could attempt
excessive allocation or CPU use before post-load object validation. The
artifact ZIP downloader verifies the compressed artifact digest and size, but
the generic ZIP/TAR extractors do not impose a universal total uncompressed
byte/member-count cap. The base64 proof-bundle extractor has explicit caps;
the workflow-artifact extraction path should still gain producer-specific
uncompressed limits.

Remaining provenance and integrity work:

1. Replace pickle maps with a versioned, non-executable numeric format
   (NPZ with `allow_pickle=False` or a validated tensor format) when practical.
2. Add a content-addressed per-shard manifest binding producer repository,
   source commit, workflow run/attempt, shard index, artifact ZIP digest,
   extracted map/summary digests, and operator version.
3. The S4-sign aggregator pins the exact triggering producer SHA. The mixed
   Stage A/Stage B workflows now pin the frozen source commit
   `f6bfcd775b8c3cd33c6c235d48293e40657baf9d` and verify both
   `git rev-parse HEAD` and the mixed-bundle Git blob SHA
   `cbbaa80f6b6b61d353458db1fa1f03f824730532`. The earlier mutable-branch
   checkout concern for these two stages is therefore closed. The remaining
   gap is the end-to-end per-shard manifest tying each output map and summary
   to that pinned producer, run/attempt, shard index, and operator version.
4. The mixed-master aggregator now calls
   `validate_shard_coverage` on summary shard IDs and
   `validate_index_coverage` on both record and map indices. The shared guard
   requires exact zero-based coverage with no duplicates, missing indices,
   booleans, or noninteger values. Regression tests cover duplicate/missing
   shards and record/map gaps. This previously open coverage gap is closed.
5. Preserve read-only permissions, disabled checkout credentials, and
   producer-branch restrictions in artifact-processing jobs.

The risk should therefore be described as a **remaining conditional
supply-chain and resource-exhaustion concern**, not as unrestricted
pickle-code-execution still being present on main. No evidence has been found
that any existing artifact is malicious.

### 3. Exception handling

The branch-wide scan found 85 `except Exception` handlers across 64 paths and one bare `except:`.

Many broad handlers are intentional worker boundaries: they serialize an exception into a failed shard result or wrap it with a more specific failure. Those are not automatically defects. The review should distinguish them from handlers that silently continue or convert computation failures into apparently valid defaults.

The bare handler is in `aggregate_existing/bqg_s4sign_aggregate_capacity.py:61`; it catches every `BaseException` while parsing a block filename and continues. It should be narrowed to `ValueError` (and any other explicitly expected parsing exception).

### 4. Process execution and network access

The scanner found 48 subprocess call sites across 15 paths and 17 network-call sites across 14 paths.

Manual inspection of the 15 subprocess-bearing source versions found no `shell=True` and no `os.system`. Most subprocess calls invoke Git with argument arrays to attest the source commit/worktree; one launches a local llama.cpp server. This is materially safer than shell-string execution, although the launched binary and arguments still need provenance controls.

Network calls are concentrated in GitHub API artifact/metadata retrieval and local model-server communication. Network access itself is not a defect; the relevant controls are destination allowlisting, TLS verification, response-size/time limits, and digest validation. The current audit tool's own GitHub API call is included in this count.

### 5. Assertions and numerical gates

There are 610 `assert` statements across the distinct source versions. Assertions are useful in tests, but Python removes them under `python -O`. Acceptance criteria that must remain active in every execution mode should use explicit checks that raise exceptions, not bare assertions. This is a review priority for scientific certification and artifact-validation paths, not a claim that all 610 assertions are unsafe.

### 6. Dependency reproducibility

`requirements.txt` uses lower bounds rather than a lock or exact constraints. The current CI resolver selected newer major/minor versions (including NumPy 2.x and pandas 3.x) and `torch` pulled a large CUDA dependency stack. The repository does use `torch`, `pandas`, `networkx`, and `sklearn` in some branch versions, so simply deleting these dependencies would be incorrect.

**Recommended hardening:** maintain a pinned lock/constraints file for reproducible proof runs; separate the minimal numerical-gate dependencies from optional ML/courier dependencies; and pin CPU-only Torch where GPU execution is not required.

## CI observations

- The branch-wide audit workflow completed successfully and produced the artifact linked above.
- `numerical-continuation` completed successfully on the audit head.
- A recent `core-regression` run failed at `scripts/audit_core_scope.py` because five Markdown files contain retired vocabulary; the subsequent `Compile all Python sources` step was skipped. This is a documentation-scope gate failure, not a Python syntax failure, but it means that CI run did not execute the compile step.
- Higher-shell and L1-q4 jobs have also failed while downloading artifacts from the legacy `Shtenco/binary_quantum_theory` repository (HTTP 301). These are artifact-source/redirect failures, not evidence of a Python numerical regression.

## Limits of this audit

This is a branch-wide static source audit: it establishes source inventory, syntax parsing, and locations of potentially risky APIs. It does not execute all 752 parseable source blobs. Many are heavyweight numerical jobs, artifact aggregators, or scripts that launch processes or contact external services; indiscriminately executing them would be unsafe and scientifically meaningless without their exact inputs and provenance.

The next semantic phase should prioritize (1) the current-main archive/pickle boundary, (2) the V21/V22/V26 Lorentzian worker and collector chain, (3) the depth-6 artifact aggregation/recovery family, and (4) numerical gates with acceptance logic based on assertions or broad exception handling. Every reported numerical result must remain tied to its source commit, inputs, operator version, and artifact hashes.


## Phase 2 — semantic review of the current-main rank-certification path

A focused review of both current-main streaming rank aggregators found a missing matrix-shape invariant. They computed an SVD and accepted a block when rank equals m and the smallest returned singular value exceeded tolerance, but did not first require the stacked image matrix to have exactly m columns. A matrix with more than m columns can have rank m while still having a nontrivial right kernel; in that case rank equals m is not a full-column-rank certificate for the claimed m-column orbit block.

This was a fail-open validation gap for malformed or mismatched artifact maps. It does not establish that any existing artifact had the wrong shape, but the code did not enforce the invariant needed by its own certificate.

**Correction committed on main:**

- Added proof_tools/rank_matrix_guard.py, which requires a 2D matrix, exact expected column count, and finite entries.
- Applied the guard before SVD in both proof_tools/bqg_mixed_master_aggregate_streaming.py and proof_tools/bqg_s4sign_aggregate_streaming.py.
- Added tests/test_rank_matrix_guard.py, including a regression where a 2×3 matrix has rank 2 but must not pass as a two-column full-rank certificate.
- The full numerical-continuation workflow passed on commit e371f1c6d2ebece81abe760e50507d5c31f7949f, including python -m unittest discover -s tests and the existing finite-evidence checks.

No operator formula, coefficient, basis, threshold, or scientific acceptance tolerance was changed. The new check only rejects malformed matrix dimensions and nonfinite values before the existing SVD criterion is applied.

### Remaining rank-pipeline provenance issues

The same two aggregators still consume tarballs and pickle maps without a content-addressed manifest binding every shard, summary, and pickle map to its producer commit. Summary counts and orbit coverage are checked, but those checks do not prove that the loaded numerical maps are the exact outputs of the intended producer. Since pickle deserialization is executable, these inputs must remain trusted until a digest-bound, non-executable interchange format is introduced.

The SVD certificate remains a **finite numerical rank certificate**, not an exact symbolic rank theorem. The implementation uses floating-point singular values and a numerical tolerance. Its own proof_logic string correctly calls the master image numeric; downstream documentation must preserve that limitation.

The output metadata field unique_rows / unique_q is populated from the geometric unique-support set before filtering to map keys. The actual matrix row count is separately recorded as rows. These fields should be renamed or split in a later schema revision so that geometric support counts cannot be mistaken for the number of loaded nonzero row blocks.


## Phase 3 — scope reconciliation and artifact trust-boundary review

### Correct interpretation of the 1,500+ count

The authoritative artifact from run 42802977103 reports 18,674 Python-file appearances across 157 branch snapshots. These are **branch-snapshot appearances**, not 18,674 distinct modules. After deduplication there are 674 distinct Python paths and 757 path/content versions (754 distinct blob contents). The current main head at this audit refresh has 728 tree entries and 198 `.py` files.

The source inventory includes `.py`, `.pyi`, `.pyx`, `.pxd`, and `.ipynb` blobs. It does not yet count Python source embedded in Markdown/YAML code blocks or source bundles stored as archives/base64. In particular, the `research/depth6-s4sign-proof` branch contains base64-encoded `.tar.gz` bundles that its workflows unpack before executing Python workers. Their contents are outside the 754-blob source count. Therefore this audit resolves the large count as branch-snapshot scale, but it does **not** certify the total number of embedded Python modules until those bundles and inline code blocks are inventoried separately.

### Artifact ingestion is a high-priority security boundary

The current-main `bqg-depth6-s4sign-aggregate.yml` workflow is triggered by `workflow_run`, checks out the orchestration branch and the `research/depth6-s4sign-proof` source branch, downloads artifacts from the triggering run, and then invokes an aggregator that extracts tarballs and deserializes pickle maps. The mixed-master aggregation workflow follows the same general pattern.

The current-main aggregator paths still contain:
- `tarfile.extractall` without explicit member-path/link rejection or the `data` extraction filter;
- `pickle.load` on block maps extracted from workflow artifacts;
- no per-shard content-addressed manifest that binds the loaded map and summary to the producer commit and exact artifact digest.

The downloader also uses the GitHub artifact name directly in `args.out / name` and then calls `ZipFile.extractall`. Artifact names should be treated as untrusted path components and validated before creating directories. The downloader currently does not compare the downloaded ZIP bytes against the artifact digest returned by GitHub.

This is a **conditional supply-chain/code-execution risk**, not evidence that any existing artifact is malicious. A party able to influence the producer branch or artifact contents could potentially cause code execution through pickle deserialization in the downstream aggregation job. Read-only workflow permissions reduce impact but do not make deserialization safe; checkout credentials should also not persist in workspaces that process untrusted artifacts.

Recommended remediation order:
1. Stop loading pickle from artifacts; migrate numeric matrices to a non-executable format such as NPZ with `allow_pickle=False`, or a validated binary tensor format.
2. Until migration, restrict producer branch/run provenance, verify artifact digests, bind each shard to its producer commit and workflow run, and keep the artifact-processing job free of checkout credentials and unrelated secrets.
3. Replace tar extraction with explicit member validation, reject absolute/traversal paths, symlinks, hardlinks and device entries, and use `filter="data"` on supported Python versions.
4. Validate artifact names as single safe path components and validate ZIP members before extraction; stream downloads to disk with size limits rather than buffering an entire archive in memory.
5. Make the mixed aggregator verify the same exact shard-index coverage that the S4-sign aggregator already checks.

### Additional semantic observations

- `proof_tools/bqg_s4sign_aggregate_streaming.py` verifies that shard IDs equal `range(expected)`; the mixed aggregator checks the number of summaries and global block/map coverage but does not independently enforce unique, complete shard IDs. Add that invariant for symmetry and provenance clarity.
- Both streaming aggregators report `unique_q` / `unique_rows` from the geometric unique-support candidate set, while `rows` is the actual stacked matrix row count. Preserve the existing result schema for old artifacts, but rename or split these fields in a versioned schema so that support candidates are not confused with loaded matrix rows.
- The current-main numerical continuation tests pass on the recorded audit lineage. The branch-wide audit remains static: it does not prove numerical correctness of all source versions, and no scientific acceptance threshold should be inferred from AST parsing.


### Audit-tool reproducibility correction

The first branch-wide scanner revision deduplicated by blob SHA but attached each finding to the first path encountered for that blob. Since branch trees are fetched concurrently, that representative path could depend on completion order and could omit other branch/path locations for the same source blob.

The scanner was corrected on main at `86dd015d5aefd6181a7301f85c53383e2312f86c`:
- every finding now carries all known `(commit, path, extension)` locations for its blob;
- the representative path is selected deterministically;
- the report includes the exact branch-head map and a separate current-main AST result;
- the audit exits nonzero if current main has syntax errors, while historical branch syntax errors remain reported without blocking main.

The refreshed workflow passed. It confirmed `main_ast_status=PASS`, `main_syntax_errors=0`, 198 main Python files, and the same branch-wide inventory: 157 branches, 146 unique heads, 18,674 file appearances, 674 paths, 757 path/content versions, 754 distinct blobs, 752 parseable sources, and two historical syntax-error blobs.


## Phase 4 — pin the exact depth-6 S4-sign producer commit

The `workflow_run` aggregation workflow previously checked out the mutable `research/depth6-s4sign-proof` branch tip while downloading artifacts from a particular completed producer run. Those are not necessarily the same revision if the branch advances before aggregation. This is a source-to-evidence provenance gap independent of whether the numerical matrices are correct.

**Correction on main:** `.github/workflows/bqg-depth6-s4sign-aggregate.yml` now (a) requires a successful triggering run from the same repository and the expected proof branch, (b) checks out the triggering run's exact `head_sha`, and (c) verifies that `git -C proofsrc rev-parse HEAD` equals that SHA before unpacking its proof bundle. This changes no scientific operator, numerical tolerance, or certificate claim.

The aggregator still needs a per-shard manifest binding producer SHA, shard index, raw artifact digest, extracted matrix-map digests, and operator-version metadata. GitHub artifact ZIP digest verification and producer checkout pinning are necessary but not sufficient for end-to-end provenance. The source checkout change has been committed, but the next producer-triggered aggregation has not yet been observed to complete; do not label it a verified end-to-end PASS.


## Phase 5 — finite RG completion and CI/provenance follow-up

### Five-step multiplicity-channel ancestry

The final full channel-flow shard, (s_2=7), completed in workflow run
[42797583866](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42797583866).
Its matrix is persisted at
`docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_7.json`. The full
sequence (s_2=3,4,5,6,7) now has five consecutive passed matrices. The
maximum-overlap dynamic-programming path is

[
a_0^{(3/2)}	o a_1^{(2)}	o a_1^{(5/2)}
	o a_1^{(3)}	o a_2^{(7/2)}	o a_2^{(4)}.
]

Its cumulative overlap product is 0.49282787562801145; the all-low-index
shortcut gives 4.694570240506031e-5, a ratio of 10497.8273. The path, five
matrices, and regression test are persisted in
`BQG_FULL_CHANNEL_FLOW_MATRIX.md`,
`BQG_RG_ANCESTRY_PATH.md`,
`docs/numerical_continuation/BQG_RG_ANCESTRY_PATH_5STEP.json`, and
`tests/test_bqg_rg_ancestry_path.py`.

This remains a finite representation-RG projector-overlap result. It is not
a physical transition probability, a physical-time history, a continuum
fixed point, or a physical-projector convergence theorem.

### Mixed depth-6 workflow checkout correction

The aggregate jobs in
`.github/workflows/bqg-depth6-mixed-stage-a.yml` and
`.github/workflows/bqg-depth6-mixed-stage-b-final.yml` had duplicate
`path` keys in a single checkout `with` mapping. The correction leaves the
main orchestration checkout at the workspace root and places the research
source checkout explicitly at `proofsrc`. The existing bundle and aggregator
paths then resolve unambiguously. This is a workflow-layout correction only;
it changes no operator, matrix, rank tolerance, or scientific criterion.

The aggregate checkout blocks in both workflows were corrected to keep the
main orchestration checkout at the workspace root and the frozen proof source
at `proofsrc`; duplicate `path` keys were removed. Stage A run
[42819682938](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42819682938)
was queued when the initial path fix was made. The first Stage B push-triggered
workflow file-validation run occurred before its corresponding path correction
and failed without jobs; it is not evidence of a numerical failure. The later
workflow definitions now pin and verify the frozen source commit and bundle
hash in each shard and aggregation job.

### W5 tolerance ablation: timeout and one-pass replacement

The original eight-pass W5 tolerance-ablation run
[42795399991](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42795399991)
was cancelled at the 35-minute job limit before its result JSON was created.
It produced no ablation evidence.

To address the runtime bottleneck without relaxing the scalar guard, a new
one-pass runner was added:
`scripts/w5_scalar_tolerance_ablation_onepass.py`. It generates each
T-sequence branch once, evaluates both projection tolerances on that same
branch, and feeds production/zero accumulator variants in parallel. It
preserves the captured key, operator path, sequence order, coefficients,
internal basis cutoff and scalar-preservation guard. The workflow now invokes
this runner; its run
[42820210450](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42820210450)
was queued when this note was written. Until that run completes and its
per-variant full-map comparisons are inspected, the pruning/aggregation cause
of the W5 residual remains unresolved.

### Current finite carrier-leakage computation

The compressed-return (L,D) calculation for the 243-dimensional coarse
carrier and 32-dimensional embedded fine carrier is still running in
[42797344820](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42797344820).
No value of the off-carrier leakage, Galerkin defect, or Hermiticity
consistency check is claimed before its artifact is available. Even when
completed, its scope is the compressed local-singlet return operator; it is
not the full graph-changing habitat residual (E_a,F_a).

### Critical numerical gates under Python optimization

A separate regression now checks that the canonical refinement-residual gate
uses explicit runtime checks rather than optimization-sensitive `assert`
statements, and that it still passes under `python -O`. This is a focused
guard improvement; it does not change the algebraic theorem or numerical
acceptance thresholds.


### Direct multi-step composition gate

The five-step local ancestry product raises a separate mathematical issue:
the product of scalar one-step channel matrices is not automatically the
overlap of the composed microscopic refinement map. The matrix product
(F_jF_{j+1}) corresponds to an intermediate channel measure-and-prepare
map that discards off-diagonal channel coherence and traceless structure
inside each multiplicity block.

This is now specified in
`BQG_MULTISCALE_CHANNEL_COMPOSITION_GATE.md`. The executable diagnostic
`scripts/bqg_two_step_channel_composition.py` composes the actual generated
rank-2 carrier subspace through two adjacent microscopic refinement tensors,
compares the direct endpoint matrix (G_{j:j+2}) with (F_jF_{j+1}), and
records intermediate off-diagonal and internal-anisotropy block norms. It
recomputes both one-step matrices and checks them against the persisted JSON
references before interpreting the two-step defect. Workflow
[42821288997](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42821288997)
was queued at the time of this update; no factorization result is inferred
before its artifacts arrive. The maximum-product five-step path therefore
remains a local ancestry candidate until this stronger composition gate is
measured.


## Phase 6 — expanded embedded-source audit and direct two-step composition

### Branch-wide scanner scope expansion

The branch-wide audit runner now scans more than standalone Python files. It
also:
- inventories Python source members in tracked TAR/ZIP/WHL/EGG archives,
  including base64-encoded archives, with explicit archive/member/total-source
  size and member-count limits;
- extracts explicitly marked Python fences from Markdown/RST and Python
  heredocs from YAML/shell sources, dedents them, and parses them as independent
  AST units;
- reports main-only finding counts and inline-source syntax errors separately
  from branch-wide findings.

The current main snapshot contains eight tracked ZIP/archive files, including
the compact-holonomy source archives and saved-vector bundle. These are now
inside the scanner's supported archive scope. The scan remains static: it
never imports or executes the extracted source. Its scope explicitly excludes
nested archives, unmarked prose fragments, unsupported archive formats, and
non-head historical commits.

The expanded audit runs were repeatedly cancelled while main continued to
advance. The latest observed run,
[42822233592](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42822233592),
was cancelled before producing an artifact. No new archive-member or
inline-Python counts are claimed until one run completes on a stable head.

### Two-step microscopic channel composition

The five-step ancestry score is the maximum product of local one-step
overlaps. The matrix product (F_jF_{j+1}) is a different object: it sums over
intermediate channels and corresponds to a channel measure-and-prepare map
that discards off-diagonal channel coherence and traceless structure within
each multiplicity block. The maximum-product path score instead selects one
intermediate channel. Neither is automatically the direct overlap of the
composed microscopic refinement.

The new gate
`BQG_MULTISCALE_CHANNEL_COMPOSITION_GATE.md` defines the direct two-step
endpoint matrix (G_{j:j+2}), the scalar matrix-product defect
(G_{j:j+2}-F_jF_{j+1}), and the intermediate off-diagonal/internal
anisotropy diagnostics. The runner
`scripts/bqg_two_step_channel_composition.py` composes the actual generated
rank-2 carrier subspace through two adjacent microscopic blocking tensors,
recomputes and verifies the one-step matrices against their persisted JSON
references, and separately records the local max-product path score matrix.
Its workflow run
[42821558786](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42821558786)
was pending at this update. No factorization outcome is inferred before the
two matrix artifacts are inspected.

This is the next decisive refinement test: it determines whether the
five-step local ancestry score is supported by coherent multi-scale
composition, or whether the scalar channel matrices only describe a
channel-twirled surrogate.


## Phase 7 — inline-source audit implementation and workflow verification

### Current scanner revision

The current branch-wide scanner revision is on main commit
`033ca9779cdd61299d048d9dfcbe884f0d7f9692`. In addition to standalone
Python files, notebooks and direct archive members, it extracts labelled
Python fences from Markdown/RST and Python heredocs from YAML/shell files.
It now parses the same Git blob once per distinct source suffix, so identical
content stored as both Markdown and YAML is not interpreted using only the
first path encountered. Findings are mapped back to the containing file's
original line numbers and restricted to locations with the matching suffix.

Regression tests in `tests/test_python_audit_inline_blocks.py` cover:
- Markdown Python fences with language metadata;
- YAML `python - <<'PY'` heredocs;
- indentation normalization;
- separate parsing of multiple code fences.

The audit reports both repository-wide and main-only finding counts, and
separates inline syntax errors from standalone-module syntax errors. The
current audit run is
[42822081116](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42822081116);
it was pending when this note was written. Earlier expanded-scan runs
42821334886, 42821558904 and 42821910697 were cancelled by successive pushes before
producing a completed artifact. Do not treat their missing artifacts as
audit failures or as successful scans.

### Optimization-sensitive gate correction

The canonical refinement residual gate previously used nine Python
`assert` statements for the isometry, projector, residual decomposition,
orthogonality, norm bounds and registered scalar check. Python removes
these statements under `-O`, but the script still printed `PASS`.

The gate now uses explicit `require(condition, message)` checks. The new
test `tests/test_critical_gate_optimization.py` verifies that the gate's AST
contains no `ast.Assert` node, that a false condition still raises under
optimization, and that the full gate completes under `python -O`.
`numerical-continuation.yml` invokes the gate. No operator, identity,
coefficient, basis, or tolerance changed.

A second active numerical-continuation gate,
`docs/numerical_continuation/verify_volume_blocks.py`, also had six
optimization-sensitive assertions covering the recoupling-basis isometry,
volume-operator carrier preservation, characteristic-polynomial identities,
and agreement with frozen exact blocks. These are now explicit `require`
checks as well. The optimization regression test parses this script's AST
and executes it under `python -O`. The exact SU(2) operators, matrices,
and symbolic equalities are unchanged.

Five additional standalone scientific gates still retain assertion-based
checks and print `PASS`. Rather than mechanically rewriting their symbolic
and numerical expressions, each now refuses to run when `__debug__` is
false: `bqg_three_node_reduced_chain_gate.py`,
`bqg_local_information_resistance_gate.py`,
`bqg_two_node_reduced_gluing_gate.py`,
`bqg_projector_graph_transition_gate.py`, and
`bqg_minimal_entanglement_volume_gate.py`. Their formulas and assertion
tolerances are unchanged; this guard prevents a false success under
`python -O`. The optimization test suite checks for the guard in each
module. The remaining assert-bearing main files are not yet all migrated.

This is one corrected gate, not a blanket clearance of the other
assert-based gates. The branch-wide report previously counted 171
production assertion nodes on its main snapshot. Each remaining
main-only gate must be reviewed for whether its assertions are
certification criteria, internal debugging checks, or ordinary tests.

### Depth-6 workflow provenance and trigger behavior

The mixed Stage-A, Stage-B and keymeta-salvage workflows now pin the
research source to commit
`f6bfcd775b8c3cd33c6c235d48293e40657baf9d` and verify the mixed bundle
blob SHA
`cbbaa80f6b6b61d353458db1fa1f03f824730532`. The keymeta-salvage
workflow additionally verifies the exact harvester, extractor and safe
archive/pickle helper blobs. The regression test
`tests/test_depth6_workflow_source_pinning.py` locks these values.

A push to the keymeta-salvage workflow file had previously invoked a
fallback to stale run ID `36899125190`. That run failed only after the
pinned source and helper checks succeeded; the harvester correctly
rejected the stale run as not being a successful same-repository main
Stage-A run. The fallback has been removed, the push trigger has been
removed, and salvage now requires a genuine Stage-A `workflow_run` or an
explicit `workflow_dispatch`.

The mixed Stage-A workflow's push-to-own-YAML trigger was also removed.
Previously, merely editing that workflow launched its 30-job matrix.
The upstream Stage-A `workflow_run` trigger remains. The queued push
matrix from before that trigger correction is not evidence of a
completed proof run; its source-pinning checks have not yet been
confirmed by a completed matrix job.

### CI state at this update

- Branch-wide Python audit: run 42822313742 pending.
- Numerical continuation: latest run 42821910672 queued.
- Core regression: latest run 42821910573 queued.
- Mixed Stage-A push matrix: run 42821076917 queued; it predates removal
  of the self-trigger and has not completed its source-pinning step.
- Keymeta salvage: run 42821051288 failed at the stale run-ID validation;
  source commit and bundle verification steps succeeded.

The CI queue has not yet supplied completed results for the latest
scanner and optimization-mode tests. The code changes are committed,
but their CI status must remain pending until the workflows complete.


### Core-regression compile ordering

The core workflow previously ran the documentation scope audit before
`python -m compileall -q .`. A documentation-only failure therefore
skipped the Python compilation step. The order is now reversed: the full
main-tree compile runs first, followed by the documentation scope gate.
This keeps syntax/bytecode compilation evidence independent of
documentation-policy failures. It changes no scientific test or acceptance
criterion.


## Phase 8 — strict artifact-index validation in both depth-6 aggregators

A second semantic review of the two active streaming aggregators found a
coercion gap in artifact metadata. Both paths used `int(value)` on JSON
record indices and `b<index>.pkl` filenames. That accepts values such as
`true`, numeric strings, signed indices and noncanonical leading-zero
filenames as if they were ordinary integer indices. Exact set coverage
after coercion is not enough to prove that the producer supplied the
required schema.

The rank guard now provides:
- `validate_index_value`: accepts only a nonnegative JSON integer and
  explicitly rejects booleans and numeric strings;
- `validate_positive_integer`: requires strictly positive dimensions;
- `parse_indexed_map_filename`: accepts only canonical ASCII
  `b0.pkl`, `b1.pkl`, … index stems and rejects leading zeros, signs,
  non-ASCII digits and malformed names.

Both `bqg_mixed_master_aggregate_streaming.py` and
`bqg_s4sign_aggregate_streaming.py` now validate raw record indices and
positive dimensions before storing records, parse map filenames through
the strict helper, and enforce exact zero-based shard/record/map coverage.
The S4-sign path now uses the same shared coverage guard as the mixed
aggregator rather than its earlier local integer-coercion check.

During integration, a manual review caught that the script-mode fallback
imports initially omitted the new helper names. Those imports are now
corrected in both aggregators; this matters because the workflow executes
these files as scripts, where relative package imports can fall back to
the top-level import branch. Regression tests in
`tests/test_rank_matrix_guard.py` cover booleans, strings, negative and
leading-zero indices, non-ASCII digits, and nonpositive dimensions.

These checks reject malformed artifact metadata only. They do not change
the master operator, support graph, SVD threshold, rank criterion, or
scientific result. The first CI run containing the complete strict-index
change is
[42823223090](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42823223090);
it was pending when this phase was recorded.


## Phase 9 — ZIP Unix file-type validation

The archive helper's ZIP path previously rejected symlinks, character/block
devices and FIFOs, but did not reject every other special Unix file type
(e.g. a socket), and it did not verify consistency between Unix mode bits
and the ZIP filename's directory marker. The extractor would write such an
entry as an ordinary file rather than create a socket, so this was not a
path-escape or code-execution primitive; it was nevertheless inconsistent
with the helper's documented regular-file/directory-only contract.

`safe_extract_zip` now uses an allowlist of ZIP Unix file types: unknown or
special types are rejected; only absent type metadata, regular files, and
directories are accepted. It also rejects a directory mode without a
directory name and a regular-file mode with a directory name. Regression
tests cover socket metadata and type/name mismatch. Normal Windows ZIP
entries with no Unix type bits remain accepted.

The fix changes no legitimate file payload or extraction destination. The
branch-wide audit run containing this change is
[42823223090](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42823223090);
it was pending when this phase was recorded. Numerical continuation and
core regression for the same commit were queued as runs
[42823223134](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42823223134)
and
[42823223142](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42823223142).
