# Branch-wide Python audit — 2026-11-10

## Scope and reproducibility

Repository: `Shtenco/synergy_binary_quantum_theory`  
Audit head on `main`: `c48941f1e0034232402aefa731ef285116ca48a2`  
Audit workflow: [python-all-branches-audit](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42802111888)  
Machine-readable artifact: [PYTHON_ALL_BRANCHES_AUDIT.json](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42802111888/artifacts/13756647023)

The audit enumerated every branch head returned by GitHub, deduplicated identical commit heads, fetched every distinct Python-source blob reachable from those snapshots, and parsed the sources with Python's AST. The scan did not execute repository code.

| Measure | Result |
|---|---:|
| Branches | 157 |
| Unique branch-head commits | 146 |
| Python-file occurrences across branch snapshots | 18,672 |
| Unique Python paths | 672 |
| Unique path/content versions | 753 |
| Unique source blobs fetched | 750 |
| Sources parsed successfully | 748 |
| Syntax errors | 2 |
| Total tree entries across unique branch heads | 46,179 |

The current `main` snapshot is smaller: 723 tree entries and 195 `.py` files. The larger repository-wide count comes from the many historical and experimental branches. The 18,672 figure counts appearances across branch snapshots, not distinct modules; after deduplication there are 672 distinct paths and 753 path/content versions.

## Findings

### 1. Syntax errors in a historical research branch

Two distinct Python blobs fail AST parsing. Both occur in the branch snapshot `research/bqg-v23-exact-finite-certificates-20261108` (commit `5b6d53bfc826193f5008633a753e76cb2458862a`), not in the current `main` snapshot:

- `scripts/pl_lorentzian_group_worker_v2.py:49`: the metadata assignments following `d=json.loads(...)` are indented one level too far, causing `IndentationError: unexpected indent`.
- `tests/test_v21_c1_checkpoint_repair.py:455`: a literal escaped newline was inserted inside a source line, causing `SyntaxError: unexpected character after line continuation character`.

The later V26 branch has a different, parseable worker version; this does not retroactively make the V23 branch executable. Do not treat the V23 branch as a valid runnable checkpoint without repairing or explicitly freezing it as archival.

### 2. Archive extraction and pickle trust boundary

Across branch-head source versions the scanner found:

- 14 archive-extraction call sites using `extractall`;
- 31 `pickle.load/loads` calls across 23 paths.

The most important current-`main` cases are:

- `proof_tools/bqg_mixed_master_aggregate_streaming.py`: tar extraction without a data filter or member-path validation, followed by pickle loading of extracted block maps.
- `proof_tools/bqg_s4sign_aggregate_streaming.py`: the same extraction/deserialization pattern.
- `scripts/download_workflow_artifacts_paginated.py`: ZIP artifact extraction with no explicit validation of the artifact name used to construct the destination directory.

Untrusted pickle data can execute arbitrary code during deserialization. Untrusted tar members can write outside the intended extraction directory, including through link entries if not rejected. These are conditional vulnerabilities when input artifacts can be supplied or influenced by an untrusted party; they are not proof that the current trusted CI artifacts are malicious.

One older `mixed_parallel/bqg_mixed_master_recovery.py` version checks member paths before extraction, but does not fully reject symlink/hardlink entries. That check is not a complete substitute for safe extraction.

**Recommended hardening:** use `tarfile.extractall(..., filter="data")` on supported Python versions; fail closed on unsupported runtimes; explicitly validate ZIP artifact names and member paths; reject link/device members; verify artifact digests/manifests before deserialization; and document that pickle inputs must come only from trusted producers. Prefer a non-executable interchange format for long-lived artifacts where practical.

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

This is a branch-wide static source audit: it establishes source inventory, syntax parsing, and locations of potentially risky APIs. It does not execute all 748 parseable source blobs. Many are heavyweight numerical jobs, artifact aggregators, or scripts that launch processes or contact external services; indiscriminately executing them would be unsafe and scientifically meaningless without their exact inputs and provenance.

The next semantic phase should prioritize (1) the current-main archive/pickle boundary, (2) the V21/V22/V26 Lorentzian worker and collector chain, (3) the depth-6 artifact aggregation/recovery family, and (4) numerical gates with acceptance logic based on assertions or broad exception handling. Every reported numerical result must remain tied to its source commit, inputs, operator version, and artifact hashes.
