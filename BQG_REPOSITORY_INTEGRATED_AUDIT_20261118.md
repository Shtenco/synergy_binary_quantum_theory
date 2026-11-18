# Integrated BQG repository audit — 2026-11-18

## Executive result

The research branch is substantially healthier as a software artifact than the previous Q2-only CI view suggested: its Python sources compile, the canonical structural/physical ledgers validate, the documentation checks pass, the complete `tests/` suite has a green run, the Q2 export workflow produces valid artifacts, and the all-branch Python AST audit completes.

That does **not** mean the physical theory is closed. The machine-readable ledgers still report an open structural candidate and an open physicalization frontier. The W5 production calculations also remain incomplete: several expensive jobs were interrupted before emitting the artifacts needed for a production/shadow comparison.

Audited branch: `work/q2-two-sided-singlet-qfi-20261118`  
Main branch: not modified; no merge performed.

## 1. Reproducible software evidence

### 1.1 Full branch regression

Workflow: [bqg-branch-integrity-audit, run 44382012617](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44382012617)

- `python -m compileall -q .`: PASS.
- `audit_core_scope.py`: PASS.
- `verify_theory_gates.py`: PASS.
- `verify_depth6_frontier.py`: PASS.
- `verify_physicalization_gates.py`: PASS.
- `verify_readme_links.py` and the selected GitHub-LaTeX checks: PASS.
- Full `pytest -q tests`: **509 passed, 59 subtests passed, 0 failures, 0 errors, 0 skips** (354.75 seconds).

The successful test run is on commit `a39fe4ef5335a6bab31c347599b777a5895ff954`. Subsequent commits up to the audit snapshot changed documentation/workflow configuration, not the tested Python implementation. A second full branch-integrity run was started for the documentation-updated snapshot.

### 1.2 Q2 graph/QFI export

Workflow: [q2-pl-graph-export, run 44382474317](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44382474317)  
Artifact: [q2-pl-graph-export](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44382474317/artifacts/14464941457)

The artifact contains 23 JSON exports. All 23 were parsed independently; none has `passed=false`.

The exhaustive reduced-boundary sweep reports:

- 144 boundary subsystems;
- 2,880 partition-order evaluations;
- maximum reduced-density error: `3.8860812576573896e-16`;
- maximum reduced-derivative error: `1.6667795129466247e-16`;
- maximum mixed-QFI error: `9.43689570931383e-16`;
- minimum eigenvalue of pure-minus-mixed QFI: `0.0`.

The Haar pushforward control covers seven finite configurations. The product-Haar support is a Segre variety inside the fusion image; the measured complex projective codimensions are 1, 4, 11, 4, 20, 9, and 16 for the tested `(j, vertex_count)` configurations. This establishes non-closure of the specified independent local Haar toy ensemble under the tested fusion map. It does **not** identify the physical BQG measure.

### 1.3 All-branch Python static audit

Workflow: [python-all-branches-audit, run 44383335720](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44383335720)  
Artifact: [PYTHON_ALL_BRANCHES_AUDIT.json](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44383335720/artifacts/14464687471)

- 174 branch heads; 163 unique head commits;
- 22,809 Python-file occurrences across branch snapshots;
- 814 unique Python paths and 1,006 path/content versions;
- 1,003 unique Python blob contents; 1,564 sources parsed, including embedded and inline sources;
- 10 embedded archives scanned, 983 members and 335 embedded Python sources;
- 0 unsafe archive members and 0 archive extraction errors;
- `main_ast_status=PASS`, with 0 syntax errors;
- current work-branch HEAD at scan time: 0 syntax errors.

Two syntax errors remain only in historical branch snapshots:

1. `.github/workflows/lorentzian-pruned-first-column.yml#heredoc-PY-2:86`, on `recovery/bqg-physical-spectral-history-2026-09-06` and `research/physical-scalar-kernel-dynamics`.
2. `scripts/pl_lorentzian_group_worker_v2.py:49`, on `research/bqg-v24-finite-index-constraint-gates-20261108`, `research/bqg-v25-local-HDA-charge-20261108`, and `audit/bqg-independent-proof-20261108`.

The current branch has heuristic AST findings (17 broad-exception handlers, 2 network-access sites, and 9 process-execution sites). The network sites are in archived documentation code; process calls occur in audit/test/archived contexts and local worker orchestration. These are review markers, not confirmed vulnerabilities. The current branch has no detected pickle deserialization, bare `except`, or unguarded archive-extraction finding.

## 2. Current mathematical status from the canonical ledgers

The current `theory_gates.json` contains 50 core gates:

| Status | Count |
|---|---:|
| `proved` | 23 |
| `tested_finite` | 21 |
| `conditional` | 5 |
| `no_go` | 1 |

`core_theory_closed_declared=false`; therefore `structural_candidate_closed=false`. The core no-go is `PURE_ANCILLA_NODE_REFINEMENT`: a fixed pure four-qubit Gauss-singlet ancilla does not induce the canonical local (S_4)-equivariant refinement map. The separate bilinear refinement construction has a passing finite algebraic gate and now has explicit regression coverage. That does not prove dynamical preparation of the required physical state.

The physical ledger still has all 9 mandatory physicalization gates and all 4 supplemental physical gates open. `physicalization_complete=false` and `experimentally_confirmed=false`.

## 3. Q2 fusion: what is established and what is not

The finite vertex-fusion implementation contracts internal singlet edges and preserves the resulting boundary tensor (up to the declared global-phase convention) and the QFI of retained inter-block edge parameters in the tested finite fixtures. The test family includes `path_3`, `path_4`, `path_5`, `star_4`, and `cycle_4`, with the declared spin coverage and exhaustive finite partition controls.

The reduced-boundary controls also verify equality of the reduced density matrices, their parameter derivatives, and SLD mixed QFI before and after exact fusion. The pure-minus-mixed QFI gap is positive semidefinite within numerical tolerance.

Important claim boundary:

- Internal edges are contracted exactly; this is not a trace over internal degrees of freedom.
- No physical spin sum or BQG history measure is derived by this operation.
- Collective parameter pullbacks are explicit restrictions of the microscopic parameter space; their spectra are diagnostics, not universal coarse-graining invariants.
- The Segre/Haar result concerns a deliberately specified independent local Haar ensemble, not the physical BQG ensemble.
- No continuum limit or physical RG flow follows from these finite identities.

The causal-to-PL correspondence and finite Feshbach controls provide finite graph/compression evidence. The quotient fibres are not generally equitable, so this is not an exact graph isomorphism. The resulting Galerkin/compressed operator is not yet a QFI Laplacian with conductances derived from the same physical BQG state and measure. That common-data bridge remains open.

## 4. W5 production and Lorentzian frontier

The W5 branch contains substantial test and instrumentation work, but the production evidence is not closed.

- The SU(2) pruning-coverage inventory run [43745695837](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43745695837) succeeded. It establishes source-level coverage inventory, not production/shadow equivalence or negligible truncation error.
- The direct all-(M) run [43876911486](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43876911486) passed its 20 direct-gate unit tests and 7 comparator unit tests, then the actual `sine-H_E` operator computation ended with exit code 143 after a runner shutdown signal. No complete numerical result was emitted.
- The production/shadow aggregate [43749858391](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43749858391) reported `INCOMPLETE_PRODUCTION_SHADOW_EQUIVALENCE`, `comparison_performed=false`, `input_columns_compared=0`, and `max_relative_l2_error=null`, because no shard artifacts were available. This is no comparison, not a measured numerical mismatch.
- The full shadow/no-prune matrix run [43743640638](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43743640638) was interrupted before a complete matrix result; the V11 common-carrier pilot [43745695902](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/43745695902) was likewise cancelled before its production computation.

Additional documented blockers remain:

- the default (J_{\max,2}=10) spin wall is below the conservative one-step adequacy bound (20); this invalidates it as a completeness guarantee, although it does not by itself prove a nonzero omitted physical amplitude;
- per-update threshold pruning has a synthetic coherent-sum counterexample, so production loss needs a global/error-controlled treatment;
- full Lorentzian W5 joint-Gram closure still lacks the complete raw (A/B) production amplitude package and cutoff-converged adjoint evidence;
- the external C8 clock remains a positive control, not a derived physical BQG clock.

The full physical causal cone, physical rigging map, physical Born measure, retarded response, and common photon/matter/gravity cone therefore remain open.

## 5. Repository history and integration safety

At the updated audit snapshot (work-branch commit `6ef322664c5ef64889d82100952ef789245b71b8`), the work branch and `main` remain historically diverged from common base `5dc5789a0511c7cae6183bd824939e5c20cbd798`:

- work branch: 1,166 commits ahead of `main`;
- `main`: 12 commits not contained in the work branch;
- direct recursive Git-tree comparison of file/blob entries: `main` has 724 files and the work branch has 871;
- the work tree adds 147 files (18 workflows, 39 scripts, 45 tests and 45 documentation files), removes 0 files, and has 33 shared files with different blob contents (4 workflows, 18 scripts, 1 test, 8 Markdown documents and 2 JSON ledgers).

This direct tree comparison corrects the earlier shorthand that described two documentation paths as main-only: both `BQG_L2_PRODUCTION_PRIMITIVE_LINEARITY_AUDIT_20261111.md` and `BQG_RESEARCH_CONTINUITY_CHECKPOINT_20261111.md` exist on both branches, but the work-branch blobs are expanded revisions (18,580 vs 9,047 characters and 40,850 vs 11,788 characters, respectively). There are no files present only on `main` in this tree snapshot. The 33 differing shared files and divergent commit histories still require a deliberate integration review; do not merge or force-update `main` automatically.

## 6. Audit changes made on the work branch

- Repaired the broken README link to a non-existent `BQG_FULL_SYSTEM_AUDIT_20261113.md`.
- Marked the 10 October `AUDIT_CURRENT.md` as a historical snapshot and updated the README status navigation.
- Corrected the stale `structural_candidate_closed=true` declaration in `OPEN_PROBLEMS.md` to match the canonical core ledger and documented the single core no-go.
- Added a current machine-ledger snapshot to `THEORY_STATUS.md`.
- Added regression tests for the fixed-pure-ancilla no-go and the bilinear refinement construction.
- Fixed the missing `sys` import in `q2_relational_higher_spin_edge_qfi_gate.py`.
- Added a branch-only integrity workflow for compile, ledger, documentation, and complete test-suite validation.
- Enabled the all-branch Python AST audit workflow to run from this research branch.

## 7. Priority order for remaining work

1. Obtain a complete, resumable W5 all-(M) production run with every shard artifact present; then perform the production/shadow comparison on the complete declared input set.
2. Resolve W5 pruning loss and spin-wall adequacy before making any zero-loss or completeness claim.
3. Complete the Lorentzian W5 raw-amplitude, adjoint, cutoff, and joint-Gram evidence.
4. Derive the physical history measure, rigging map, relational clock/response, and common causal cone from the actual BQG constraints.
5. Derive QFI edge conductances and the QFI Laplacian from the same physical state/measure used by the causal-to-PL construction.
6. Define an actual RG map that treats internal degrees of freedom by the required trace/sum/effective action, then test scale dependence across larger spins and graph families.
7. Repair or explicitly archive the two historical syntax-error branch families.

## Claim boundary

This audit establishes software/test status and records finite mathematical controls. It does not certify the full BQG theory, physical Lorentzian dynamics, a physical RG fixed point, a continuum limit, or experimental confirmation.
