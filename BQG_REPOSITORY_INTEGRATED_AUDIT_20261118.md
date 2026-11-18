# Integrated BQG repository audit — 2026-11-18

## Executive result

The research branch is substantially healthier as a software artifact than the previous Q2-only CI view suggested: its Python sources compile, the canonical structural/physical ledgers validate, the documentation checks pass, the complete `tests/` suite has a green run, the Q2 export workflow produces valid artifacts, and the all-branch Python AST audit completes.

That does **not** mean the physical theory is closed. The machine-readable ledgers still report an open structural candidate and an open physicalization frontier. The W5 production calculations also remain incomplete: several expensive jobs were interrupted before emitting the artifacts needed for a production/shadow comparison.

Audited branch: `work/q2-two-sided-singlet-qfi-20261118`  
Main branch: not modified; no merge performed.

## 1. Reproducible software evidence

### 1.1 Full branch regression

Most recent completed full branch-integrity run before the graph-topology Gram addition: [run 44393440176](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44393440176), on commit `36336bd6d8f37ba0e94963eb14d0d8d406377552`.

- Python compilation: PASS.
- Canonical structural, depth-6 and physicalization ledgers: PASS.
- README/local-link validation: PASS.
- Full `pytest -q tests`: **509 passed, 59 subtests passed, 0 failures, 0 errors, 0 skips** (292.02 seconds).

After that run, the graph-topology fusion module and four regression cases were added. The focused Q2 workflow [44394058368](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44394058368) completed successfully: the existing vertex-fusion suite reports 36/36 passing, the new graph-Gram module reports 4/4 passing, all export steps succeeded, and the 26-file JSON artifact was uploaded.

A fresh full-repository regression run [44394789020](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44394789020), on commit `bdb3fd15a007f7e6894ea5fb9a6c43607604e543`, was in progress at the time of this audit update. Do not treat the earlier 509-test run as a full-suite verdict for the newly added graph-Gram module; that module is currently covered by the successful focused Q2 run.

### 1.2 Q2 graph/QFI export

Latest workflow: [q2-pl-graph-export, run 44394058368](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44394058368)  
Artifact: [q2-pl-graph-export](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44394058368/artifacts/14468838975)

The artifact contains **26 JSON exports**. I downloaded the archive and independently parsed every JSON file; all 26 are valid, every top-level `passed` field is true, and no nested `passed=false` record occurs.

The exhaustive reduced-boundary sweep reports:

- 144 boundary subsystems;
- 2,880 partition-order evaluations;
- maximum reduced-density error: `3.885780586188048e-16`;
- maximum reduced-derivative error: `1.4033044459902308e-16`;
- maximum mixed-QFI error: `6.106226635438361e-16`;
- minimum eigenvalue of pure-minus-mixed QFI: `0.0`.

The independent local-Haar pushforward control covers seven finite spin/block configurations. The product-Haar support is a Segre variety inside the fusion image; its measured complex projective codimensions are 1, 4, 11, 4, 20, 9 and 16. This establishes non-closure of that specified toy ensemble under the tested fusion maps, not of the physical BQG measure.

The new topology-dependent Gram controls are also in the artifact:
- `star_4` at j=1/2: rank 16/16, scalar Gram map, condition number 1;
- `cycle_4` at j=1/2: rank 12/16, rank defect 4;
- `cycle_4` at j=1: rank 61/81, rank defect 20.

All three report `passed=true`; boundary Gauss residuals remain below 3e-15. The complete finite result and its interpretation limits are recorded in [BQG_Q2_VERTEX_FUSION_GRAPH_GRAM_RESULT_20261119.md](BQG_Q2_VERTEX_FUSION_GRAPH_GRAM_RESULT_20261119.md).

### 1.3 All-branch Python static audit

Latest scan after adding the graph-Gram script and tests: [python-all-branches-audit, run 44393987545](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44393987545)  
Artifact: [python-all-branches-audit](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44393987545/artifacts/14469902365)

- 174 branch heads; 163 unique head commits;
- 22,811 Python-file occurrences across branch snapshots;
- 816 unique Python paths and 1,008 path/content versions;
- 1,005 unique Python blob contents; 1,566 sources parsed, including embedded and inline sources;
- 10 embedded archives scanned, 983 members and 335 embedded Python sources;
- 0 unsafe archive members and 0 archive extraction errors;
- `main_ast_status=PASS`, with 0 syntax errors;
- the work-branch snapshot containing the new graph-Gram module has no syntax-error finding.

Two syntax errors remain only in historical branch snapshots:

1. `.github/workflows/lorentzian-pruned-first-column.yml#heredoc-PY-2:86`, on two historical branch heads.
2. `scripts/pl_lorentzian_group_worker_v2.py:49`, on three historical branch heads.

The scan also reports one inline-Python syntax error across the historical snapshot corpus; the main branch has zero inline syntax errors. These are historical-source findings, not failures in the current work-branch Python source.

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

A graph-topology Gram follow-up is now implemented in `scripts/q2_vertex_fusion_graph_gram.py`, with report [BQG_Q2_VERTEX_FUSION_GRAPH_GRAM_RESULT_20261119.md](BQG_Q2_VERTEX_FUSION_GRAPH_GRAM_RESULT_20261119.md). For the four-node tree `star_4` at spin j=1/2, the fusion map has rank 16/16 and C†C = I/64 within 1.39e-17. For the four-cycle, the Gram map is non-scalar and rank-deficient: 12/16 at j=1/2 and 61/81 at j=1. Boundary Gauss residuals remain below 3e-15. This means a normalized loop-block ensemble needs state-dependent norm reweighting; the tree Haar-pushforward covariance formula cannot simply be reused. The focused module reports 4/4 passing, and Q2 workflow [44394058368](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44394058368) completed successfully, exporting 26 JSON controls including the three new Gram records.

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

At the updated audit snapshot (work-branch commit `bdb3fd15a007f7e6894ea5fb9a6c43607604e543`), the work branch and `main` remain historically diverged from common base `5dc5789a0511c7cae6183bd824939e5c20cbd798`:

- work branch: 1,173 commits ahead of `main`;
- `main`: 12 commits not contained in the work branch;
- direct recursive Git-tree comparison of file/blob entries: `main` has 724 files and the work branch has 874;
- the work tree adds 150 files (18 workflows, 40 scripts, 46 tests and 46 documentation files), removes 0 files, and has 33 shared files with different blob contents (4 workflows, 18 scripts, 1 test, 8 Markdown documents and 2 JSON ledgers).

This direct tree comparison corrects the earlier shorthand that described two documentation paths as main-only: both `BQG_L2_PRODUCTION_PRIMITIVE_LINEARITY_AUDIT_20261111.md` and `BQG_RESEARCH_CONTINUITY_CHECKPOINT_20261111.md` exist on both branches, but the work-branch blobs are expanded revisions (18,580 vs 9,047 characters and 40,850 vs 11,788 characters, respectively). There are no files present only on `main` in this tree snapshot. The 33 differing shared files and divergent commit histories still require a deliberate integration review; do not merge or force-update `main` automatically.

## 6. Audit changes made on the work branch

- Repaired the broken README link to a non-existent `BQG_FULL_SYSTEM_AUDIT_20261113.md`.
- Marked the 10 October `AUDIT_CURRENT.md` as a historical snapshot and updated the README status navigation.
- Corrected the stale `structural_candidate_closed=true` declaration in `OPEN_PROBLEMS.md` to match the canonical core ledger and documented the single core no-go.
- Added a current machine-ledger snapshot to `THEORY_STATUS.md`.
- Added regression tests for the fixed-pure-ancilla no-go and the bilinear refinement construction.
- Reconciled the README's gate registry count with the canonical ledger (56 total records; 50 core gates) and added the missing Lorentzian full-node conditional gate to the reader-facing atlas.
- Recomputed branch divergence from recursive Git trees at blob level; corrected the report to distinguish history divergence from actual missing files.
- Added and CI-wired the tree-versus-loop fusion Gram diagnostic, its four regression cases, three JSON exports, and a dedicated mathematical report.
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
