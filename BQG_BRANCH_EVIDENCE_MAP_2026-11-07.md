# BQG unified branch evidence map — 2026-11-07

## Snapshot and scope

Repository: `Shtenco/binary_quantum_theory`  
Canonical main: `90a5acde961c09cde565f942c9e21a7d7db1376b`  
Branch refs fetched from GitHub on 2026-11-07: **147** total, **146** non-main. Unique head commits: **135**. Refs sharing identical heads: 4 groups.

This is an evidence map, not a declaration that BQG is physically complete. The preceding audit structurally compared all 146 non-main refs with main. Large compare responses can truncate the per-file list at 300. Deep semantic review focused on the named scientific families and the V21/V22 execution path; it is not yet a line-by-line semantic audit of every file on every branch. No fresh complete local replay of the full physics pipeline was performed.

## Executive verdict

There are exact finite mathematical theorems, finite numerical computations, conditional IR results, and finite no-go/necessary-bound results. **There is no complete proof that BQG derives physical GR/Einstein dynamics from its binary microstructure.**

The nine mandatory physical gates in `physicalization_gates.json` remain `open_physical`: physical projector/history; connected inter-block history; physical TT kernel; six physical IR coefficients; common physical scale; dynamical Maxwell kernel; physical background cosmology; physical scalar cosmology; lensing/dynamics closure.

## Evidence taxonomy

- **proved_scoped**: theorem is exact only under its stated construction and hypotheses.
- **tested_finite / computed_saved_data**: finite carrier or saved calculation; no automatic continuum inference.
- **conditional**: conclusion follows if hypotheses hold; BQG derivation of those hypotheses remains open.
- **proved_conditional_finite**: exact finite necessary bound under declared common-space assumptions.
- **defect_confirmed_code_review**: source-code acceptance/provenance flaw; not by itself a refutation of the mathematical model.
- **open_physical / incomplete**: required physical bridge or numerical witness absent.

## Core evidence map

| ID | Result | Status | Evidence | Scope / limitation |
|---|---|---|---|---|
| BQG-DIM3-001 | Frozen q=2 refinement fixed point (d_*=3) | Scoped theorem | `main: theory_gates.json`, dimension-3 closure gate | Applies to the specified refinement rule, not all possible information systems |
| BQG-MAN3-001 | Selected q=2 PL completion and manifold link checks | Finite tested | `main: GLOBAL_MANIFOLD_Q2_COMPLETION.md` | Selected PL completion and tested refinements |
| BQG-WALSH-001 | Walsh tetrahedral flux frame and finite Gauss-singlet geometry | Scoped theorem | `main: MICRO_WALSH_QGEOM_BRIDGE.md` | Finite algebraic bridge, not physical dynamics |
| BQG-FESHBACH-001 | Exact finite Feshbach/Schur and block-Krylov identities | Scoped theorem | `main: FESHBACH_INTERBLOCK_EFFECTIVE_KERNEL.md` | Constraint spectral parameter is not physical frequency |
| BQG-MASTER-001 | `ker(C†GC)=intersection ker(C_A)` for positive finite constraint metric | Scoped theorem | `main: MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md` | Does not construct continuum rigging map |
| BQG-L1-024 | First-refinement L1 rank 24 | Finite theorem | `research/physicalization-l1-depth2-metric` | Not 24 established physical graviton modes |
| BQG-METRIC-CAL-001 | Metric-map singular values (1/\sqrt6,1/\sqrt6,1/\sqrt3), condition number \(\sqrt2\) | Finite theorem | `research/physicalization-l1-depth2-metric` | Physical kinetic ratio not yet measured |
| BQG-IR-UNIV-001 | ADM/Einstein IR universality under H1–H4 | Conditional theorem | `research/bqg-v16-first-leg-20261107` | H1–H4 not derived from BQG |
| BQG-CARRY-001 | Rank-one micro-carry identity and one-channel lower bound | Finite theorem | `research/bqg-micro-carry-rank-20261031` | Kinematic, not physical history |
| BQG-HIST-001 | (D_\infty=\mathbb Z\rtimes\mathbb Z_2) history/orientation algebra | Scoped algebraic theorem | `research/q2-history-gravity-orientation-current` | Not a derivation of physical time |
| BQG-W5-GRAM-001 | W5 Euclidean Gram saved matrix: Gershgorin lower bound ≈30.9860; diagonalized minimum ≈32.4136 | Saved-data computation | `research/bqg-node1-lorentzian-columns-20261031` | K5 W5 only; printed-decimal bound does not enclose SU(2) numerical errors |
| BQG-HDA-NEC-001 | Ten-channel PL16 HDA compensation rank/norm lower bound | Conditional finite bound | `research/bqg-v16-first-leg-20261107`, V19 report | Actual K/D not computed; not global no-go |
| BQG-HDA-Q-NEC-001 | K5 W5 common-kernel compensation norm lower bound | Necessary condition | Same V19 report | Actual Q not computed; K5 W5 only |
| BQG-V21-GUARD-001 | Generic stage diagnostics, empty-state acceptance, trusted guard flag | Code defect | V21 guard/collector/tests | Blocks use as physical acceptance certificate |
| BQG-V21-INPUT-001 | Input provenance does not encode direct per-stage input checkpoint digest | Code/provenance defect | V20 worker + V21 collector | NPZ itself is hash-checked, but chain attestation incomplete |
| BQG-V21-ZERO-001 | Exact-zero V4/V2 can compare as relative error zero; final pair nonzero not required | Semantics defect | V21 collector | Exact zero may be valid but must be classified separately |
| BQG-PHYS-GATES-001 | Nine physical gates open | Open physical | `main: physicalization_gates.json` | No physical completion claim |
| BQG-DEPTH6-001 | Three depth-6 numerical witnesses pending; finite target kernel theorem not proved | Incomplete | `main: depth6_frontier.json` | Structural support closure ≠ numeric master-kernel closure |

## V21/V22 code audit

Reviewed files on `research/bqg-v21-c1-checkpoint-repair-20261107`:
- `scripts/v21_checkpoint_guard.py`
- `scripts/v21_stage_logic.py`
- `scripts/pl_v21_pair_collector.py`
- `scripts/pl_v20_stage_worker.py`
- `scripts/pl_lorentzian_triple_worker.py`
- `scripts/pl_lorentzian_physical_pair_worker_v4.py`
- `tests/test_v21_c1_checkpoint_repair.py`

Findings:
1. Guard uses one generic diagnostic key set for all stages and requires only one relevant key.
2. Empty stage state can pass the guard; tests explicitly encode that.
3. Collector trusts `meta.guard.stage_physical_basis_pass` instead of recomputing from original diagnostics.
4. Stage metadata carries the original seed input hash through later stages; it does not provide a direct, independently checked input-checkpoint digest chain.
5. V4/V2 exact-zero states can yield zero relative mismatch. Final pair norm is recorded but nonzero is not a condition for a physical signal.
6. V21/V22 have identical head SHA `2c28123720db59306890daa0dfba1ad21c2e7eb1`; V22 is not an independent code/science increment.

Required repair:
- Stage- and operation-specific mandatory diagnostic schemas.
- Collector recomputes guard from raw diagnostics; never trusts a boolean copied from metadata.
- Direct input NPZ SHA, output NPZ SHA, metadata SHA, backend SHA, source commit, and stage signature bound in a manifest.
- Separate `EXACT_ZERO`, `VERIFIED_NONZERO`, and `INVALID/INCOMPLETE` outcomes.
- Add adversarial tests for omitted diagnostics, wrong-stage keys, modified JSON, broken predecessor links, stale backend hashes, empty state, and zero-pair acceptance.
- Only after independent replay of all 32 paths and four independent V2 references may C1 be considered complete.

## Carrier / graph separation — do not combine

- **K5 W5 Euclidean Gram** and **PL16 ten-channel HDA probes** are different carriers and output spaces. Their matrices, spectra, ranks, and norm bounds must never be concatenated or compared as if they were one operator.
- **W5 Euclidean calculation** is not a Lorentzian HDA calculation.
- **Saved-data replay/post-processing** (V17–V19) is not a new SU(2) amplitude calculation.
- **Heat/master spectral parameter** is not physical time or physical frequency.
- **Finite C8 Page–Wootters clock** is a positive control, not a clock derived from BQG gravitational dynamics.
- **Internal SU(2) spin** is not automatically spacetime spin-2.
- **L1 rank-24 subspace** is not automatically a physical TT sector.
- **Kinematic graph-changing isometry** is not a physical rigging map.

## Superseded / historical results

- Old `FULL_OPERATOR_FIRST_HDA_CERTIFICATE.md` using (G_v=-\frac23E_v+\frac{32i}{9}L_{raw,v}) is superseded/conditional because (L_{raw}) is not globally anti-Hermitian. Correct Hermitian structural block uses (S_v=-\frac i2(L_{raw}-L_{raw}^\dagger)), then (G_v=-\frac23E_v-\frac{32}{9}S_v).
- `release/canonical-theory-package-v1` is a historical snapshot; older closure language does not override current physical gates.
- `research/bqg-v22-local-c1-e2e-20261107` duplicates V21 head.
- `research/v19-canonical-k5-replay-20261015` and `research/v20-160-canonical-20261015` share a head SHA.
- Ten `research/bqg-2t-*` refs share one head SHA; count as one line of evidence, not ten confirmations.
- `cleanup/gr-qm-bit-continuum` is historical background; preserve it, do not wholesale merge.
- Workflow failures with zero job steps are infrastructure/admission failures, not physical no-go results. Exact cause remains unproven.

## Branch inventory and family disposition

The machine-readable registry `BQG_PROOF_REGISTRY_2026-11-07.json` records all 147 branch names and head SHAs, same-head aliases, family labels, evidence IDs, and dispositions. Branch family labels are triage categories, not a claim that every file in that branch has been read.

## Next gates, in order

1. Repair V21 acceptance/provenance and zero-result semantics.
2. Demonstrate runner admission before launching expensive physical paths.
3. Complete 32 genuine V4 paths and compare against four independent V2 primitives.
4. Build same-carrier full HDA (E,S,R,D) and evaluate the necessary bounds on matching spaces.
5. Construct nonzero stable physical projector/history across refinement.
6. Derive connected correlators and interacting physical TT pole.
7. Only then test the H1–H4 hypotheses required by the conditional IR universality theorem.

**Conclusion:** BQG has substantial finite mathematical content, but the complete physical derivation of GR remains open. Preserve negative results and scoped theorems; do not inflate a finite or conditional result into a physical proof.


## V21 repair specification

A separate implementation-ready audit document is included: `V21_C1_ACCEPTANCE_PROVENANCE_REPAIR_SPEC.md`. It defines operation-specific diagnostic schemas, direct checkpoint digest chaining, trusted manifest requirements, zero-result semantics, adversarial regression cases, and strict C1 closure criteria. It is a specification only; no V21 worker code has been modified.

## Conditional IR gate interpretation — clarified

The root-level theorem is `BCQG_IR_UNIVERSALITY_CLOSURE_THEOREM.md` (blob SHA `0cf465e954300d452d3d7bd2c11f7304086aa560`). The associated gate script `scripts/bcqg_ir_universality_closure_gate.py` explicitly says it checks theorem-level reductions and frozen finite prerequisites only. Its stored result `verification_results/BCQG_IR_UNIVERSALITY_CLOSURE.json` has `passed=true`, but lists H1–H4 as open. In particular, the gate sets the algebraic target values and rank tuple to verify the conditional reductions; it does not measure the collective HDA or prove those ranks from BQG microdynamics. Therefore its status remains **conditional theorem + finite algebraic/prerequisite checks**, not physical IR closure.
