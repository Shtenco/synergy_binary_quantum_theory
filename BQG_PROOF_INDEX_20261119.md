# BQG proof and evidence index — 2026-11-19

**Integration base:** `main` at `401b006d368caf4bfd7578c0d312795cf2f2323a`  
**Integration branch:** `integration/bqg-proof-index-mirror-spectral-20261119`  
**Merge policy:** selective file-level integration only; no wholesale merge of divergent research branches.

This index is a provenance map. A source blob SHA identifies exact file bytes; a commit SHA identifies the source revision where a result was produced or the integration commit that imported it. A CI PASS is scoped to its exact head, workflow and finite test domain.

## 1. HDA / collective gravity

| Evidence | Source ref / blob | Exact provenance | Status and boundary |
|---|---|---|---|
| Collective HDA intertwining theorem | `research/q2-history-gravity-orientation-current`; blob `02cd270fcfa97ef2890ae37581853fb96bc6f0a7` | Imported as `COLLECTIVE_HDA_INTERTWINING_THEOREM.md` | Exact finite-dimensional compression identity; conditional collective-HDA inheritance via leakage and structure-function defect. Production collective residual is not yet measured. |
| IR universality closure theorem | same ref; blob `0cf465e954300d452d3d7bd2c11f7304086aa560` | Imported as `BCQG_IR_UNIVERSALITY_CLOSURE_THEOREM.md` | Conditional on H1–H4: regular 3D metric, local two-derivative ADM form, first-class HDA, no second-class remainder. Not an unconditional derivation of GR. |
| Physical sine two-node HDA | `research/bcqg-core-candidate-v1`; blob `9e6743a1c048ff01562779f30545e52f26bca051` | Workflow `two-node-sine-hda`, run `31855735615`, tested head `c02c1892cd9d18b9a62413d6c06abf1211e6b7f6`; artifact id `9239198320`, recorded SHA256 `21e2da508fd583d9007a5bd400d074e8cee39990656e6c75e5968d2601323526` | Finite preregistered PASS: p_cross=1.0056949, p_GG=2.0074904, p_joint=1.0076444; endpoint 0.02003034. Artifact endpoint currently returns no downloadable artifact; provenance remains in the result document. |
| Full signed operator-first HDA certificate | `research/bcqg-core-candidate-v1`; blob `6aa8f2dc4f067949d66daf6601c516158950eea7` | Spin-changed route workflow run `31858615323`, tested head `905358f76699370aa13e017bd852c45b696b3e5f`; artifact id `9244277324`, recorded SHA256 `c1af8de00183fddf328f6bdfba386e2320b842e10d3de98d90ad150b0876213c` | Conditional asymptotic operator-first certificate plus finite spin-changed route evidence. EL/LE/LL full signed channel calculation was cancelled before artifacts; not a complete finite production PASS. |
| HDA compression regression added in this integration | `tests/test_bqg_hda_compression_identity.py` | Integration commit `36dc99f72c9c548eaed9ffb5a5a586377b6169c2` | New deterministic NumPy tests for exact compression identity and leakage bound. Tests algebraic infrastructure, not the full BQG production bracket. |

**HDA branch comparison:** `research/bcqg-core-candidate-v1` was 411 commits ahead / 632 behind main; `research/q2-history-gravity-orientation-current` was 422 ahead / 632 behind at comparison. These branches were not merged wholesale. Four high-value theorem/result documents were copied individually.

## 2. Mirror mass theorem and Wilson–Dirac carrier

| Evidence | Source blob / commit | Integration | Status and boundary |
|---|---|---|---|
| General mirror mass-derivative theorem | blob `0da78562204a5d196d39a32b19be6c5ae7e95919` | `MIRROR_MATTER_MATRIX_ELEMENT.md` | Mirror covariance implies equal aligned partner masses and opposite static sigma derivatives. The actual microscopic beta_m remains to be derived. |
| Wilson–Dirac carrier | blob `afbd2d5e43f59086b8c75d7da73fc4678da3408c`; source commit `743abc6fba0f1c11cbe55c17ebed119c9ed5f6dc` | `MIRROR_WILSON_MATTER.md` | Finite regular-lattice carrier, not realistic irregular PL/Peter–Weyl chiral matter. |
| Executable Wilson gate | blob `6ba84564cc9bb00767da23a36475d988f1db2c2b`; source commit `c703b7ce1059956f93baff68a2fa56ddb0ad7060` | `scripts/mirror_wilson_matter_gate.py` | Checks mirror spectra, positive masses, opposite Hellmann–Feynman sources, moving-state derivative and one physical Wilson corner. beta=0.37 is an input, not a derived constant. |
| Matrix-element gate | blob `561037f46852ad828e3cb32b4bbd0b1c635aeceb` | `scripts/mirror_matter_matrix_element_gate.py` | Separates axial density (not a universal static charge) from the mirror-odd rest-mass derivative. |
| New integration tests/workflow | integration commits `c20b2c8ad16c9c8d909700440c3210912e50bacc`, `6ba8c8d20b5f6bb29703f8fbea52abef0136fca9` | `tests/test_mirror_wilson_matter_gate.py`, `.github/workflows/bqg-mirror-wilson-dirac.yml` | New PR CI will run both executable gates and upload JSON evidence. |

Source mirror branch was 3 commits ahead / 633 behind main. Only the two matter-gate documents and their two small executable gates were selected; unrelated mirror-orientation work was not merged.

## 3. Spectral-history

| Evidence | Source blob | Integration | Status and boundary |
|---|---|---|---|
| Operator-derived spectral graph theorem | `d02d4a83bba4c9a6cee93ad29754535fdb781b67` | `BQG_SPECTRAL_HISTORY_GRAPH.md` | Exact finite block-Lanczos/Jacobi closure theorem and actual finite Euclidean first-edge protocol. Constraint spectral dimension is not automatically physical spacetime dimension. |
| Spectral graph gate | `735a0f037cacce3b7037a164313500ed5b8240bc` | `scripts/bqg_spectral_history_graph_gate.py` | Regression mode is deterministic and synthetic; actual Peter–Weyl mode is a distinct heavier calculation and is not run by the new regression workflow. |
| B3 preregistration | `e32f873a3b9ddfed8120a5177cc5a734d96a9825` | `PETER_WEYL_SPECTRAL_HISTORY_B3_PREREGISTRATION.md` | Next finite Euclidean shell test; no B3 closure outcome assumed. |
| Claim/provenance audit | `1a75d709dccd5fe2eeb909d93fdb772a4ec803c0` | `REPOSITORY_AUDIT_2026-09-06_SPECTRAL_HISTORY.md` | Correctly distinguishes q=2 route-family assumptions, selected PL S3 completion, finite constraint spectrum and open physical history. |
| New integration regression | integration commits `8e5ed5ba4b1479362ba1da1603bbed9719d17e50`, `f4f7df5dbfcfbebfd2da547665487d9df2a9d44f` | `tests/test_bqg_spectral_history_regression.py`, `.github/workflows/bqg-spectral-history-algebra.yml` | CI runs only the deterministic algebraic block-Lanczos regression; it does not claim the expensive actual K5 first-edge calculation passed. |

The source branch `spectral-history-closure` was 13 commits ahead / 591 behind main at comparison. The B3 workflow from that branch was deliberately not copied: it downloads source artifacts from the separate `binary_quantum_theory` repository and has external run dependencies. Its preregistration is preserved, but execution remains a separate task.

The integration also imports the actual sharded K5 master first-edge executor: source blob `c30b4c9274d3755fb59153f3662e1bc1b0b44557`, imported as `scripts/bqg_euclidean_master_first_edge_sharded.py` (integration commit `491b44416f2177091bd81c4b939584e81cb8d82c`). New workflow `.github/workflows/bqg-spectral-history-actual-k5.yml` (integration commit `0196b7da0b90f2f4ff6b94f368adcd8526c66322`) computes all 32 exact `M_E|i>` columns, assembles `A0`, `mu2`, `R1†R1=B1†B1`, checks the recurrence/PSD/cutoff guards, and reports whether the finite Euclidean master history is open or closed after its first edge. A PASS of the algebraic guards does not require B1 rank zero; the scientific outcome is reported separately.

The B3 preregistration's external source artifacts were checked directly: run IDs `31844567559` and `31852849936` in `Shtenco/binary_quantum_theory` currently return no downloadable artifacts. Therefore the B3 extension cannot be reproduced from those archived CI downloads at present. Do not run or claim B3 closure until the exact source columns are restored or regenerated with new pinned provenance.

## 4. Physical history and causal cone

| Evidence | Source blob | Integration | Status and boundary |
|---|---|---|---|
| Physical probability/causality frontier | `3302f4e16c181dde8a6f07dd660a5b537331c906` | `BQG_PHYSICAL_PROBABILITY_AND_CAUSALITY_FRONTIER_20261112.md` | Conditional Gram/Born/decoherence mathematics; physical inner product, clock and retarded kernels are not supplied by the theory yet. |
| Physical causal-cone audit | `0f210379029fc7795da25d73fe162dcccb9b9a81` | `BQG_PHYSICAL_CAUSAL_CONE_AUDIT_20261111.md` | Distinguishes graph support, relational order and physical retarded cone. |
| Lemaitre reference benchmark | `bc780eee16b2ddc5f0f955c745fa334f7905ba64` | `LEMAITRE_BQG_CAUSAL_CONE_BENCHMARK.md` | External finite reference benchmark, not a physical BQG light-cone proof. |
| Reference engine and tests | engine blob `f5e7d6ae13ad03e956686b8783ec0d637938aa8d`; tests blob `a2ebf1100aca8bbfb656dc55eae93f851b9ec579` | `scripts/bqg_physical_history_decoherence_engine.py`, `tests/test_bqg_physical_history_decoherence_engine.py` | Finite reference-only diagnostics for Gram positivity, decoherence, coherent support and sampled channel fronts. |
| New CI | integration commit `4c7f1e08a00a8101cb167c197a2d02f89f9a3bb0` | `.github/workflows/bqg-physical-history-reference.yml` | Explicitly enforces `PASS_REFERENCE_ONLY` and `OPEN_PHYSICAL`; cannot promote finite diagnostics to physical closure. |

The physical-history branch was 378 commits ahead / 26 behind main at comparison. Only the physical-history reference engine, its tests and claim-boundary documents were selected; no full 378-commit merge was attempted.

## 5. W5 V09/V10 live execution

- PR: [#86](https://github.com/Shtenco/synergy_binary_quantum_theory/pull/86), head `bae033314cd682d307e460a8e99fc3bf2c1aa31c`.
- Frozen V09 original commit: `6ae145ae61a19ed5dfb5a7dae5abf8d79d79eff9`; frozen file blob: `a2cc3d5eaee8377e94fd7d24a1a279d2399df133`.
- Workflow run: `44482668149`.
- In run `44482668149`, the two Jmax2=10 shards (`shadow_no_prune` and `production`) completed their unit-test stages (5 paired tests, 5 V09 adapter tests, 23 V10 adjoint tests) but each runner received a shutdown signal during the heavy numerical calculation. The remaining four jobs stayed queued; no JSON artifacts were uploaded. This is **INCOMPLETE / runner interruption**, not a mathematical FAIL and not a PASS.
- The paired harness was then hardened on PR #86: it now shards by one frozen seed K per job (30 jobs = 5 seeds × 3 cutoffs × 2 pruning modes), captures all 24 complex epsilon-term maps before and after signed accumulation, stores raw and serialized Gauss coefficients, verifies term-map reconstruction and the exact `1e-13` serializer cutoff, and aggregates all 30 reports fail-closed. New run `44496003089` was created for head `303cea6d291506041db602262b25996b1459bf68`; at last inspection all 30 jobs were queued and no artifacts existed yet.
- The old run's failed jobs cannot be restarted while that workflow remains active. The new seed-sharded run is the current execution path; do not call it PASS until all 30 reports and the aggregate JSON are inspected.

## 6. Merge and claim rules

1. Keep every source branch and its original commit history intact.
2. Integrate only files whose dependencies and scope have been checked.
3. CI status is valid only for the exact tested head SHA; never transfer an older green run to a newer PR head.
4. Distinguish `proved`, `tested_finite`, `conditional`, `reference_only`, `open_physical`, `incomplete` and `experimental`.
5. No BIP39/seed encoding step is evidence for a physical holographic principle. Such a connection needs a separately defined physical observable and falsifiable prediction.
6. Do not merge this integration branch into main until all five focused workflows are green on the final head and the PR diff is reviewed. The physical gates remain fail-closed.


## 7. Current integration PR and CI observation

- Selective integration PR: [#87](https://github.com/Shtenco/synergy_binary_quantum_theory/pull/87).
- PR head immediately before this index-status update: `8f72ef584d8d13d31023657364fb895fee853d67`; the status update itself is commit `0ac49d573cdbc966feb183900c1e3595d3d5592c`. Base `main` remains `401b006d368caf4bfd7578c0d312795cf2f2323a`.
- Five focused workflow definitions are present. At the audited head `2cca65c400dedd56b6be0a14a4263bb1f892ab8f`, HDA, mirror/Wilson, spectral-history algebra and physical-history reference completed successfully (runs `44496628707`, `44496629155`, `44496628696`, `44496628597`). Their finite scopes and artifacts were inspected. The actual 32-column K5 first-edge workflow `44496628504` remained queued with no column artifacts. The head has since received documentation/workflow-trigger corrections; all five checks must be observed again on the resulting final head.
- PR #87 remains Draft and must not be merged until the four workflows have run and their JSON artifacts have been inspected.


## 8. Current PR #87 head

The selective integration branch contains five focused workflows. At head `2cca65c400dedd56b6be0a14a4263bb1f892ab8f`, four completed successfully: HDA (2 tests), mirror/Wilson (both JSON gates passed), synthetic spectral regression (7-dimensional closure), and physical-history reference (`PASS_REFERENCE_ONLY`, with `OPEN_PHYSICAL`). Actual K5 run `44496628504` was queued and had no column artifacts at the last poll. Because the documentation and K5 dependency-trigger filters are being corrected, the final head must receive fresh runs before merge. W5 is separate: run `44496003089` still has 30 queued jobs and no artifacts; status remains `INCOMPLETE`.
