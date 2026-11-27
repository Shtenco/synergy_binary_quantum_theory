# BQG V23 — unified proof bridge and finite K5 evidence

**16.10.2026. Scope: research audit and reproducible finite calculations; NOT a full physical BQG proof.**

## Theorems already present in the repository

- MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md — exact finite positive master kernel/projector identity.
- BQG_REFINEMENT_PROJECTOR_STABILITY_THEOREM.md — projector-intertwining bound from master residual / minimum master gap; summable residual-to-gap ratios imply convergence of embedded projectors.
- BQG_CANONICAL_REFINEMENT_RESIDUAL_DECOMPOSITION.md and BQG_CONSTRAINT_LEVEL_REFINEMENT_BOUND.md — exact single-hit residual identities. These are the efficient alternative to blind K5 720-shell enumeration when studying *cross-scale* refinement.
- BQG_MICROSCOPIC_BILINEAR_REFINEMENT_THEOREM.md — exact bilinear ancilla refinement; a fixed pure singlet ancilla is an invalid symmetric shortcut.
- BQG_SERIAL_MULTIPLICITY_MASTER_RESULT.md — measured selection-channel gaps gamma_j for j=2,5/2,3,7/2,4. These are **not** full physical master gaps Delta_j.
- Q2_RELATIONAL_HISTORY_PROJECTOR.md — external-clock finite positive control, NOT physical time derived from the full BQG projector.
- TT_PROPAGATOR_FIRST_PASS.md — reduced free massless TT positive control, NOT interacting physical BQG 1PI.
- research/audit_2026_10_10/CONDITIONAL_PROOF.md — Einstein–Hilbert endpoint conditional upon physical Lorentzian metric, locality, dynamics and IR estimates.

## Numerical verification using saved canonical K5 output data from V22

Data source: BQG_V22_K5_Nested_Kernel_Audit.zip, its 101×101 V21_RESTRICTED_MASTER.npz and accompanying V22_STATUS.json. Original K5 backend SHA-256 596c4179cb2c8af20b41c18b3ab87b711c3850c5c5e05b1bc91a1f87c9f11cc4. All the listed calculations recompute the **finite matrix**, not the operator columns.

Results: 32 original even + 69 completed first-shell odd input states, full rank 101/101; lowest Gram eigenvalue 5.883687739952585; highest 14.483850719096251; condition number 2.4616960; forbidden even-to-odd master block norm 0.

For nested input dimensions 32+n, with n = 0,1,4,8,16,24,32,41,49,65,69, we recomputed R=M_big*I-I*M_small and the minimum positive Gram gaps. The **sum of ratios over only these ten finite inclusion steps** is 0.7560569651310478. It is not an asymptotic theorem, not a physical refinement sequence, and not a replacement for cross-scale BQG residual data.

**Key nontrivial insight:** for the initial 32-dimensional even sector embedded into the 101-dimensional even+odd space, R=0 because the mixed parity block vanishes, and both exact-zero projectors are identically zero. Hence exact projector consistency can be **vacuous**. A nonzero normalized anchored physical state or controlled positive rigging-map limit is an independent mathematical obligation.

## Physical witness gates still outstanding

1. One full H_E+H_L, Gauss+diffeo+HDA constraint family with matched domains, metric and ordering.
2. Verified microscopic refinement isometry in the full graph-changing habitat including new ancillas.
3. Cross-scale R_j, or one-hit L_j/D_j norms, with rigorous errors for this same family.
4. Positive **full master** gaps Delta_j (not multiplicity-channel gamma_j).
5. Summable epsilon_j and a nonzero anchored physical projector or normalized rigging-map form.
6. Connected physical boundary/history generating functional Z[J] from that same projector.
7. Positive interacting frequency-dependent TT pole and Ward relations extracted from that history.
8. Controlled nonlinear/local Einstein–Hilbert IR remainder.

The existing theorem proofs retain their valid scopes. The above physical witness gates are **not** automatically passed by a completed Euclidean K5 matrix, positive controls, reference propagator, or status flag.

## Run

Extract the separate full V23 ZIP provided in this conversation. It contains input NPZ and JSON, two standalone Python auditors, 17 passing tests, numerical CSV/PNG, a Russian README and SHA256 manifest.

From this GitHub checkout, after extracting V23/data files:

    python research/v23_unified_proof/check_v22_projection.py --v22-data /PATH/TO/V23/data --out /tmp/bqg_v23.json

The full package also contains verify_refinement_family.py, to evaluate actual future full-BQG scale matrices and embedding operators; it checks SHA256, isometry, gaps, projector residual bound and anchor, while refusing to pronounce the full physical theory proven without independently validated analytical and physical evidence.

This is a draft research continuation. No modification to main.
