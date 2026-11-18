# Higher-spin graph relational QFI (2026-11-18)

The higher-spin path calculation is extended to direct contractions on both
the three-node path and four-cycle fixtures. Each edge uses the normalized
spin-\(j\) invariant bilinear tensor, and each internal edge receives an
independent three-axis twist at its first endpoint. The complete boundary
pure-state QFI is computed, including all cross-edge blocks.

For the path, the diagonal edge block is
\[
F_{ee}=\frac{4j(j+1)}3 I_3,
\]
and the cross block is checked against the exact middle-node Casimir:
\[
F_{12}=\frac43\langle\mathbf J_{B0}\cdot\mathbf J_{B3}\rangle I_3.
\]
For the cycle, no path formula or edge-additivity assumption is imposed;
the full contracted-state QFI is computed directly and checked for symmetry,
PSD, and nonzero inter-edge correlations.

The \(j=1/2\) specialization is cross-checked against the existing graph QFI
implementation for both graph fixtures. Tests cover \(j=1/2\) and \(j=1\).

Files:
- \`scripts/q2_relational_higher_spin_graph_qfi_gate.py\`
- \`tests/test_q2_relational_higher_spin_graph_qfi_gate.py\`
- \`BQG_Q2_RELATIONAL_HIGHER_SPIN_GRAPH_QFI_RESULT_20261118.md\`

CI exports \`verification_results/Q2_RELATIONAL_HIGHER_SPIN_GRAPH_QFI.json\`.
This remains a finite pure-state graph computation, not a spin sum or
history-weighted continuum conductance.
