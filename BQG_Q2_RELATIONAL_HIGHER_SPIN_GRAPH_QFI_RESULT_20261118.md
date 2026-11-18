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

## Extended topology and QFI spectral diagnostics (2026-11-18)

The direct-contraction fixture family now includes:
- `path_3`: 3 vertices, 2 edges, 8 boundary half-links;
- `path_4`: 4 vertices, 3 edges, 10 boundary half-links;
- `star_4`: 4 vertices, 3 edges, 10 boundary half-links;
- `cycle_4`: 4 vertices, 4 edges, 8 boundary half-links.

For each spin and trial, the gate now records the smallest and largest eigenvalues
of the complete edge-twist QFI, its numerical rank and nullity, and verifies
the rank-nullity partition against the parameter dimension (3E). The
eigenvalue threshold is relative to `max(1, max(abs(eigenvalues)))), with
factor (10^{-9}). These diagnostics allow degeneracies to be observed rather
than silently assuming full rank.

The expanded suite is configured for (j=1/2) and (j=1), with two random
intertwiner samples per graph and spin. The test checks the exact boundary-leg
counts and spectral bookkeeping. This is an added test specification; the
expanded test run and CI status are not yet confirmed.

**Scope distinction:** these are spectra of the (3E\times3E) edge-twist
pure-boundary-state QFI. They are not the spectrum of the full (12V\times12V)
half-link QFI and do not by themselves classify vertex-Gauss or edge-frame
null modes.

The higher-spin graph contraction validator was also hardened to reject
fractional, boolean, string, and otherwise non-integer endpoint/half-link
indices before tensor contraction. Regression cases were added for these
invalid inputs; they remain subject to the pending CI execution.



## Coarse-graining-oriented spectral diagnostics and path-5 extension (2026-11-18)

The finite graph family now includes `path_5` (5 vertices, 4 internal edges,
12 boundary half-links). Its direct dense boundary-state calculation is
restricted to (j\leq1): the output vector has dimension
((2j+1)^{12}), so extending this fixture to (j=3/2,2) would create
unnecessarily large dense state vectors. The (j=3/2,2) scans remain on
the existing fixtures with at most 10 boundary half-links.

For each graph/trial, the edge-twist QFI spectrum is normalized by the exact
isolated-link scale (4j(j+1)/3). The report now includes:
- normalized spectral moments of orders 1 through 4;
- mean empirical spectral quantiles at probabilities 0, 0.25, 0.5, 0.75, 1;
- normalized spectral variance;
- numerical rank fraction (rank divided by (3E));
- explicit graph-by-spin coverage, so skipped high-dimensional cases cannot
  be mistaken for tested cases.

These quantities are diagnostics for comparing finite fixtures of different
parameter dimensions. No invariance under graph blocking or coarse-graining is
assumed. In particular, the current QFI is a pure-boundary-state QFI for
independent first-endpoint edge twists, not a renormalized effective action
or a QFI after tracing/blocking boundary degrees of freedom. A genuine
coarse-graining map and its induced parameter pullback must be specified
before claiming spectral flow or continuum universality.

The path-5 test and new spectral summaries are added to the workflow suite;
execution remains pending until a matching Actions run is observed.
