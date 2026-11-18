# Higher-spin two-node relational QFI (2026-11-18)

## Analytic prediction

For an invariant four-valent node in the irreducible representation \(V_j\),
Schur's lemma implies that the reduced density matrix of any one leg is
\[
\rho_{\rm leg}=\frac{I_{2j+1}}{2j+1}.
\]
Consequently
\[
\langle J_a\rangle=0,\qquad
\frac12\langle\{J_a,J_b\}\rangle
=\frac{j(j+1)}{3}\delta_{ab}.
\]
For a normalized two-node gluing through the normalized spin-\(j\) singlet,
the contraction probability is
\[
p_{\rm glue}=\frac{1}{(2j+1)^2}.
\]
The pure-state relational QFI for a twist at either endpoint is therefore
\[
\boxed{F_{ab}=\frac{4j(j+1)}{3}\delta_{ab}.}
\]
The simultaneous endpoint action remains a singlet stabilizer, so
\(D_v=-D_u\).

For \(j=1/2\), this reproduces \(p=1/4\) and \(F=I_3\), the previously
derived two-block result.

## Explicit finite-dimensional implementation

The script constructs the four-valent invariant subspace as the common
kernel of the three total generators \(\sum_{\ell=0}^3J_a^{(\ell)}\).
Its expected dimension is \(2j+1\). Random normalized complex states are
sampled from this subspace, two nodes are glued through \(C_j\), and the
normalized six-leg boundary state and its three endpoint derivatives are
computed directly.

The gate checks:
- Gauss invariance and orthonormality of the invariant basis;
- gluing probability \(1/(2j+1)^2\);
- both endpoint QFI matrices against \(4j(j+1)I_3/3\);
- opposite endpoint projective derivatives;
- positive definiteness and JSON serialization.

Files:
- \`scripts/q2_relational_higher_spin_edge_qfi_gate.py\`
- \`tests/test_q2_relational_higher_spin_edge_qfi_gate.py\`
- \`BQG_Q2_RELATIONAL_HIGHER_SPIN_EDGE_QFI_RESULT_20261118.md\`

The workflow runs the tests and exports
\`verification_results/Q2_RELATIONAL_HIGHER_SPIN_EDGE_QFI.json\`.

## Scope

This is a higher-spin two-node channel, not a multi-edge higher-spin graph
calculation. It does not yet sum over spins, include spin-foam amplitudes,
or derive a causal-PL coarse-grained physical conductance.
