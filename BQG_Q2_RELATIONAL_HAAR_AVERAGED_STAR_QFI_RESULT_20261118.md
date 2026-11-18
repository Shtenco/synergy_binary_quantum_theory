# Haar-averaged four-valent spin QFI and parametric spin ensemble (2026-11-18)

## Exact average over pure intertwiners

For fixed \(j\), the four-valent invariant subspace has dimension \(d_I=2j+1\)
and projector \(P_j\). A Haar-distributed normalized pure intertwiner
\(|\iota\rangle\) obeys
\[
\mathbb E_{\rm Haar}|\iota\rangle\langle\iota|=P_j/d_I.
\]
Every invariant pure intertwiner has \(\langle\mathbf J_\ell\rangle=0\), so
the pure-state QFI is linear in the two-generator expectations:
\[
\mathbb E F_{(\ell a),(m b)}
=\frac4{d_I}\operatorname{Tr}_{\mathrm{Inv}}
(P_j J_\ell^aJ_m^bP_j).
\]

For one leg, Schur's lemma gives
\[
q_j=F_{\ell\ell}^{aa}=\frac{4j(j+1)}3.
\]
For distinct legs, recoupling the pair \((\ell,m)\) through \(k=0,\ldots,2j\)
gives eigenvalues
\[
\mathbf J_\ell\cdot\mathbf J_m
=\frac12[k(k+1)-2j(j+1)].
\]
Their uniform trace average is \(-j(j+1)/3\), hence
\[
\mathbb E F_{\ell m}^{ab}=-\frac{q_j}{3}\delta_{ab}.
\]
The exact \(12\times12\) ensemble-averaged QFI is therefore
\[
\boxed{\mathbb E F_j
=q_j\left(\frac43 I_4-\frac13\mathbf1\mathbf1^T\right)\otimes I_3.}
\]
Its spectrum is zero with multiplicity 3 (the nodewise Gauss rotations) and
\(4q_j/3\) with multiplicity 9. The averaged metric is positive definite on
the 9-dimensional quotient.

## Weighted spin sectors

For a specified classical ensemble of orthogonal spin sectors with
parameter-independent normalized weights \(p_j\), the sector-QFI average is
linear:
\[
\bar F=\sum_jp_j\mathbb E F_j,\qquad
\bar q=\sum_jp_j\frac{4j(j+1)}3.
\]
The result has the same star structure with \(q_j\) replaced by \(\bar q\).

The code includes an illustrative distribution \(p=(0.1,0.2,0.3,0.4)\) on
\(j=(1/2,1,3/2,2)\) solely to verify normalization, linearity and positivity.
These probabilities are not inferred from BQG dynamics.

## Important distinction

This is the average of pure-state QFIs over a classical Haar ensemble of
intertwiners. It is not generally equal to the QFI of the Haar-averaged mixed
density matrix; QFI is nonlinear under mixing. A physical spin sum must
specify whether intertwiner labels are classical and accessible, whether
spin sectors are orthogonal/superselected, and how BQG history amplitudes
determine their weights.

Files:
- \`scripts/q2_relational_haar_averaged_star_qfi_gate.py\`
- \`tests/test_q2_relational_haar_averaged_star_qfi_gate.py\`
- \`BQG_Q2_RELATIONAL_HAAR_AVERAGED_STAR_QFI_RESULT_20261118.md\`

CI exports \`verification_results/Q2_RELATIONAL_HAAR_AVERAGED_STAR_QFI.json\`.
