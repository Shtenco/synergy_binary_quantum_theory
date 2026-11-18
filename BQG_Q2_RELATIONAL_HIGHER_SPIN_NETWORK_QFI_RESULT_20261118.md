# Higher-spin four-valent path relational QFI (2026-11-18)

## Construction

For \(j\geq\tfrac12\), let \(d=2j+1\). The four-valent invariant subspace is
constructed numerically as the joint kernel of
\[
J^a_{\rm tot}=\sum_{\ell=0}^{3}
I^{\otimes\ell}\otimes J^a\otimes I^{\otimes(3-\ell)},
\qquad a=x,y,z.
\]
Its expected dimension is \(2j+1\), corresponding to the intermediate
recoupling spin \(k=0,1,\ldots,2j\). The code checks the nullspace residual.

The normalized invariant two-leg tensor is
\[
(C_j)_{mn}=\frac{(-1)^{j-m}}{\sqrt{2j+1}}\delta_{n,-m}.
\]
Two copies contract the A-B and B-C edges of a three-node path. The retained
boundary state is normalized and its full six-parameter pure-state QFI is
computed, retaining cross-edge blocks.

## Exact predictions

Schur's lemma gives the one-leg reduced density matrix \(I_d/d\), hence
\[
F^{(e)}_{ab}=4\,\mathrm{Cov}(J^a,J^b)
=\frac{4j(j+1)}{3}\delta_{ab}.
\]
For two links meeting at a four-valent middle node, SU(2) covariance gives
\[
F^{(12)}_{ab}
=\frac43\langle\iota_B|
\mathbf J_{B0}\cdot\mathbf J_{B3}|\iota_B\rangle\,\delta_{ab}.
\]
The pair-Casimir spectrum is
\[
\operatorname{spec}(\mathbf J_0\cdot\mathbf J_3)
=\left\{\frac12[k(k+1)-2j(j+1)]:
k=0,\ldots,2j\right\}.
\]
Thus network non-additivity persists at arbitrary tested spin; the
correlation coefficient is state-dependent and bounded by this exact
recoupling spectrum.

## Numerical gate

The new module uses only NumPy and the repository's Python standard-library
workflow. It tests \(j=1/2,1,3/2,2\), with:
- invariant-subspace dimension and Gauss residual;
- singlet normalization and Lie-algebra invariance;
- isotropic single-link QFI with the exact \(4j(j+1)/3\) scale;
- isotropic cross-link block and pair-Casimir spectral bounds;
- Gauss pullthrough on both shared middle-node legs;
- positivity of the complete 6x6 network QFI.

Files:
- \`scripts/q2_relational_higher_spin_network_qfi_gate.py\`
- \`tests/test_q2_relational_higher_spin_network_qfi_gate.py\`
- \`BQG_Q2_RELATIONAL_HIGHER_SPIN_NETWORK_QFI_RESULT_20261118.md\`

CI exports \`verification_results/Q2_RELATIONAL_HIGHER_SPIN_NETWORK_QFI.json\`.

## Scope

This is a finite pure-state three-node path calculation for individual
representations. It does not yet sum over spins, use BQG history amplitudes,
trace boundary legs, or derive continuum graph conductances.
