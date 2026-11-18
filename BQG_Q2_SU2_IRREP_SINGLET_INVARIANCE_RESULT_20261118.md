# Arbitrary-irrep local SU(2) singlet invariance (2026-11-18)

## Construction

For spin \(j\in\{0,\tfrac12,1,\ldots\}\), dimension \(d=2j+1\), and
magnetic basis \(m=j,j-1,\ldots,-j\), define
\[
(C_j)_{m n}=\frac{(-1)^{j-m}}{\sqrt{2j+1}}\,\delta_{n,-m}.
\]
This is the normalized invariant bilinear tensor. It is antisymmetric for
half-integer \(j\) and symmetric for integer \(j\).

With Hermitian generators \(J_a\) in the spin-\(j\) representation,
\[
J_a C_j+C_jJ_a^T=0.
\]
Consequently, for \(U_j=\exp(-i\theta^aJ_a)\),
\[
U_j C_j U_j^T=C_j.
\]
Thus the common-endpoint cancellation underlying the two-sided singlet
quotient is not peculiar to the fundamental representation.

## Numerical control

The new module constructs \(J_x,J_y,J_z\) directly from the ladder
operators and constructs \(C_j\) without an external symbolic package.
It checks normalization, Lie-algebra invariance, finite group invariance,
and exchange parity for \(j=1/2,1,3/2,2,5/2,3\). The tests additionally
include \(j=4\) and invalid input.

Files:
- \`scripts/q2_su2_irrep_singlet_invariance_gate.py\`
- \`tests/test_q2_su2_irrep_singlet_invariance_gate.py\`
- \`BQG_Q2_SU2_IRREP_SINGLET_INVARIANCE_RESULT_20261118.md\`

The workflow runs these tests and exports
\`verification_results/Q2_SU2_IRREP_SINGLET_INVARIANCE.json\`.

## Scope

This proves the local invariant-tensor identity for finite-dimensional
SU(2) irreducible representations. It does not construct higher-spin
four-valent intertwiners, spin sums, graph contractions, higher-spin QFI,
or BQG history-weighted conductances. Those require a separate network-level
implementation and validation.
