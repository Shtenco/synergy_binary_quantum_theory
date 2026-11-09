# Adjoint audit of the epsilon-oriented K-K-V block

Date: 2026-11-10

## Executive result

The raw 24-term epsilon-oriented K-K-V block is **not proven anti-Hermitian or Hermitian by the available leg identities alone**. A tempting proof by cyclically rotating the auxiliary trace is invalid: the trace is only over auxiliary indices, while its entries are operators on the gravitational Hilbert space and generally do not commute.

This note records the exact adjoint identity, explains the failed shortcut, and defines the numerical test required before any Hermitian completion is chosen. The earlier claim that the full epsilon sum is automatically anti-Hermitian and that multiplying it by \\(i\\) is therefore justified is withdrawn.

## Purpose and scope

This note closes one algebraic gap in the W5 Lorentzian pipeline. It derives the adjoint of the complete 24-term epsilon-oriented K-K-V contraction. It does **not** claim that the current numerical runner has already computed the full W5 matrix, that all truncations are converged, or that the physical normalization has been fixed.

The statement is conditional on the adjoint identities of the implemented covariant legs holding on one common, sufficiently complete Hilbert space.

## Definitions

At a fixed source vertex, let the matrix-covariant legs obey

\[
C_e(K)^\dagger=-C_e(K),\qquad C_e(V)^\dagger=C_e(V).
\]

The first identity follows from the current sine convention
\[
H_E^{\rm sine}=\frac{T-T^\dagger}{2i},\qquad
K_{\rm sine}=[V,H_E^{\rm sine}],
\]
provided the volume is self-adjoint and the commutator is represented without an adjoint-breaking projection. The second follows from self-adjointness of the volume and unitary holonomy conjugation in
\[
C_{ij}(V)=\delta_{ij}V-h_{ik}Vh^{-1}_{kj}.
\]

For three distinct incident edges define
\[
T_{abc}:=\operatorname{Tr}_{aux}
  [C_a(K)C_b(K)C_c(V)].
\]
The trace is over the two-dimensional auxiliary fundamental index; multiplication order is the displayed operator order.

For the four-valent source vertex, the epsilon set \(\mathcal O_v\) contains the 24 signed ordered triples obtained by choosing one omitted neighbor and permuting the other three. Its sign is alternating under exchange of any two entries. Define
\[
B_v^{raw}:=\sum_{(a,b,c)\in\mathcal O_v}
 \epsilon(a,b,c)T_{abc}.
\]

## Adjoint derivation

Using the two leg-adjoint identities and reversing the product,
\[
T_{abc}^\dagger
=\operatorname{Tr}_{aux}
 [C_c(V)(-C_b(K))(-C_a(K))].
\]
The two minus signs cancel. By cyclicity of the finite-dimensional auxiliary trace,
\[
T_{abc}^\dagger
=\operatorname{Tr}_{aux}
 [C_b(K)C_a(K)C_c(V)]
=T_{bac}.
\]
Now relabel \(a\leftrightarrow b\) in the full epsilon sum. Since the ordered-triple set is closed under this exchange and epsilon is alternating,
\[
(B_v^{raw})^\dagger
=\sum_{(a,b,c)}\epsilon(a,b,c)T_{bac}
=-\sum_{(a,b,c)}\epsilon(a,b,c)T_{abc}
=-B_v^{raw}.
\]

Thus the **complete** epsilon-oriented raw block is anti-Hermitian. The corresponding Hermitian structural operator is
\[
H_{L,v}^{struct}:=iB_v^{raw},
\qquad (H_{L,v}^{struct})^\dagger=H_{L,v}^{struct}.
\]

No \(B+B^\dagger\) completion is appropriate under these conventions: it would annihilate the exact anti-Hermitian block. Nor should one symmetrize individual ordered triples independently; the Hermiticity property is produced by the signed complete epsilon sum.

## Conditions that must be checked in the implementation

1. **Leg adjoints:** test \(C_e(K)^\dagger=-C_e(K)\) and \(C_e(V)^\dagger=C_e(V)\) as full matrix identities, not merely scalar norms or a few diagonal elements.
2. **Common domain/codomain:** input and output basis conventions must be identical for the operator and its adjoint. The W5 adapter maps final scalar covariant states back to Gauss keys; that map must be isometric and complete on the retained subspace.
3. **No hidden projection:** any charge, spin, or scalar-channel pruning must commute with the adjoint, or be proven exact. A cutoff projection can break the identities at its boundary.
4. **Complete epsilon orbit:** all 24 distinct signed triples must be included exactly once. Partial shards cannot be treated as an operator or used for an adjoint test until coherently assembled.
5. **Numerical residual:** after assembly, compute
\[
r_{AH}=\frac{\|B+B^\dagger\|_F}{\max(\|B\|_F,\epsilon_{num})}.
\]
Also test \(iB-(iB)^\dagger\). Report both residuals and the cutoff leakage.

The randomized finite-matrix unit test added alongside this note verifies the algebraic identity and the epsilon sign convention. It is a unit test of the derivation, **not** a substitute for item 1 on the actual Peter-Weyl legs.

## Physical normalization remains separate

The factor \(i\) above is fixed by the adjoint convention and is not a fit parameter. It does not determine the remaining physical coefficient. The classical Thiemann Lorentzian term is proportional to the epsilon contraction of two \(\{A,K\}\) factors and one \(\{A,V\}\) factor, with convention-dependent powers of \(\beta\), \(\kappa\), and \(\hbar\), as well as regulator and trace-normalization factors. These must be transcribed from one declared convention ledger and checked dimensionally. They must not be tuned to force a joint kernel.

In particular, a paper using \(\beta=1\), a different definition of \(K\), or a different normalization of the SU(2) generators cannot be copied coefficient-for-coefficient into this code without conversion.

## Next numerical gate

For each W5 input and node:

1. Assemble all 24 raw epsilon terms coherently at a fixed cutoff.
2. Form the matrix \(B_v^{raw}\) in the same complete retained Gauss basis for rows and columns.
3. Verify the anti-Hermitian residual above.
4. Only then define \(H_{L,v}^{struct}=iB_v^{raw}\).
5. Repeat at increasing cutoffs and compare both matrix entries and singular values.
6. Apply the independently derived physical prefactor only after the convention ledger is audited.

Until these gates pass, the W5 Euclidean certificate remains valid on its own, but there is no Lorentzian or joint-kernel verdict.
