# W5 K-commutator normalization: exact adjoint and hbar bookkeeping

Date: 2026-11-10

## 1. What the current v10 runner actually computes

The sine-Euclidean engine defines

\[
\widehat H_E^{sine}=\frac{T-T^\dagger}{2i},
\qquad
\widehat K_{comm}=[\widehat V,\widehat H_E^{sine}].
\]

Because both \(\widehat V\) and \(\widehat H_E^{sine}\) are Hermitian,

\[
\widehat K_{comm}^\dagger=-\widehat K_{comm}.
\]

The v10 W5 closure runner uses this anti-Hermitian commutator as its K-leg object. Its raw map is therefore a structural map built from two anti-Hermitian commutators and one Hermitian volume leg:

\[
B_{comm}=\sum_{v}\sum_{(a,b,c)\in\mathcal O_v}
\epsilon(a,b,c)\operatorname{Tr}_{aux}
[K_{a,comm}K_{b,comm}V_c].
\]

This is internally consistent for the direct-adjoint test, but it is not yet the map built from two Hermitian K operators.

## 2. Exact conversion to a Hermitian commutator operator

For any anti-Hermitian commutator \(K_{comm}\), the standard quantum-to-Hermitian conversion is

\[
\boxed{
K_H=\frac{1}{i\hbar}K_{comm}=-\frac{i}{\hbar}K_{comm}.
}
\]

Indeed,

\[
K_H^\dagger
=\left(\frac{1}{i\hbar}K_{comm}\right)^\dagger
=-\frac{1}{i\hbar}K_{comm}^\dagger
=\frac{1}{i\hbar}K_{comm}
=K_H.
\]

For the two K legs, the scalar conversion is exact:

\[
\boxed{
K_{H,a}K_{H,b}V_c
=-\frac{1}{\hbar^2}K_{a,comm}K_{b,comm}V_c.
}
\]

Consequently, if the only difference between the two maps is this conversion,

\[
\boxed{B_H=-\hbar^{-2}B_{comm}.}
\]

The factor is real, so the direct adjoint relation is preserved:

\[
B_H^\dagger=-\hbar^{-2}B_{comm}^\dagger.
\]

This sign must not be lost when later combining Euclidean and Lorentzian maps in a joint Gram calculation. In units \(\hbar=1\), the conversion is still a minus sign.

## 3. What this does NOT fix

The conversion above only makes the commutator-derived K-leg Hermitian. It does not prove that the chosen \([V,H_E^{sine}]\) is the correctly normalized canonical extrinsic-curvature operator of the full Ashtekar-Barbero Hamiltonian.

The full identification can carry additional factors involving the Poisson-bracket convention, \(\kappa\), \(\beta\), the precise Thiemann identity, edge orientation, regulator cell, and volume normalization. Those factors must be derived from one frozen classical Hamiltonian convention and then quantized consistently.

In particular, the classical homogeneous kinetic-sector target already registered in LORENTZIAN_BETA_CANCELLATION.md is

\[
H_E^{kin}=-\beta^2 H_{DW},\qquad
H_L^{corr}=(1+\beta^2)H_{DW}.
\]

That identity fixes the relative classical kinetic-sector coefficient. It does not, by itself, determine the absolute coefficient multiplying the finite Peter-Weyl KKV matrix elements.

## 4. Ordering consequence

Let \(B_{comm}\) denote the raw v10 epsilon-oriented map and let \(B_{comm}^\dagger\) be the independently computed reversed-order adjoint. The minimal symmetric completion of this chosen raw ordering is

\[
B_{sym}=\frac12(B_{comm}+B_{comm}^\dagger).
\]

It is Hermitian by construction. The corresponding Hermitian-K version is

\[
B_{H,sym}=-\frac{1}{2\hbar^2}
(B_{comm}+B_{comm}^\dagger).
\]

This is a well-defined symmetric completion, not yet a proof that it is the unique physical Thiemann ordering. Uniqueness requires matching the regulated classical expression, its commutator/Poisson-bracket convention, and the chosen placement of volume and holonomy factors. Other symmetric orderings can differ by genuine quantum terms.

The anti-Hermitian component

\[
B_{comm,anti}=\frac12(B_{comm}-B_{comm}^\dagger)
\]

must also be reported as a diagnostic; it is not to be silently discarded before the adjoint-closed calculation is complete.

## 5. Frozen bookkeeping for the next W5 stages

Every numerical artifact must identify separately:

1. B_comm: raw map using anti-Hermitian commutators;
2. B_comm_dagger: direct reversed-order adjoint;
3. B_H: map after the exact factor -hbar^-2, if that is the only K conversion;
4. B_H_sym: symmetric completion for the declared raw ordering;
5. additional classical/regularization coefficient multiplying the Lorentzian correction;
6. any trace, edge, node-lapse, and regulator normalization.

Do not fold these into one undocumented number. The full physical W5 Gram is blocked until these entries and the common output basis are frozen.

## 6. Scope

The -hbar^-2 factor is an exact algebraic consequence of converting two anti-Hermitian commutators into Hermitian operators. It is not a completed derivation of the full Lorentzian Hamiltonian prefactor, nor a W5 numerical result, nor a continuum statement.
