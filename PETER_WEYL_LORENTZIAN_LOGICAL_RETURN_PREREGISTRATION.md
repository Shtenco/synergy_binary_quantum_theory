# Preregistration: full epsilon-oriented Peter-Weyl Lorentzian logical return

Status before first run: **frozen finite-regulator test protocol; physical outcome unknown**.

This gate addresses the next operator question left open by `PETER_WEYL_LORENTZIAN_PARITY_FRONTIER.md`:

\[
\boxed{P_{\rm logical} H_L P_{\rm logical}=0\;?}
\]

The calculation is deliberately performed **before** any physical-history interpretation, source fit, cosmology fit, or dark-sector claim.

## 1. Frozen operator ordering

Use the sine-Hermitian Euclidean ordering already tested by the repository,

\[
H_E^{\rm sine}=\frac{T-T^\dagger}{2i},
\qquad
K=[V,H_E^{\rm sine}],
\]

and the exact matrix-covariant legs

\[
C_e(O)=h_e[h_e^{-1},O].
\]

For one ordered triple,

\[
T_{abc}=\operatorname{Tr}_{\rm aux}
\left[C_a(K)C_b(K)C_c(V)\right],
\]

with state action read right-to-left. No cyclic reordering of the auxiliary partial trace is allowed because the matrix entries are noncommuting geometry operators.

The single-Hamiltonian cutoff is frozen at

\[
\boxed{J_{\max}=7/2}.
\]

## 2. Frozen epsilon assembly

For source node `v=0`, use the four neighbors in the repository `PW.NEIG[0]` ordering. For every omitted neighbor `r`, let `(x,y,z)` be the remaining oriented face and use the previously tested face sign `(-1)^r`. Sum all six permutations of the face with their permutation parity:

\[
\boxed{
H_L^{\rm raw}|\psi\rangle
=\sum_{r=0}^{3}(-1)^r
\sum_{\pi\in S_3}\operatorname{sgn}(\pi)
T_{\pi(xyz)}|\psi\rangle.
}
\]

Thus the production object contains exactly **24 ordered triples**. This convention is frozen from the independent noncommuting-operator assembler identity in `scripts/peter_weyl_lorentzian_triple_algebra_gate.py`.

No overall Lorentzian normalization, factor of `1+beta^2`, lapse normalization, continuum scale, or phenomenological coefficient is inserted in this witness. Such factors cannot decide whether a matrix element is structurally zero.

## 3. Frozen input

The first witness column is the first exact all-`j=1/2` Gauss/logical K5 basis state returned by `PW.basis_full_jhalf()`, which in the canonical ordering is the all-`K=0` state.

The physical question for this run is not fitted:

- if the projected logical norm is numerically nonzero, then one tested column already proves that the finite operator satisfies `P H_L P != 0`;
- if that column is zero, **no global zero claim is allowed**. The remaining 31 logical columns must then be evaluated.

The nonzero detection threshold is frozen at an absolute projected norm of `1e-10`. It is a reporting threshold, not a pass/fail target.

## 4. Hard integrity acceptance

The gate passes as a computation-integrity test if all of the following hold:

1. exactly 24 signed ordered triples are assembled;
2. every ordered triple is computed from the sine-Hermitian `K=[V,H_E^sine]`, never from the older plus-ordering K;
3. the final nonscalar source-`J` weight fraction of every ordered triple is below `1e-8` whenever that term is nonzero;
4. the final nonscalar source-`J` weight fraction of the complete epsilon sum is below `1e-8` whenever the sum is nonzero;
5. `C(V)` complete-basis leakage is below `1e-9`;
6. `C(K)` outer complete-basis leakage is below `1e-9`;
7. `C(K)` internal volume-sector leakage is below `1e-9`;
8. no reached Peter-Weyl spin exceeds `Jmax=7/2`;
9. conversion of final `J=0` covariant states back to the Gauss basis is one-to-one on the projected sector.

The historical primitive fixed-index charge diagnostic is **reported but is not a hard acceptance criterion**, consistently with the already-frozen sine-ordered `C(K)` gate: the physical complete gauge-invariant sum, not an intermediate fixed-index branch, is the accepted charge test.

A completely vanishing epsilon sum is a scientifically allowed outcome and therefore does not fail the integrity gate.

## 5. What a nonzero result would mean

A nonzero logical return establishes only

\[
\boxed{P_{\rm logical}H_L^{\rm raw}P_{\rm logical}\neq0}
\]

for the declared finite habitat, ordering and cutoff.

It would justify the next finite calculation:

1. compute the complete `32 x 32` logical matrix;
2. form and audit the Hermitian Lorentzian completion required by the chosen constraint convention;
3. decompose the logical matrix under the exact `S4`/Pauli dictionaries;
4. combine the declared Euclidean and Lorentzian constraints into the finite master constraint;
5. only then insert relational metric/volume sources in the physical zero sector.

It would **not** by itself be a propagating scalar degree of freedom, dark matter, dark energy, a physical Hamiltonian, or a physical `omega`.

## 6. What a zero result would mean

A zero first column means only that this particular all-`K=0` logical state has no direct return under the full raw epsilon sum at the tested cutoff. It does not imply `P H_L P=0` on all 32 logical states and does not imply the absence of collective `j=1` volume dynamics.

## 7. Connection to the scalar physicalization programme

The legal downstream chain remains

\[
\boxed{
H_L\;\text{amplitudes}
\to \{C_A\}
\to \mathbb M
\to P_{\rm phys}
\to Z[J_{\rm scalar}]
\to W[J_{\rm scalar}]
\to \Gamma^{(2)}_{\rm scalar}
\to \text{constraint/gauge reduction}
\to \{\Phi,\Psi,\mu,\Sigma\}.
}
\]

No constraint spectral parameter is identified with physical frequency.