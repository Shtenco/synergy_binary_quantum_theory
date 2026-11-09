# Adjoint audit of the epsilon-oriented K-K-V block

Date: 2026-11-10

## Executive result

The raw 24-term epsilon-oriented K-K-V block is **not proven anti-Hermitian or Hermitian by the available leg identities alone**. A tempting proof by cyclically rotating the auxiliary trace is invalid: the trace is only over auxiliary indices, while its entries are operators on the gravitational Hilbert space and generally do not commute.

The earlier claim that the full epsilon sum is automatically anti-Hermitian and that multiplying it by i is therefore justified is withdrawn.

## Scope and definitions

With the sine convention,

`H_E^sine = (T - T†)/(2i),   K_sine = [V, H_E^sine]`

the gravitational operator K_sine is anti-Hermitian when V and H_E^sine are Hermitian and the commutator is represented on a common domain.

The expected matrix-covariant leg identities are:

`C_e(K)† = -C_e(K),   C_e(V)† = C_e(V).`

For one ordered triple,

`T_abc = Tr_aux[C_a(K) C_b(K) C_c(V)].`

The raw node block is

`B_raw = Σ_(a,b,c in O_v) ε(a,b,c) T_abc,`

where O_v contains the 24 signed ordered triples.

## Exact adjoint identity

The adjoint reverses the product and contributes two minus signs from the two K legs:

`T_abc† = Tr_aux[C_c(V) C_b(K) C_a(K)].`

This is the safe identity. It does **not** imply `T_abc† = T_bac`.

For ordinary complex-number matrices, the trace is cyclic. Here the matrix entries are operators on the gravitational Hilbert space. The partial auxiliary trace

`Tr_aux(XY) = Σ_(i,j) X_ij Y_ji`

is not cyclic when the entries do not commute. Cyclicity is available only for a full trace over auxiliary and gravitational spaces (when trace-class conditions hold), or under additional commutation hypotheses not established here.

Consequently, the alternating epsilon sign under a↔b is insufficient to infer `B_raw† = -B_raw`. The actual adjoint reverses the order of all three legs and leaves the auxiliary trace in that reversed order.

## What a valid completion requires

1. **Adjoint matrix construction.** Assemble the complete raw map in a common orthonormal basis and independently assemble its adjoint using the reversed-order formula above. Check that the computed adjoint agrees with the conjugate transpose of the assembled matrix.
2. **Physical ordering choice.** Identify the precise Hermitian Lorentzian operator ordering from the chosen quantization prescription or primary reference. If that prescription is symmetric ordering of the raw operator, its Hermitian part is `(B_raw + B_raw†)/2`. If it specifies a different ordering, implement that ordering instead. Do not assume that `i B_raw`, `B_raw+B_raw†`, or any other completion is physically correct without this derivation.

A Hermitian part can always be formed algebraically for a finite matrix, but that fact alone does not establish that it is the intended quantization of the classical Lorentzian constraint.

## Required numerical gates

1. Verify the full matrix adjoints of C_e(K) and C_e(V), including auxiliary indices and Hilbert-space matrix elements.
2. Ensure input and output basis conventions are identical and orthonormal. Check the covariant-scalar-to-Gauss map as an isometry on the retained subspace.
3. Prove that all charge/spin/scalar pruning is adjoint-compatible, or include the omitted sectors in the matrix.
4. Assemble all 24 epsilon terms coherently. Partial shards are not operators and cannot be adjoint-tested as complete blocks.
5. Compute the residual between the adjoint of the assembled raw matrix and the independently assembled reversed-order expression:

`r_adj = ||B_assembled† - B_reversed†||_F / max(||B_assembled||_F, eps_num).`

6. Compare increasing cutoffs and report leakage, matrix residuals, and singular-value stability.

## Literature cross-check (not a coefficient transplant)

Yang and Ma, *New Hamiltonian constraint operator for loop quantum gravity*, Phys. Lett. B 751 (2015), DOI: https://doi.org/10.1016/j.physletb.2015.10.062, display the classical Lorentzian term with an explicit (1+β²) factor and derive a regulated operator with convention-specific β, κ, ħ, representation-trace, and inverse-volume factors (their Eqs. (1), (9)–(12)). They separately define a symmetric Hamiltonian by ½(H + H†) (their Eq. (19)). This supports the general requirement to specify an adjoint-based symmetric ordering; it does not determine the coefficient or precise ordering for this repository's KKV implementation. Their regulator and operator definitions must not be copied without a full convention conversion.

## Physical normalization remains open

The factor multiplying the Lorentzian term, including the dependence on β, κ, ħ, trace normalization, and regulator conventions, remains a separate audit. A literature formula must be translated into the exact definitions of K, SU(2) generators, holonomy orientation, and commutator normalization used by this code. No coefficient may be fitted to force a joint kernel.

## Current conclusion

The prior suggestion that the 24-term epsilon sum is automatically anti-Hermitian was too strong and is withdrawn. The accompanying unit test checks the reversed-order adjoint and demonstrates that an auxiliary partial trace of operator-valued matrices is generally not cyclic. The actual W5 Lorentzian matrix, its physical Hermitian ordering, cutoff convergence, and the joint-kernel test remain uncomputed.
