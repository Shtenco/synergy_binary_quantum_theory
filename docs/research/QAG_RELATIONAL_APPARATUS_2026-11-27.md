# QAG relational apparatus toy model — 2026-11-27

## Purpose

This increment adds a source and a position-readout apparatus with explicit
finite-ring position registers. Their interaction depends only on the relative
coordinate. It tests whether an apparatus can become correlated with a source
without breaking the joint translation symmetry.

It is deliberately a finite-dimensional quantum-information model. It is not
a gravitational Hamiltonian, a detector calibration, or a continuous-space
mechanics model.

## Hilbert space and interaction

The tensor order is source position first and apparatus position second:
\[
\mathcal H=\mathbb C^N_S\otimes\mathbb C^N_A.
\]
For sites \(x_s,x_a\in\mathbb Z_N\), define
\[
r=(x_s-x_a)\bmod N,\qquad
V(r)=\lambda\cos(2\pi q r/N),
\]
and
\[
U=\exp[-i\tau V(\hat r)/\hbar].
\]
The implementation stores this diagonal unitary in the source-major basis.
The joint position probability matrix is the direct two-register readout
distribution. Reduced density matrices and the entanglement entropy quantify
source-apparatus correlations.

## Exact symmetry result

Let \(T_N|x\rangle=|x+1\bmod N\rangle\), and
\[
T_{\rm tot}=T_N^{(S)}\otimes T_N^{(A)}.
\]
Because simultaneous translation leaves \(x_s-x_a\) invariant,
\[
[U,T_{\rm tot}]=0
\]
exactly (up to floating-point roundoff). Consequently
\[
\langle T_{\rm tot}\rangle_{\rm in}
=\langle T_{\rm tot}\rangle_{\rm out}.
\]
The test suite checks unitarity, the commutator norm, normalization, this
expectation-value identity, and nonzero entanglement for a nontrivial
relative-coordinate coupling.

## Important limits

1. A finite cyclic position register has no globally valid self-adjoint
   canonical momentum satisfying \([X,P]=i\hbar I\). The trace of any
   finite-dimensional commutator is zero, whereas
   \(\mathrm{Tr}(i\hbar I)\ne0\).
2. The exact conserved quantity here is the joint cyclic translation
   operator (or its modular momentum-sector label), not an unbounded
   mechanical momentum with units calibrated to a laboratory apparatus.
3. The cosine interaction is a generic relational coupling. No gravitational
   constant, source mass, detector mass, energy budget, or physical length
   scale is calibrated.
4. Correlation and entanglement are not by themselves a complete measurement
   model. A detector's amplification, pointer stability, readout noise,
   back-action, and energy/momentum accounting remain to be modeled.
5. This model does not establish repulsive gravity, antigravity, propulsion,
   or a violation of momentum conservation.

## Next physical closure

A stronger model would require a continuous or controlled large-ring limit,
a detector Hamiltonian and operational readout, explicit source/apparatus
kinetic terms, a calibrated interaction switched in time, and a full
conserved generator with all external supports included. Only after those
checks should the toy model be compared quantitatively with the weak-value
postselection channel. It must remain separate from the BQG physical-projector
and TT-kernel gates until those are independently closed.
