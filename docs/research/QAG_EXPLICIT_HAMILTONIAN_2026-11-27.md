# QAG-1 explicit Hamiltonian checkpoint

**Date:** 2026-11-27  
**Scope:** finite-grid quantum-information toy model; not a gravity simulation.

## New layer

- `qag_hamiltonian.py` constructs the branch Hamiltonian from an explicit linear potential:
  `V_s(X)=-F_s X), hence (mathcal F_s=-dV_s/dx=+F_s).
- It computes the exact diagonal propagator (U_s(t)=exp[-iH_st/hbar]).
- `branch_unitaries_from_linear_potential` creates both controlled branch unitaries from the same declared convention.

For the right branch,
[
U_R(t)=exp!left(+iF_RtX/hbaright)
      =exp!left(iDelta p X/hbaright),
qquad Delta p=F_Rt.
]
This agrees exactly with the existing momentum-kick convention and produces a positive probe momentum shift for positive (F_R).

## Tests added

`tests/test_qag_hamiltonian.py` covers:
1. sign and Hermiticity of the potential Hamiltonian;
2. zero-duration identity;
3. unitarity and exact equivalence to the momentum-kick propagator;
4. momentum shift (F_Rt);
5. rejection of non-diagonal input by the diagonal propagator helper;
6. rejection of negative duration.

The dedicated workflow now compiles and executes this module and its tests as part of the QAG suite.

## Important finite-grid caveat

The periodic grid has a seam. A globally linear potential is not periodic, so this implementation is valid as a local packet-support toy model only when the wavepacket remains well away from the seam and boundary wraparound is negligible. It must not be presented as a globally well-defined periodic Hamiltonian or as a Newtonian/Einstein gravitational field.

## What this closes — and what it does not

Closed at toy-model level: the branch unitary is no longer an unexplained primitive kick; it is derived from a declared potential Hamiltonian with an explicit force sign and duration.

Still open:
- a dynamical source degree of freedom and reciprocal source recoil;
- a Stinespring measurement apparatus and its momentum/energy bookkeeping;
- a non-periodic or boundary-controlled continuum limit;
- a gravity-specific source–probe Hamiltonian derived from a stated Newtonian or relativistic approximation;
- BQG relational position observables and physical projector.

No claim of physical antigravity, net propulsion, free energy, or violation of momentum conservation follows from this toy model.
