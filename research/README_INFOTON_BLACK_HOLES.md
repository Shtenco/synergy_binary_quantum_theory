# BQG infoton and black-hole research extension

This branch consolidates two related but distinct physicalization programmes:

1. [Graviton–infoton bridge](infoton_graviton/PROGRAM.md)
2. [Black-hole sector](black_holes_bqg/PROGRAM.md)

They share the same prerequisite: a physical Lorentzian sector obtained from the existing BQG constraints and relational history. Neither programme may bypass that dependency.

## Reproduce finite controls

From the repository root:

\`\`\`bash
python scripts/bqg_infoton_bh_foundation_gate.py
\`\`\`

The gate checks:
- four-qubit spin decomposition and extremal \(m=\pm2\) states;
- TT projection selection rules and scalar/vector negative controls;
- conditional Mathieu/Floquet resonance and numerical convergence;
- a finite unitary scrambling-and-release Page-curve control with norm conservation and final radiation purity.

A passing result is a finite mathematical/numerical control only. It is **not** a physical BQG infoton, black-hole solution, Hawking derivation, or BQG Page curve.

## Frozen status at branch creation

| Component | Status |
|---|---|
| Four-qubit \(j=2\) carrier | Exact finite representation theory |
| Massless helicity subspace | Two-state kinematic identification; physical reduction still required |
| TT selection rules | Exact kinematic projector control |
| Mathieu/Floquet instability | Conditional toy dynamics; coupling not derived |
| Infoton physical pole | Open |
| Physical Lorentzian metric/history | Open |
| Dynamical horizon | Open |
| BQG area operator and \(1/4\) entropy coefficient | Open |
| Hawking spectrum and greybody factors | Open |
| Evaporation from BQG amplitudes | Open |
| Page curve from BQG radiation density matrix | Open |

## Critical consistency rules

- Do not identify the QFI resistance graph metric with a Lorentzian spacetime metric.
- Do not identify the \(j=2\) carrier with an observed particle.
- Do not set the graviton–infoton coupling \(\xi\) by hand.
- Do not identify a tensor-network cut-entropy upper bound with black-hole entropy.
- Do not insert the Page curve by taking the minimum of two entropy curves.
- Do not promote observer-smoothing exponents to vacuum spectra. The explicit reduced Gaussian TT vacuum result is \(P_{TT}(k)\sim k^{-1}\); the full interacting physical-history result remains open.

A physical completion requires one common frozen BQG microscopic model to pass the constraint, continuum, horizon, entropy, Hawking-flux, energy-balance and radiation-unitarity gates without changing models between stages.
