# Infoton physical-pole audit — 2026-11-26

This audit checks the current graviton–infoton programme and its physical pole criteria.

## Verified finite controls

Fresh command:

```powershell
python scripts/bqg_infoton_bh_foundation_gate.py
```

The script returned `passed: true` for its declared finite controls:

- the finite spin-carrier decomposition and extremal-state projectors;
- TT projection, including zero response for the tested scalar and longitudinal-vector controls;
- the prescribed Mathieu/Floquet example, with `mu_fine = 0.0049999414`, leading estimate `0.005`, and zero measured off-resonance growth in the two controls;
- a separate finite unitary scrambling-and-release Page-curve toy.

The Floquet input is imposed: `xi_derived_from_BQG` is `null`, and the script itself labels the result conditional. Its README and programme state that the finite tests do not establish a physical infoton pole or black-hole dynamics.

## Physical-operator audit

`research/infoton_graviton/PROGRAM.md` describes the needed (Q_I^{ij}) as a candidate operator to be built from route variables. It does not provide a completed gauge-invariant operator surviving the physical constraints and history projector. It likewise leaves the physical connected two-point function, positive spectral density, stable pole, mixed graviton–infoton propagator, and the common-action vertex Γ^(3)_(hII) open.

Therefore the representation and TT controls establish a kinematic channel only. A particle claim would require a physical spectral representation with positive weight and an isolated pole (or a precisely defined resonance), after regulator removal.

## Shared-kernel formulation of the requested bridge

The right common object, if it is derived from one physical history, is a retarded inverse kernel on the TT graviton and candidate information channels,

\[
\mathbb K^R(\omega,\mathbf k)=
\begin{pmatrix}
K_T^R & K_{TI}^R\\
K_{IT}^R & K_I^R
\end{pmatrix}.
\]

Its normal-mode condition is

\[
\det\mathbb K^R(\omega,\mathbf k)=0.
\]

This determinant alone is not a particle proof: the poles must have positive physical spectral weight, and both channels and the off-diagonal terms must be obtained from the same constrained (Z[J]) or Γ. The graviton endpoint additionally requires the massless TT pole with two gauge-reduced helicities. A derived nonzero (K_{TI}) would support mode conversion; a derived Γ^(3)_(hII) could instead source parametric production in a specified background. The toy Mathieu result proves neither coupling nor particle existence.

## Current theorem boundary

The repository supports finite representation/projector facts and a conditional Floquet calculation. It does not yet prove an infoton, its coupling ξ, or a common physical graviton–infoton quantum field theory. The missing chain is: physical operator → projected correlator → positive spectral measure → regulator-stable pole → mixed/three-point functions from the same physical effective action.