# BQG graviton–infoton bridge: physicalization programme

**Branch:** \`research/bqg-infoton-black-holes-20261120\`  
**Status:** candidate-theory programme; exact representation-theory input, conditional dynamics, physical coupling open.

## 1. Compatibility with the frozen BQG architecture

This extension does not add a new fundamental particle or bypass the existing physicalization ladder:

\`binary/qubit → SU(2) quantum geometry → 3D → 3+1D → GR/ADM/HDA → relational history → physicalization\`.

The graviton–infoton sector is downstream of the same microscopic route/intertwiner degrees of freedom, the same constraint images, and the same physical history measure. It must eventually be derived from the constrained theory, not appended as a phenomenological oscillator.

The current exact finite input remains

\[
(1/2)^{\otimes4}=2(0)\oplus3(1)\oplus1(2).
\]

The unique \(j=2\) irrep has dimension five. Its extremal states have \(m=+2,-2\), and form a two-dimensional helicity subspace after a propagation direction and massless TT projection are specified. The full \(j=2\) multiplet is **not** itself a massless graviton: the \(m=0,\pm1\) components must be removed by the correct gauge/constraint reduction.

## 2. Required microscopic operator

Let \(\mathcal H_{\rm route}\) be the information/route Hilbert space already present in BQG. A candidate bosonic collective coordinate \(Q_I^{ij}(x)\) must be constructed from its operators and carry a symmetric rank-two representation. Define the TT projection at nonzero momentum:

\[
\Pi_{ij,kl}^{TT}(\hat k)
=\frac12(P_{ik}P_{jl}+P_{il}P_{jk})-\frac12P_{ij}P_{kl},
\quad P_{ij}=\delta_{ij}-\hat k_i\hat k_j.
\]

Then \(Q_{I,ij}^{TT}=\Pi_{ij,kl}^{TT}Q_I^{kl}\). A scalar or vector information operator has zero overlap with this channel by representation theory. This is the first hard negative control.

The operator must be gauge invariant under the microscopic SU(2) constraints and survive the physical-history projector. A kinematic \(j=2\) state is necessary but not sufficient.

## 3. Coupling: no free \`xi\`

The only admissible leading interaction is obtained by expanding the *same BQG microscopic action* around its physical background:

\[
S_{\rm int}^{(2)}
=-\frac12\int d^4x\,\sqrt{-g}\,
h^{TT}_{ij}\,T_I^{ij,TT}.
\]

The coefficient is fixed by the normalized three-point response

\[
\Gamma^{(3)}_{hII}
=\left.
\frac{\delta^3\Gamma_{\rm phys}}
{\delta h^{TT}\delta Q_I\delta Q_I}
\right|_{\rm physical\ background}.
\]

After canonical normalization of both propagators, the dimensionless mixing/coupling \(\xi\) is a derived ratio of this vertex to the physical pole residues. It may not be set to one or fitted to obtain a desired resonance. If the vertex vanishes in the physical limit, the bridge is falsified at this order.

Two channels must remain distinct:

* **Mode conversion:** a bilinear off-diagonal inverse propagator \(G^{-1}_{hI}\neq0\), with an avoided crossing near \(\omega_h=\omega_I\).
* **Parametric production:** a cubic/background-induced modulation of the information-mode quadratic kernel, with pair production near \(\Omega_{\rm GW}=2\omega_I\).

A Mathieu equation tests the second channel only. It cannot establish mode conversion or the existence of an infoton.

## 4. Required physical calculation

The closure chain is:

1. Construct \(Q_I^{ij}\) from the existing microscopic route variables.
2. Project it through SU(2), diffeomorphism and Hamiltonian constraints and the relational-history projector.
3. Compute its connected two-point function and spectral density.
4. Prove positive spectral weight and identify an isolated or controlled resonance pole.
5. Compute the mixed \(h\)-\(Q_I\) two-point function and the \(hQQ\) three-point vertex from one common effective action.
6. Canonically normalize using the same \(Z_T\) and physical TT pole already used by the BQG propagator programme.
7. Run continuum/refinement and regulator-removal tests before interpreting a pole as a particle.

A candidate particle requires a positive-residue pole (or a precisely defined resonance in a continuum), not merely a finite-dimensional spin label.

## 5. Numerical gate specification

The existing \`oldgate.py\` is retained as a finite control, not promoted to a physical prediction. The new gate must contain:

* exact spin decomposition and projector leakage checks;
* helicity-subspace orthogonality;
* TT representation negative controls (scalar/vector source gives zero TT response);
* Floquet convergence under tighter ODE tolerances and longer integration;
* resonant and off-resonant controls;
* an explicit input field \`xi_derived\` that remains null until computed from \(\Gamma_{\rm phys}\).

The test may verify the conditional equation

\[
\ddot q+\omega_I^2[1+\xi h\cos(\Omega t)]q=0,
\quad \mu\simeq|\xi h|\omega_I/4,
\]

but must label every nonzero \(\xi h\) used in a toy test as an externally chosen dimensionless control parameter.

## 6. What counts as closure

**Finite mathematics closed:** four spin-1/2 constituents contain one \(j=2\) irrep; extremal \(m=\pm2\) states are exact and orthogonal.

**Conditional dynamics closed:** the Mathieu/Floquet instability follows for a prescribed periodic modulation.

**Physical bridge open:** gauge-invariant \(Q_I^{TT}\), physical positive spectral density, physical pole, mixed propagator, and derived \(\xi\).

**No-go criteria:** no TT projection; zero physical mixed vertex; negative spectral weight; no regulator-stable pole; or a required coupling that diverges/vanishes in the continuum limit.

## 7. Compatibility warning

The earlier heuristic inference \(P_{\delta g}(k)\sim k^{1.003414}\) from the observer-smoothing exponent is not a vacuum prediction. The explicit reduced Gaussian TT calculation instead gives \(P_{TT}(k)\sim k^{-1}\). The former remains rejected as an identification of the reduced TT vacuum. A different exponent may only be claimed after direct evaluation of the full interacting physical-history two-point function.

## References inside this repository

* \`docs/information_infoton_review/oldbridge.md\`
* \`docs/information_infoton_review/oldgate.py\`
* \`TT_PROPAGATOR_FIRST_PASS.md\`
* \`TT_VACUUM_TWO_POINT_RESULT.md\`
* \`TT_TO_REAL_PHYSICS_OBSERVABLES.md\`

This programme is a continuation of those files; it does not supersede their status labels.
