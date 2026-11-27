# BQG discrete null-optical operator: Jacobi/Riccati audit

**Date:** 2026-11-27  
**Status:** finite-dimensional numerical/mathematical gate passes; derivation from the full BQG constraints remains OPEN.  
**Branch:** \`research/null-cone-optical-riccati-gate-20261127\`

## 1. Purpose and scope

This note tests the first concrete bridge between the repository's reduced graph-causal results and the null-cone geometry used in arXiv:2611.34924.

The tested object is the optical deformation matrix on a two-dimensional screen orthogonal to an affinely parametrized null generator. It is **not** yet an operator derived from the BQG master constraint, physical projector, or graph-changing histories.

The distinction is mandatory:

- proved in this note: algebraic consistency and convergence of a finite-dimensional Jacobi/Riccati discretisation;
- not proved: that the BQG graph state space supplies the screen bundle, null generator, optical tidal matrix, or Lorentzian metric.

## 2. Continuum equations and conventions

Let \(k^\mu\) be an affinely parametrized, twist-free null generator and \(e_A^\mu\), \(A=1,2\), a parallel-transported screen basis. Define the optical matrix

\[
B_{AB}=e_A{}^\mu e_B{}^\nu\nabla_\mu k_\nu.
\]

For a hypersurface-orthogonal null cone, \(B_{AB}\) is symmetric and decomposes as

\[
B_{AB}=\sigma_{AB}+\frac12\Theta\,\delta_{AB},
\qquad
\Theta=\operatorname{tr}B,
\qquad
\operatorname{tr}\sigma=0.
\]

The screen-projected optical tidal matrix is

\[
\mathcal T_{AB}
=R_{\mu\alpha\nu\beta}
e_A{}^\mu k^\alpha e_B{}^\nu k^\beta,
\]

up to the sign convention chosen for the Riemann tensor. In a parallel screen frame the Jacobi equation is

\[
J''(\lambda)=-\mathcal T(\lambda)J(\lambda).
\]

The corresponding Riccati equation is

\[
B'=-B^2-\mathcal T.
\]

In four-dimensional vacuum Einstein gravity,

\[
\operatorname{tr}\mathcal T=R_{\mu\nu}k^\mu k^\nu=0.
\]

The trace gives

\[
\Theta'
=-\operatorname{tr}(B^2)
=-\frac12\Theta^2-\operatorname{tr}(\sigma^2),
\]

which is the twist-free vacuum Raychaudhuri equation. A nonzero shear therefore focuses the null congruence even when the null-null Ricci contraction vanishes.

## 3. Discrete Jacobi recurrence

For step size \(h\), use the centered Jacobi recurrence

\[
J_{n+1}-2J_n+J_{n-1}
=-h^2\mathcal T_nJ_n.
\]

For a constant tidal matrix and initial data \(J_0=I\), \(J_1=I+hB_0\), define the forward optical matrix

\[
B_n^{(+)}
=\frac{J_{n+1}-J_n}{h}J_n^{-1}.
\]

Writing \(F_n=J_{n+1}J_n^{-1}=I+hB_n^{(+)}\), the Jacobi recurrence gives

\[
F_{n+1}=2I-h^2\mathcal T-F_n^{-1}.
\]

Consequently,

\[
\boxed{
B_{n+1}^{(+)}
=(I+hB_n^{(+)})^{-1}B_n^{(+)}-h\mathcal T.
}
\]

Expansion for small \(h\) gives

\[
B_{n+1}^{(+)}
=B_n^{(+)}
-h\left[(B_n^{(+)})^2+\mathcal T\right]+O(h^2),
\]

so the forward Riccati variable has first-order global convergence.

If \(B_n\) and \(\mathcal T\) are real symmetric and the chart is nonsingular, the update preserves symmetry: \(B_n\) commutes with \(I+hB_n\), so \((I+hB_n)^{-1}B_n\) is symmetric, and subtracting \(h\mathcal T\) preserves that property.

## 4. Second-order midpoint observable

For constant \(\mathcal T\), initialize the Jacobi field with its second-order Taylor value

\[
J_0=I,\qquad
J_1=I+hB_0-\frac{h^2}{2}\mathcal T.
\]

A midpoint optical estimator is

\[
\boxed{
B_{n+\frac12}^{\rm mid}
=\frac{2}{h}(J_{n+1}-J_n)(J_{n+1}+J_n)^{-1}.
}
\]

The code checks this against a high-resolution RK4 solution at the corresponding midpoint. This estimator converges at second order in the tested constant-\(\mathcal T\) model. It is a numerical improvement, not a new continuum identity.

The midpoint trace is also compared with the finite-difference logarithmic area rate:

\[
\Theta_A
=\frac{\log\det J_{n+1}-\log\det J_n}{h}.
\]

Before a caustic, the tested difference between \(\operatorname{tr}B^{\rm mid}\) and \(\Theta_A\) converges quadratically.

## 5. Caustics and chart failure

The Riccati variable contains \(J^{-1}\). When \(\det J=0\), the optical chart becomes singular. This is a conjugate point/caustic of the chosen congruence; it is not automatically a curvature singularity or a singularity of the spacetime.

Therefore the implementation **must not** silently invert a nearly singular \(J\). The gate explicitly detects the singular Riccati chart. A production implementation must evolve the Jacobi field through the caustic and change optical charts or use a chart-independent formulation. Pseudoinverse substitution is not accepted as a physical continuation rule.

## 6. Regression results

The script \`scripts/bqg_discrete_optical_riccati_gate.py\` was executed with NumPy. It uses a symmetric, trace-free \(2\times2\) tidal matrix, so the test is in the vacuum optical sector.

Observed results:

| Test | Result |
|---|---:|
| Jacobi/Riccati equivalence comparisons | 79 passed |
| Maximum symmetry defect | \(3.336\times10^{-16}\) |
| Vacuum \(\operatorname{tr}\mathcal T=0\) defect | 0 |
| Vacuum Raychaudhuri trace identity defect | \(1.110\times10^{-16}\) |
| First-order forward-Riccati error ratios under \(h\mapsto h/2\) | 2.0078, 2.0038, 2.0019 |
| Second-order midpoint error ratios under \(h\mapsto h/2\) | 4.0779, 4.0380, 4.0187 |
| Midpoint area/trace defect ratios | 3.7400, 3.8645, 3.9308 |
| Zero-tidal exact-solution error | \(2.719\times10^{-16}\) |
| Caustic chart singularity detection | PASS |

The RK4 reference uses 40,000 fixed steps. The tests are regression checks for the stated finite-dimensional model; they are not interval-arithmetic proofs.

## 7. What this adds to the BQG programme

The prior repository results provide two distinct finite-regulator facts:

1. local graph-changing constraints imply a finite-order support bound for graph-sector transitions;
2. a separate reduced free-fermion model has an exact analytic propagation front.

This note adds a third, still separate, result:

3. a consistent discrete optical/Raychaudhuri layer exists at the level of a two-dimensional Jacobi screen, with symmetry preservation, explicit caustic detection, and controlled continuum convergence.

The missing composition is still substantial:

\[
\text{BQG graph histories}
\longrightarrow
\text{derived null generator and screen}
\longrightarrow
\mathcal T_{AB}^{\rm BQG}
\longrightarrow
J_{AB}
\longrightarrow
B_{AB},\Theta,\sigma
\longrightarrow
\text{Raychaudhuri}.
\]

In particular, \(\mathcal T_{AB}^{\rm BQG}\) must be derived from the actual physical graph-changing dynamics, not inserted as an independent matrix. A valid continuum result must also establish refinement stability, diffeomorphism covariance, and a common causal cone for all physical sectors.

## 8. Next proof obligations

- **O1 — Screen construction:** define the two-dimensional screen bundle from physical BQG data, invariant under graph relabelling.
- **O2 — Null direction:** derive a null generator/relational flow from the physical history construction; do not identify master-flow \(T\) with spacetime time.
- **O3 — Optical tidal operator:** derive a self-adjoint/symmetric \(\mathcal T_{AB}^{\rm BQG}\) from the constraint algebra and physical histories.
- **O4 — Jacobi transport:** prove the discrete transport law and its compatibility with graph changes.
- **O5 — Refinement:** show convergence independent of regulator sequence and chart.
- **O6 — Matter universality:** show that matter and gravitational excitations share the same emergent null cone.
- **O7 — Quantum clock:** establish a self-adjoint relational clock with a controlled domain and finite fluctuations.

Until O1–O7 are met, no full-BQG null cone, anti-gravity effect, or quantum-gravity closure is claimed.

## 9. External comparison and caution

The null-cone preprint explicitly uses the expansion, shear, null constraints, and Raychaudhuri equation to determine the area element along each generator (see its Section 3.2, especially equations (3.22)–(3.25)). Its clock is an affine parameter normalized by gravitational boundary data on the initial cap. These are the correct continuum objects to compare against, but their existence in the preprint does not supply the missing BQG derivation.

A separate anonymous, unrefereed note released on 2026-11-26 claims ultraviolet/domain obstructions for particular linearized prescriptions of the null-cone clock. This is not an established refutation and has not been independently reproduced here. It motivates making self-adjointness, domain, and regulator removal explicit proof obligations rather than assuming the clock is a well-defined operator.
