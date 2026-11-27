# From BQG graph histories to null optical geometry

**Date:** 2026-11-27  
**Status:** derivation programme / proof obligations; no graph-to-Riemann map claimed as proved.

## 1. Objective

The discrete Jacobi/Riccati gate establishes that a two-dimensional optical screen admits a stable finite-step description. The next problem is not to tune its parameters, but to derive every input from BQG:

\[
\text{physical graph histories}
\to
\text{relational null direction}
\to
\text{screen geometry}
\to
\text{optical tidal operator}
\to
\text{Jacobi transport}.
\]

The graph Laplacian spectrum alone is insufficient. It does not determine a Lorentzian signature, a physical time normalization, a null direction, or a four-dimensional curvature tensor.

## 2. Required data and their provenance

| Object | Candidate BQG source | Required proof |
|---|---|---|
| Spatial metric/area | Fluxes, spins, intertwiners, cell volumes | Nondegenerate semiclassical reconstruction and refinement stability |
| Connection/parallel transport | Holonomies on graph edges | Correct continuum connection and gauge covariance |
| Physical time direction | Derived relational history generator | Must follow from constraints; master-flow parameter is not physical time by definition |
| Null normal | Reconstructed Lorentzian tetrad plus relational evolution | Null norm and normalization must be derived |
| Screen metric | Two-dimensional subspace orthogonal to a null pair | Positive-definite screen metric and graph-relabeling invariance |
| Spacetime curvature | Four-dimensional connection or a valid Gauss-Codazzi reconstruction | Correct Riemann symmetries and continuum limit |
| Optical tidal tensor | Screen projection of spacetime curvature | Symmetric/self-adjoint operator on physical states |

Kinematical spin-network labels may supply candidate geometric data, but they do not alone solve the Hamiltonian constraint or construct a Lorentzian spacetime.

## 3. Geometric construction to be derived

Let \(k^\mu\) be the future null generator and \(n^\mu\) an auxiliary null normal normalized by

\[
k^\mu k_\mu=0,\qquad
n^\mu n_\mu=0,\qquad
k^\mu n_\mu=-1.
\]

The screen projector is

\[
q^\mu{}_\nu
=\delta^\mu{}_\nu+k^\mu n_\nu+n^\mu k_\nu.
\]

Choose an orthonormal screen dyad \(e_A^\mu\), \(A=1,2\), satisfying

\[
e_A\cdot k=e_A\cdot n=0,\qquad
e_A\cdot e_B=\delta_{AB}.
\]

The optical tidal matrix is

\[
\mathcal T_{AB}
=R_{\mu\alpha\nu\beta}
e_A^\mu k^\alpha e_B^\nu k^\beta.
\]

With the Riemann convention fixed in the companion audit, it enters

\[
J''=-\mathcal T J,\qquad
B'=-B^2-\mathcal T.
\]

For a twist-free null hypersurface, \(B\) is symmetric. In four-dimensional vacuum Einstein gravity,

\[
\operatorname{tr}\mathcal T
=R_{\mu\nu}k^\mu k^\nu=0,
\]

because the cosmological term drops out when contracted with a null vector. The trace-free part of \(\mathcal T\) can nevertheless be nonzero and source shear.

This gives a stringent vacuum test for any proposed BQG curvature operator: its screen trace must vanish on physical vacuum states, while its trace-free components must retain radiative tidal information.

## 4. Candidate discrete construction on changing graphs

A graph history is a sequence of physical sectors and transitions,

\[
(\Gamma_0,\psi_0)\xrightarrow{\mathsf R_0}
(\Gamma_1,\psi_1)\xrightarrow{\mathsf R_1}\cdots.
\]

The transition maps \(\mathsf R_n\) must be derived relationally from the constrained theory. They must not be identified with the positive master heat-kernel parameter \(T\).

For each history slice, the theory must construct a two-dimensional screen fibre \(\mathscr S_n\). If the graph changes, adjacent screen fibres need an isometric transport map

\[
U_{n+1\leftarrow n}:\mathscr S_n\to\mathscr S_{n+1}.
\]

The covariant discrete Jacobi equation is then a candidate:

\[
U_{n\leftarrow n+1}J_{n+1}
-2J_n
+U_{n\leftarrow n-1}J_{n-1}
=-h^2\mathcal T_nJ_n,
\]

where all terms are transported into the same screen fibre before addition. The precise direction and composition of the \(U\) maps must be fixed by the chosen bundle convention; this equation is a design target, not yet a theorem.

Required properties:

1. \(U\) is an isometry and respects the physical inner product.
2. \(\mathcal T_n\) is symmetric on a stated dense domain.
3. The update is covariant under graph relabelling and internal gauge transformations.
4. The screen dimension remains two in the macroscopic limit, even if microscopic graph valence changes.
5. The recurrence preserves the discrete symplectic/Wronskian form.
6. The construction is independent of arbitrary choices of graph refinement up to controlled errors.

## 5. How curvature might be obtained from graph data

A possible route begins with holonomy and flux variables:

\[
h_e\in SU(2),\qquad E(S)\in\mathfrak{su}(2)^*.
\]

Small-loop holonomies can encode spatial curvature, while fluxes and intertwiners can reconstruct a spatial triad and area data. However, spatial curvature alone is not the four-dimensional Riemann tensor required by \(\mathcal T_{AB}\).

A valid spacetime reconstruction must also recover extrinsic-curvature information and the normal/time evolution. In a canonical formulation this requires the Hamiltonian and diffeomorphism constraints, their correct hypersurface-deformation algebra, and a controlled continuum limit. Gauss-Codazzi relations may then combine intrinsic and extrinsic data into the required spacetime curvature.

Therefore the following shortcut is forbidden:

\[
\text{graph Laplacian}\not\equiv
\text{spacetime Riemann tensor}.
\]

Nor is a plaquette holonomy by itself a complete null optical tidal operator.

## 6. Area clock as a relational observable

For a classical Jacobi map before a caustic,

\[
A(\lambda)=A_0\det J(\lambda).
\]

Given a target area \(A_*\), a classical area clock can be defined on a monotone segment by

\[
\lambda_*(A_*)
=\inf\{\lambda\geq0:A(\lambda)=A_*\}.
\]

In BQG, this becomes a quantum-observable problem. One must construct a self-adjoint area/clock operator with a specified domain and prove that its spectral projections are meaningful. Replacing the classical determinant by a formal determinant of noncommuting operator entries is not sufficient.

The clock also depends on boundary normalization of the null tangent. This normalization must be carried through graph changes and gluing; otherwise the affine parameter is not a relationally invariant observable.

## 7. Constraint and projector discipline

The existing BQG no-go result shows that the physical projector alone does not reconstruct the graph-transition distance. Consequently, the derivation must preserve the history generator and its local support structure while imposing the physical constraints.

The target order is

\[
\{C_A\}
\to M
\to \text{derived relational histories}
\to \text{screen transport}
\to \mathcal T_{AB}
\to J_{AB}.
\]

It is not permissible to infer causal order from \(P_{\rm phys}\) alone, nor to identify master-flow time with affine spacetime time.

## 8. Falsifiable gates

| Gate | Pass condition | Failure condition |
|---|---|---|
| G1: Lorentzian reconstruction | A nondegenerate metric with signature \((-+++)\) is reconstructed on physical semiclassical states | Only a positive graph metric or arbitrary time normalization is available |
| G2: Null direction | \(k^2=0\) and affine normalization is relationally fixed | Null vector inserted by hand |
| G3: Screen | Positive two-dimensional screen bundle with gauge-covariant transport | Screen depends on arbitrary vertex labels |
| G4: Curvature | Riemann symmetries and Gauss-Codazzi consistency hold | Only spatial plaquette curvature is available |
| G5: Vacuum trace | \(\operatorname{tr}\mathcal T=0\) on physical vacuum states | Nonzero trace without a physical source |
| G6: Jacobi transport | Discrete Wronskian/symplectic structure is preserved and refinement converges | Drift, non-isometric transport, or regulator-dependent cone |
| G7: Matter universality | Gravity and matter share the same continuum null cone | Species-dependent limiting cones without an accepted physical mechanism |
| G8: Clock operator | Self-adjointness, domain, normalization, and finite relevant moments are established | Formal classical clock only or regulator-dependent divergent operator |
| G9: Entropy | Physical boundary-state counting independently yields the area law | Area law inserted as a normalization or inferred from kinematical dimension alone |

## 9. Anti-gravity is a downstream test, not an input

Only after G1–G8 can one examine whether quantum corrections modify the optical tidal operator:

\[
\mathcal T_{AB}^{\rm eff}
=\mathcal T_{AB}^{\rm cl}+\Delta\mathcal T_{AB}^{\rm Q}.
\]

A positive expansion derivative is not by itself anti-gravity. One must distinguish initial-condition-driven defocusing, ordinary Weyl/Ricci effects, and a genuine invariant change in geodesic deviation or probe motion.

No sign or magnitude for \(\Delta\mathcal T^{\rm Q}\) is assumed here.

## 10. Relation to the null-cone preprint and clock-domain challenge

The Chandrasekaran preprint's Section 3.2 uses the conformal metric of null cuts to determine shear and integrates the focusing equation to recover the area element along a generator. Its affine clock is normalized by gravitational boundary data on the initial cap. This is the continuum target for the map above.

An anonymous, unrefereed note dated 2026-11-26 claims that several natural operator prescriptions for the linearized focusing clock fail in a particular Fock realization because of ultraviolet/domain problems. That claim has not been independently reproduced in this repository and is not treated as established. It nevertheless makes G8 non-negotiable: the clock cannot be accepted merely because its classical expression is finite.

References:
- Chandrasekaran, *An algebraic quantization of gravity from the null cone*, arXiv:2611.34924: https://arxiv.org/abs/2611.34924
- Unrefereed clock-domain challenge: https://zenodo.org/records/24563185
- Existing BQG projector-only causality no-go: https://github.com/Shtenco/synergy_binary_quantum_theory/blob/main/BQG_PROJECTOR_ONLY_CAUSALITY_NOGO.md
- Existing graph-changing causal-order theorem: https://github.com/Shtenco/synergy_binary_quantum_theory/blob/main/BQG_GRAPH_CHANGING_CAUSAL_ORDER_THEOREM.md

## 11. Current conclusion

The discrete optical layer is now numerically controlled. The unresolved physics is the derivation of its inputs from BQG.

The next genuine theorem must not be another fit of a front velocity. It must construct, from the physical graph-history theory, a gauge-covariant two-dimensional screen bundle and a symmetric optical tidal operator whose trace and shear satisfy the vacuum null constraints and whose refinement limit yields the continuum Jacobi/Raychaudhuri system.
