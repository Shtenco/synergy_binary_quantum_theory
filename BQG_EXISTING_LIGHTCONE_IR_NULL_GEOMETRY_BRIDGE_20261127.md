# Existing BQG light cone -> Lorentzian IR null geometry

**Date:** 2026-11-27  
**Status:** Analytic bridge for the already-derived reduced cone; finite regression PASS. Full physical BQG null geometry and gravity/matter universality remain OPEN.

## 1. Correction and scope

The repository already contains two distinct causal-cone results. This note does not rebuild either one:

1. the sharp analytic propagation front of the homogeneous reduced two-band fermion model;
2. the graph-sector support bound for the graph-changing master operator, together with a separate finite relational-history locality theorem.

The task is to identify what part of the first result becomes a Lorentzian null cone in a continuum scaling limit, and how that object must be matched to the graph-history and ADM/HDA sectors.

Three objects must remain distinct:

- a finite-lattice propagation front is a bound on propagation in a regulated band;
- a continuum null cone is determined by the principal symbol of a Lorentzian field equation;
- a graph-sector support cone is combinatorial and does not become a spacetime cone until a physical clock and metric distance are derived.

## 2. Existing reduced dispersion and exact front

The existing reduced dispersion is

\[
E_\pm(k)=\pm\sqrt{4m^2+9\cos^2 k}.
\]

Its exact finite-lattice group/front speed is

\[
v_{\rm front}(m)=v_{\max}(m)
=\sqrt{4m^2+9}-2|m|.
\]

The correlation length and analytic strip obey

\[
\mu_c=\xi^{-1}=\operatorname{arsinh}(2|m|/3),
\qquad
v_{\rm front}(m)/3=e^{-1/\xi}.
\]

These are exact statements in the stated reduced model and momentum convention. They are not yet a full-BQG speed of light.

## 3. The continuum null cone is the critical scaling limit

The band crossing is at \(k_0=\pi/2\). Put

\[
q=k-\frac{\pi}{2}.
\]

Then \(\cos k=-\sin q\), and the exact dispersion becomes

\[
E_\pm(q)=\pm\sqrt{4m^2+9\sin^2 q}.
\]

For small q,

\[
9\sin^2q
=9q^2-3q^4+\frac25q^6-\frac1{35}q^8+O(q^{10}).
\]

Thus the leading continuum mass shell follows directly from the existing band:

\[
E^2=4m^2+9q^2+O(q^4).
\]

To obtain a nontrivial massive continuum limit, scale both the gap parameter and momentum:

\[
m=\varepsilon\mu,\qquad q=\varepsilon p,\qquad \varepsilon\to0.
\]

Then

\[
\frac{E_\pm}{\varepsilon}
\longrightarrow
\pm\sqrt{4\mu^2+9p^2}.
\]

Equivalently, in units \(\hbar=1\), the limiting mass-shell relation is

\[
\boxed{\omega^2-c_{\rm IR}^2p^2=(2\mu)^2,\qquad c_{\rm IR}=3.}
\]

The corresponding local 1+1 Lorentzian line element is

\[
ds^2=-9\,dt^2+dx^2,
\]

and its two null directions satisfy

\[
\frac{dx}{dt}=\pm3.
\]

This is the precise sense in which the existing reduced light cone supplies a Lorentzian IR null cone: it is the principal, linear-in-momentum cone at the critical scaling point.

### Group-velocity check

In the scaling limit,

\[
v_g(p,\mu)
=\frac{\partial\omega}{\partial p}
=\frac{9p}{\sqrt{4\mu^2+9p^2}}.
\]

For finite p and mu, this is subluminal, and

\[
\lim_{|p|/|\mu|\to\infty}|v_g|=3.
\]

The null speed is the characteristic speed of the continuum equation, not the group velocity of every massive wave packet.

## 4. Why the finite-m front is not the null speed

For every fixed nonzero lattice mass parameter,

\[
v_{\rm front}(m)<3.
\]

This is not a contradiction. The finite-m front maximizes the group velocity over the entire regulated Brillouin zone. It therefore includes lattice-scale dispersion, while the Lorentzian null cone is the principal symbol of the joint critical scaling limit.

Under \(m=\varepsilon\mu\),

\[
v_{\rm front}(\varepsilon\mu)
=\sqrt{9+4\varepsilon^2\mu^2}-2\varepsilon|\mu|
\longrightarrow3.
\]

So the two statements are consistent:

\[
\boxed{
\text{finite lattice: }v_{\rm front}(m)<3\ (m\ne0);
\qquad
\text{critical IR cone: }c_{\rm IR}=3.
}
\]

Do not use \(v_{\rm front}(m)\) as a mass-dependent spacetime light speed. It is a regulated band-front velocity. The universal candidate null speed is the critical/continuum principal-symbol coefficient.

## 5. Units and the HDA normalization problem

The value 3 is dimensionless in the reduced convention. If the microscopic hopping scale is J, lattice spacing is a, and hbar is restored, then

\[
c_{\rm IR}^{\rm phys}=\frac{3Ja}{\hbar}.
\]

Separately, the ADM/HDA branch selects

\[
c_{\rm eff}=\frac12,\qquad AB=1,
\]

and its normalized two-derivative TT equation has \(c_T^2=AB=1\) in units fixed by the hypersurface-normal time coordinate.

These values cannot simply be declared equal. The theory must derive one common conversion between the reduced model's time/length units and the ADM lapse/normalization:

\[
\boxed{
\frac{3Ja}{\hbar}
\stackrel{?}{=}c_T^{\rm phys}.
}
\]

The conversion may be fixed once by the physical clock and geometric scale. It may not be tuned independently for different fields or polarizations.

A successful joint gate must show, in the same scaling window:

1. the reduced fermionic cone has \(c_{\rm IR}=3Ja/\hbar\);
2. the Lorentzian ADM/HDA tensor principal symbol has the same physical speed after the single derived unit conversion;
3. matter and gravity share that cone;
4. regulator/refinement corrections vanish at leading derivative order.

The existing binary/HDA uniqueness chain is conditional on the microscopic dynamics actually realizing its isotropy and HDA hypotheses. It is not by itself this common-cone proof.

## 6. What the graph causal cone contributes

The separate graph-changing result states

\[
\Pi_{\Gamma'}M^n\Pi_\Gamma=0
\quad\text{if}\quad d_C(\Gamma,\Gamma')>2n,
\]

hence

\[
d_M(\Gamma,\Gamma')\ge
\left\lceil\frac{d_C(\Gamma,\Gamma')}{2}\right\rceil.
\]

This is a finite support theorem for the master operator. The master parameter T is not physical time, and the physical projector alone loses the transition-order information.

The relational-history positive control has a local step R, but its physical BQG derivation remains open. Once a genuine physical history step is derived, a candidate graph speed bound would require a physical time interval Delta-tau per step and a reconstructed spatial edge length ell-star:

\[
v_{\rm graph}\le\frac{\ell_*}{\Delta\tau}.
\]

Neither ell-star nor Delta-tau may be inserted by convention. Their common scaling must be derived from the same physical state and refinement map. The comparison target is then

\[
v_{\rm graph}^{\rm IR}
=
c_{\rm IR}^{\rm phys}
=
c_T^{\rm phys},
\]

with the understanding that a graph support bound and an amplitude front are different inequalities until their operators and clocks are connected.

## 7. Existing geometric ingredients and missing arrows

The repository already contains relevant partial ingredients:

| Layer | Existing result | What it does not yet supply |
|---|---|---|
| Reduced matter-like propagation | Exact two-band dispersion and sharp finite-lattice front | Full physical matter field and universal continuum cone |
| Graph causal support | Exact finite master support bound | Physical time; master T is not affine time |
| Relational history | Combined clock-system projector positive control | Physical BQG clock and physical system step |
| Spatial geometry | Flux-to-tetrahedral-metric reconstruction; Gauss/frame quotient | A converged Lorentzian spacetime metric on physical histories |
| Lorentzian kinetic signature | Flux pullback of DeWitt form has one negative conformal direction | Physical extrinsic curvature and normal evolution |
| ADM/HDA | Classical local ansatz selects c=1/2, AB=1 | Full microscopic multi-node off-shell HDA closure |
| K5 geometry | Exact finite SU(2) quantum-link and intertwiner/shape channels | Continuum Lorentzian Riemann tensor |

The correct geometric completion remains:

\[
\text{physical flux/shape state}
\to q_{ab},
\quad
\text{physical conjugate momentum}
\to K_{ab},
\]

\[
(q_{ab},K_{ab},N,N^a)
\to g_{\mu\nu}
\to k^\mu
\to q^\mu{}_\nu
\to \mathcal T_{AB}
\to J_{AB}.
\]

Here K_ab must come from the physical conjugate momentum and Hamiltonian evolution, not from a fitted optical matrix. The ADM reconstruction is

\[
ds^2=-N^2dt^2+
q_{ab}(dx^a+N^a dt)(dx^b+N^b dt).
\]

For a reconstructed unit spatial direction s^a, the local null direction is proportional to n^mu+s^mu, with normalization fixed by the relational clock. The screen is the positive two-plane orthogonal to the null pair. The four-dimensional Riemann tensor must then be reconstructed using the full connection or Gauss-Codazzi relations; a spatial plaquette holonomy alone is insufficient.

## 8. Null-optical equations are downstream, not inputs

After the above geometry exists, define

\[
\mathcal T_{AB}
=R_{\mu\alpha\nu\beta}
e_A^\mu k^\alpha e_B^\nu k^\beta.
\]

Then

\[
J''=-\mathcal T J,\qquad
B'=-B^2-\mathcal T,
\]

and for a twist-free vacuum null congruence in four dimensions,

\[
\operatorname{tr}\mathcal T=R_{\mu\nu}k^\mu k^\nu=0,
\]

\[
\Theta'=-\frac12\Theta^2-\operatorname{tr}(\sigma^2).
\]

The discrete Jacobi/Riccati regression already checks the numerical layer. It does not create \(\mathcal T_{AB}\) from the BQG graph. This note therefore adds no arbitrary tidal matrix and does not infer anti-gravity from the reduced cone.

## 9. Falsifiable joint gates

| Gate | Pass condition | Fail condition |
|---|---|---|
| LC1: reduced IR scaling | E/epsilon tends to sqrt(4 mu^2 + 9 p^2) with O(epsilon^2) scaled error | No stable relativistic principal symbol |
| LC2: front distinction | v_front(m) matches the exact lattice maximum and tends to 3 in critical scaling | Treating finite-m front as universal null speed |
| LC3: unit calibration | One derived conversion matches reduced, ADM tensor and matter speeds | Independent per-sector speed fitting |
| LC4: physical clock | Relational clock and Delta-tau derived from physical constraints/boundary data | Identifying master T with physical time |
| LC5: graph geometry | Graph edge lengths and local cone converge under refinement | Bare graph distance called spacetime distance |
| LC6: metric reconstruction | Nondegenerate Lorentzian g_munu on physical semiclassical states | Only positive spatial metric or inserted null vector |
| LC7: screen/curvature | Gauge-covariant two-screen and full Riemann/optical tidal operator | Plaquette curvature or arbitrary T_AB substituted |
| LC8: universality | Matter and both tensor polarizations share the same leading cone | Species/polarization-dependent IR speeds |
| LC9: optical consistency | Jacobi/Raychaudhuri and vacuum trace identities hold on derived data | Optical equations pass only for hand-chosen tidal matrices |

## 10. Reproducibility

Run the command:

    python scripts/bqg_existing_lightcone_ir_null_bridge_gate.py

The gate checks the exact analytic finite-lattice front maximizer, the critical value c_IR=3, the joint mass/momentum scaling to the Lorentzian mass shell, quadratic convergence of scaled energy and group-velocity errors, the null norm in ds^2=-9dt^2+dx^2, and convergence of v_front(epsilon*mu) to the same critical speed.

It does not claim a full-BQG metric, a quantum-gravity proof, or an anti-gravity effect.

## 11. References in the repository

- BQG_SHARP_ANALYTIC_CAUSAL_CONE.md
- BQG_GRAPH_CHANGING_CAUSAL_ORDER_THEOREM.md
- BQG_PROJECTOR_ONLY_CAUSALITY_NOGO.md
- Q2_RELATIONAL_HISTORY_PROJECTOR.md
- FLUX_DEWITT_SIGNATURE_THEOREM.md
- ADM_HDA_PARAMETER_SELECTION.md
- BINARY_HDA_UNIQUENESS_CHAIN.md
- K5_QUANTUM_GEOMETRY_BRIDGE.md
- BQG_NULL_CONE_GRAPH_TO_OPTICAL_GEOMETRY_MAP_20261127.md

External comparison target: Chandrasekaran, An algebraic quantization of gravity from the null cone, arXiv:2611.34924 (submitted 25 November 2026). It motivates a boundary-null-cone construction with shear-sensitive affine clock, edge/frame modes, and causal-diamond gluing. Its claims are those of a new preprint and are not imported as established BQG results.

## 12. Conclusion

The previously derived BQG light cone is retained and sharpened, not replaced.

The key refinement is:

\[
\boxed{
\text{exact reduced lattice front}
\neq
\text{continuum Lorentzian null speed}
}
\]

at finite nonzero mass deformation. Their common critical scaling limit is c_IR=3 in reduced units. The graph causal cone remains a separate combinatorial support theorem. The next physical milestone is a single calibrated bridge from physical graph histories and flux/shape data to ADM metric, then to a null screen and the optical tidal operator.
