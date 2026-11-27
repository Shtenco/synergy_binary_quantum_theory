# BQG physical null-cone bridge: HDA, spectral front, graph support, and supplied cone provenance

Date: 2026-11-27

Status: **derived kinematical bridge and fail-closed comparison protocol; no claim that the full physical BQG Lorentz cone has already been derived.**

This note corrects a recurring category error: the repository contains several objects called causal/light cones, but they are not interchangeable. In particular, a spectral propagation front and a graph-sector support bound do not, by themselves, identify the null cone of a reconstructed Lorentzian metric.

## 1. Provenance discipline for the two previously supplied light-cone inputs

The research task states that two light-cone inputs were already supplied. They must remain the primary geometric inputs and must not be replaced by unrelated repository calculations. The historical conversation index available during this audit does not retain their actual image/file payloads, labels, coordinates, or equations. Therefore this note deliberately does **not** guess their contents and does **not** relabel the repository's analytic spectral cone and graph-sector bound as those two supplied objects.

When the original payloads are available in the active research workspace, register them verbatim in a provenance manifest with:
- source filename or immutable URL and content hash;
- which diagram is cone 1 and which is cone 2;
- axes, coordinate conventions, signature, units, orientation, and any marked events;
- whether each boundary is a metric null surface, a discrete causal boundary, a signal front, or a schematic drawing;
- all parameters and any stated calibration.

Until that manifest is populated, the derivations below use only explicitly identified repository equations and state their conclusions conditionally. No user-supplied geometry is silently reconstructed from memory.

## 2. Four distinct notions that must not be conflated

### 2.1 Lorentzian metric null cone

Given a reconstructed spacetime metric of signature (-+++), a tangent vector is null iff
[
g_{mu
u}v^mu v^
u=0.
]
For an ADM split,
[
ds^2=-N^2dt^2+q_{ab}(dx^a+N^a dt)(dx^b+N^b dt),
]
a null ray with coordinate velocity (u^a=dx^a/dt) satisfies
[
q_{ab}(u^a+N^a)(u^b+N^b)=N^2.
]
Equivalently, writing (u^a=-N^a+N s^a), the spatial direction obeys
[
q_{ab}s^as^b=1.
]
This is the geometric light cone. It requires a Lorentzian metric, lapse, shift, and a physically meaningful time/normal normalization.

### 2.2 HDA/TT characteristic cone

For the local two-derivative TT Hamiltonian
[
H_{TT}=int d^3xleft[Api_{TT}^2+rac B4(partial h_{TT})^2ight],
]
Hamilton's equations give
[
ddot h_{TT}=AB,
abla^2h_{TT},qquad c_T^2=AB.
]
The ADM/HDA selection documented in `ADM_HDA_PARAMETER_SELECTION.md` is (c=1/2), (AB=1), hence (c_T=1) in units set by the hypersurface normal. This is the correct leading tensor characteristic target **conditional on** the local ADM ansatz and recovered HDA. It does not yet establish that the microscopic physical theory reaches that phase.

The repository's `BINARY_HDA_UNIQUENESS_CHAIN.md` further makes the implication conditional on microscopic rotational isotropy, standard HDA closure, and a local two-derivative metric phase. Those hypotheses must not be erased.

### 2.3 Reduced free-fermion propagation front

The exact reduced-chain dispersion is
[
E_pm(k)=pmsqrt{4m^2+9cos^2k},
]
with
[
v_{front}(m)=sqrt{4m^2+9}-2|m|.
]
This is an exact result for that homogeneous reduced free-fermion model. It is not automatically the speed of light of full BQG, and it is not a Lorentzian null cone: it is a one-dimensional signal-front bound in the model's chosen coordinate/time units.

A comparison to a metric cone is admissible only after deriving a common conversion of microscopic space and relational time units, and showing that the relevant physical matter and tensor sectors share the same infrared characteristic speed. Choosing (m) to make (v_{front}=1) would be parameter fitting, not a derivation.

### 2.4 Graph-sector support order

The exact result
[
Pi_{Gamma'}M^nPi_Gamma=0quad	ext{if }d_C(Gamma,Gamma')>2n
]
is a support/order theorem for powers of the finite master constraint. It says when a graph-sector block can first become nonzero. It does not assign a metric interval, a spacetime event, or a physical time. The master-flow parameter (T) remains a projector/constraint parameter.

## 3. The legitimate bridge from the BQG variables

The repository's existing geometric ingredients have different scopes:

1. **Flux (	o) spatial metric:** `FLUX_DEWITT_SIGNATURE_THEOREM.md` gives, for a nondegenerate tetrahedral flux triple, (q(E)=|det C|C^{-1}C^{-T}), rank six, with three common-rotation/Gauss null directions. This is a local spatial reconstruction.
2. **Lorentzian kinetic signature:** the same theorem gives the pulled-back DeWitt inertia ((5+,1-,3,0)) in the finite classical flux map. This is a local signature diagnostic, not yet a four-dimensional spacetime metric.
3. **Graph/shape geometry:** `K5_QUANTUM_GEOMETRY_BRIDGE.md` supplies exact finite SU(2) quantum-link geometry, Wilson loops, intertwiners, and tensorial shape observables. It is finite and Euclidean; it is not already a Lorentzian null generator.
4. **Constraint/HDA target:** `ADM_HDA_PARAMETER_SELECTION.md` selects (c=1/2), (AB=1) within its declared local ADM ansatz. The microscopic Lorentzian HDA convergence remains open.
5. **Relational history:** `Q2_RELATIONAL_HISTORY_PROJECTOR.md` is an exact finite positive control with (R=J), explicitly not the physical gravitational evolution operator.
6. **Physical projector:** `MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md` gives an exact finite-regulator zero-sector projector. The continuum rigging map and physical boundary amplitudes remain open.

Thus the current defensible chain is
[
	ext{quantum-link flux/shape data}
	o (q_{ab},K_{ab}) 	ext{candidate coarse data}
	o 	ext{Lorentzian constraint/HDA limit}
	o g_{mu
u}
	o {v:g_{mu
u}v^mu v^
u=0}.
]
The arrow producing (K_{ab}), a common relational normal/time, and the Lorentzian 4-metric is not yet closed by the finite flux theorem alone.

## 4. What is required to turn the supplied cones into an actual BQG null-cone test

For each supplied cone, extract its boundary as a set of events/coordinates or as an explicitly defined discrete relation. Then perform the following checks without refitting after seeing the answer:

**CONE-1 / CONE-2 data gate**
- Determine whether the object is metric, discrete-order, or signal-front data.
- Preserve the original axis labels and coordinate normalization.
- If it is a metric cone, reconstruct or read off (g_{mu
u}), verify signature ((-+++)), and test (g_{mu
u}v^mu v^
u=0) on boundary tangents.
- If it is a graph cone, define the event map and test its order relation; do not infer a metric solely from reachability.
- If it is a signal front, determine the propagating operator, dispersion, and front definition; distinguish group velocity from support bounds.

**BQG prediction gate**
- Construct (q_{ab}) from physical coarse flux states, not from a target cone.
- Derive extrinsic curvature (K_{ab}) from the actual normal evolution/constraint generator.
- Recover the hypersurface-deformation algebra with the same (q^{ab}) as its structure function.
- Reconstruct lapse and shift from the same relational boundary/clock normalization.
- Derive the spacetime metric and predict its null cone before comparing with either supplied cone.
- Independently calculate the matter and TT tensor characteristic sets. Their common infrared cone must agree with the metric null cone.
- Repeat under graph refinement and physical-state projection; demonstrate convergence of the cone, not merely of a fitted speed.

## 5. Optical geometry once the metric cone exists

For a null generator (k^mu) and an orthonormal screen basis (e_A^mu), (A=1,2), define
[
B_{AB}=e_A^mu e_B^
u
abla_mu k_
u
=sigma_{AB}+	frac12Thetadelta_{AB}.
]
The optical tidal matrix is
[
mathcal T_{AB}=R_{mualpha
ueta}e_A^mu k^alpha e_B^
u k^eta,
]
and the Jacobi equation is
[
J''=-mathcal T J.
]
The Riccati equation (B'=-B^2-mathcal T) follows where (J) is invertible. In vacuum, (mathrm{tr},mathcal T=R_{mu
u}k^mu k^
u=0), while the trace-free Weyl part can still generate shear and focusing.

This optical layer is downstream of the physical metric cone. The existing discrete Jacobi/Riccati gate validates a numerical integrator for a supplied (mathcal T); it does not derive (mathcal T) from the BQG graph/flux operators. At caustics (det J=0), the Riccati chart fails; evolve (J) through the caustic and change chart rather than silently applying a pseudoinverse.

## 6. Falsifiable joint comparison

In one preregistered regulator/refinement window, require all of the following:

1. Flux reconstruction is nondegenerate and converges to a positive spatial (q_{ab}).
2. The Lorentzian kinetic form has one conformal negative direction before constraints and positive physical TT modes after reduction.
3. Quantum HDA anomaly tends to zero and the same (q^{ab}) appears in the structure function.
4. (c_{eff}	o1/2), (A_{eff}B_{eff}	o1), and the TT cone tends to (c_T=1) in the same relational normalization.
5. The matter front tends to the same (c=1), rather than merely being rescaled independently.
6. Both supplied cones are compared against the parameter-free BQG prediction using their original calibration.
7. Cone disagreement, polarization-dependent leading cones, or a refinement-dependent limiting speed is a failure, not a prompt to retune (m), lapse, or graph weights.

A quartic lattice anisotropy at order (a_*^2k^4) is allowed to survive as an irrelevant correction if the leading (k^2) cone is isotropic and the correction vanishes relative to it in the continuum limit. The complete six-dimensional TT quartic Wilson sector must be retained; do not force a two-parameter truncation before extraction.

## 7. Immediate conclusions

- The repository has an exact reduced-model front, an exact graph-sector support theorem, a finite flux-to-spatial-metric map, and a conditional HDA-to-TT-cone result.
- These are complementary pieces, not two proofs of one full BQG light cone.
- The strongest next derivation is not another fit to (v_{front}(m)). It is to derive the Lorentzian normal evolution and the common physical characteristic cone from the existing flux/shape and constraint operators.
- The supplied light-cone inputs remain reserved for a direct, provenance-preserving comparison; no missing input is fabricated in this note.

## References inside this repository

- `FLUX_DEWITT_SIGNATURE_THEOREM.md`
- `K5_QUANTUM_GEOMETRY_BRIDGE.md`
- `ADM_HDA_PARAMETER_SELECTION.md`
- `BINARY_HDA_UNIQUENESS_CHAIN.md`
- `BQG_CAUSAL_CONE_AND_W5_SPECTRAL_CERTIFICATE.md`
- `BQG_SHARP_ANALYTIC_CAUSAL_CONE.md`
- `BQG_GRAPH_CHANGING_CAUSAL_ORDER_THEOREM.md`
- `BQG_PHYSICAL_GRAPH_TRANSITION_KERNEL.md`
- `Q2_RELATIONAL_HISTORY_PROJECTOR.md`
- `Q2_RELATIONAL_METRIC_SOURCE_GENERATING_FUNCTIONAL.md`
- `MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md`
- `BQG_DISCRETE_NULL_OPTICAL_RAYCHAUDHURI_AUDIT_20261127.md`
- `BQG_NULL_CONE_GRAPH_TO_OPTICAL_GEOMETRY_MAP_20261127.md`
