# BQG physical causal cone: clock, rigging map, retarded response and IR universality

**Audit date:** 2026-11-11  
**Target branch:** `research/lemaitre-bqg-causal-cone-20261110`  
**Audited base HEAD:** `7c4bd56c7d78bdba333246a7910f6161cac7d28a`  
**Purpose:** extend the causal-cone research plan without promoting finite positive controls to physical results.

## 0. Executive status

This note formalizes the missing interfaces between the existing finite relational controls and a physical causal cone of the full binary quantum gravity (BQG) theory.

It does **not** report a new physical Lorentzian production run. It does **not** claim that a physical clock, physical rigging map, retarded response, graviton propagator, or common IR cone has been computed.

The current status must remain separated:

| Layer | Status | Exact scope |
|---|---|---|
| Frozen route rule (q+2=2^q) | proved conditional arithmetic | If the rule is assumed and (q) is a positive integer, its unique solution is (q=2) |
| Dimension-three route-volume limit | proved for frozen recurrence | Not a derivation of the recurrence or a general spectral-dimension theorem |
| Selected 16-cell completion | exact finite kinematic gluing | Does not prove dynamical or global-topology uniqueness |
| Reduced free-fermion cone | proved within reduced model | Not the physical light cone of full BQG |
| C8 combined relational projector | exact finite positive control | Uses externally declared clock and (R=J) |
| Metric-source functional | exact finite positive control | Not physical (Gamma_{m phys}^{(2)}) or a graviton propagator |
| Physical clock / rigging map | open | Must be derived from the actual constrained graph-changing theory |
| Full Lorentzian production operator | open | No completed unified production carrier is claimed here |
| Physical retarded response | open | No (G^R_{m phys}) result is claimed |
| Common IR cone | open | No demonstrated S/V/TT/matter universality |

The exact finite controls remain useful: they establish that the intended relational and source machinery can be made algebraically consistent. They do not supply the missing physical inputs.

## 1. Three causal notions must remain distinct

The audit uses three different structures.

1. **Graph transition support:** which graph sectors are connected by powers of a graph-changing operator.
2. **Relational time:** evolution of physical observables with respect to a clock or ordered boundary-history variable.
3. **Lorentzian propagation cone:** the limiting support or suppression region of physical retarded response in reconstructed spacetime.

A theorem about one layer cannot be silently promoted to another.

For example, if a positive master operator (M) has local graph support, powers (M^n) constrain the minimum number of graph transitions. But (M) is self-adjoint and its support relation is symmetric. The parameter of (e^{-sM}) is a Euclidean semigroup parameter, not automatically physical time. Likewise, the physical projector alone loses nonzero-spectrum support information and cannot reconstruct the unique graph distance.

## 2. Clock derivation gate

### 2.1. Required starting point

Let the actual constrained theory at regulator level (l) be specified by
[
C_l:mathcal D(C_l)subsetmathcal H_{{\rm kin},l}	omathcal H_{{\rm kin},l}.
]
A clock candidate (T_l) must be a degree of freedom or boundary-history variable belonging to this theory. It must not be inserted as an unrelated external C8 register merely because that gives a convenient finite construction.

If a clock momentum (p_T) exists, deparametrization may seek
[
C_l simeq A_l(T,p_T,	ext{geometry}),
igl(p_T+H_{{\rm rel},l}igr),
]
or another explicitly derived equivalent form. This factorization is a target, not an assumption.

### 2.2. Mandatory tests

The gate `CLOCK_DEPARAMETRIZATION_VALIDITY` must report:

- the exact original constraint and operator ordering;
- the domain on which solving for (p_T) is valid;
- the number of branches and the rule selecting a branch;
- whether the branch is preserved by graph-changing transitions;
- self-adjointness or a justified self-adjoint extension of (H_{{\rm rel},l});
- preservation of the physical constraint surface;
- equivalence between evolution from the original constrained system and from the reduced generator;
- possible loss of states, zero modes, or boundary sectors under deparametrization.

A formal square root or local algebraic solution for (p_T) is insufficient if its domain, branch, or physical inner product is not controlled.

### 2.3. Relational observables

For a kinematical observable (O), construct a clock-conditioned observable (O_l(\tau)). In a genuinely deparametrized sector the expected evolution is
[
O_l(\tau)=
e^{\frac{i}{\hbar}H_{{\rm rel},l}(\tau-\tau_0)}
O_l(\tau_0)
e^{-\frac{i}{\hbar}H_{{\rm rel},l}(\tau-\tau_0)}.
]

The following are separate obligations:
[
[C_l,O_l(\tau)]\approx0,
qquad O_l(\tau)^\dagger=O_l(\tau),
]
and finite, representative-independent physical matrix elements. The weak equality must be defined with respect to the actual physical domain; it must not be replaced by an unqualified matrix identity on the entire kinematical space.

A clock must also be operationally adequate: its readings must distinguish the histories needed by the causal experiment, and its resolution/dispersion must be measured. A clock with recurrences or large dispersion may fail to define a monotone macroscopic time even when its finite algebra is exact.

## 3. Rigging-map and physical-inner-product gate

The finite C8 group average
[
P_{\rm rel}=\frac18\sum_{\tau=0}^{7}G^\tau
]
is an exact finite-group projector for the declared positive-control step (G=S_{\rm clock}\otimes J_{\rm geom}). The repository explicitly identifies (J) as a positive control, not the physical gravitational evolution operator.

For a genuine constraint (C_l), a group-averaging expression such as
[
\eta_l(\phi)=\int d\alpha\,
e^{i\alpha C_l/\hbar}\phi
]
is generally distributional when zero lies in continuous spectrum. Its existence and normalization cannot be inferred from the finite C8 calculation.

The gate `RIGGING_MAP_REFINEMENT_COMPATIBILITY` must establish:

1. the test space and distributional dual on which (eta_l) acts;
2. gauge invariance in the appropriate weak sense;
3. positivity and non-degeneracy of the induced physical inner product;
4. treatment of zero modes, continuous spectrum, and gauge volume;
5. regulator/regularization independence, or a quantified controlled dependence;
6. compatibility between levels.

For refinement embeddings (iota_l), define and independently justify a physical embedding (widehatiota_l^{\rm phys}). Then test a diagram of the form
[
\eta_{l+1}\iota_l
\simeq
\widehat\iota_l^{\rm phys}\eta_l.
]
Do not identify (widehatiota_l^{\rm phys}) with the kinematical embedding by definition.

In addition to norm residuals, track physical kernel rank, near-zero eigenvalue counts, spectral gap, and conditioning. If the gap closes, an absolute residual tending to zero is not sufficient; report residual-to-gap ratios and stability of the physical subspace.

## 4. Physical retarded-response gate

### 4.1. Definition

Only after physical observables and the physical inner product are defined should the theory compute
[
G^R_{AB,l}(x,y;\tau)=
\frac{i}{\hbar}\theta(\tau)
\langle[O_{A,l}(x,\tau),O_{B,l}(y,0)]\rangle_{\rm phys,l}.
]

The source and response operators must be gauge invariant (or correctly dressed relational observables). The state, physical inner product, and operator ordering must all be specified.

Do not substitute:
- a symmetric two-point/Hadamard correlator;
- mutual information;
- an unprojected kinematical commutator;
- a heat-kernel parameter of (e^{-sM});
- a constraint resolvent whose spectral parameter has not been derived as a physical frequency.

These may be useful diagnostics, but none is automatically the retarded physical response.

### 4.2. Intervention-based causality

A nonzero correlation is not equivalent to controllable signalling. Alongside the commutator response, define an operational intervention protocol: perturb a local source coupled to (O_B), compare the expectation of (O_A) with and without the intervention, and verify that the linear-response limit agrees with the retarded correlator.

The experiment must distinguish:
- pre-existing symmetric correlations;
- gauge/constraint-induced correlations;
- actual change in a local observable caused by the intervention.

### 4.3. Front extraction

For a fixed regulator and state, record the unnormalised complex response as a function of relational separation and clock time. A threshold estimator may be defined as
[
r_{\rm front}(\tau;\epsilon)=
\sup\{r:|G^R(r,\tau)|>\epsilon R_{\rm ref}\},
]
but (R_{\rm ref}) and the threshold family must be preregistered and independent of the observed front. Also report alternative estimators, such as the maximum spatial derivative and fitted asymptotic tail.

A finite graph need not have exact compact support. The continuum claim is a controlled limiting statement: for spacelike-separated relational events, the physical retarded response must vanish or obey a rigorously quantified suppression bound as the regulator is removed and the physical volume is increased.

### 4.4. Required negative controls

- Disconnect or zero the relevant interaction channel.
- Randomly permute opaque spatial labels while preserving raw amplitudes.
- Compare retarded response with symmetric correlation and mutual information.
- Vary initial physical states and local interventions.
- Vary regulator, graph size, refinement map, and clock resolution.
- Vary the front threshold and fitting window.
- Use a deliberately nonlocal graph transition as a positive control for detecting a widened cone.

Every run must retain raw complex amplitudes, normalization, truncation error, finite-size error, regulator metadata, and source/observable definitions. Per-slice normalization that forces a visually stable front is prohibited.

## 5. Common IR cone and independent metric reconstruction

For each common physical background and common carrier, independently compute retarded responses in scalar, vector, transverse-traceless tensor, and matter sectors. Let the extracted speeds be (v_S,v_V,v_{TT},v_{\rm matter}). Define
[
\delta_{\rm univ}^{(l)}=
\frac{\max_s v_s^{(l)}-\min_s v_s^{(l)}}
{\overline v^{(l)}},
]
and for directional speeds within a sector,
[
\delta_{\rm iso}^{(l,s)}=
\frac{v_{\max}^{(l,s)}-v_{\min}^{(l,s)}}
{\overline v^{(l,s)}}.
]

The IR universality claim requires both quantities to tend to zero with controlled uncertainty and stable scaling across at least the preregistered comparable levels. A single finite-level equality is not universality.

Independently reconstruct:
- (g_{\rm causal}) from physical intervention/retarded-response fronts;
- (g_{\rm spectral}) from the principal symbol of the physical effective action (Gamma_{\rm phys}^{(2)}(\omega,\mathbf k)).

Compare their null structures up to a justified conformal normalization. The spectral metric must not be calibrated using the causal-front data it is supposed to check. Conversely, the causal metric must not use a pre-imposed target metric or cone.

Report Lorentzian signature, coordinate/relabeling covariance, directional isotropy, sector universality, held-out predictions, uncertainty, and dependence on the physical state/background.

## 6. Reduced free-fermion cone: correct scope and next tests

The existing reduced-model result is
[
E_\pm(k)=\pm\sqrt{4m^2+9\cos^2 k},
\quad
\Delta_{\rm sp}=2|m|,
]
[
\xi^{-1}=\operatorname{arsinh}(2|m|/3),
\quad
v_{\rm front}=\sqrt{4m^2+9}-2|m|.
]
The identity
[
v_{\rm front}/3=e^{-1/\xi}
]
is an exact relation within that reduced model. The contour-shift bound gives exponential suppression outside any velocity (v>v_{\rm front}); it is an asymptotic ballistic front, not strict finite-time compact support.

This result is not a physical light speed of full BQG. Static correlation length and dynamic front are linked by the same analytic spectral structure and parameter (m); their agreement is not an independent metric calibration.

Next reduced-model checks, to be implemented separately from full BQG:
1. Add controlled interactions and compare the free and interacting retarded fronts.
2. Introduce spatially varying (m(x)) and quantify reflection, mode conversion, and local front variation.
3. Test finite-size and finite-time scaling, including order of limits.
4. Check whether the front estimator converges independently of threshold and fitting window.
5. Compare independent scalar/vector/TT-like channels only where the reduced model defines them; do not label a toy channel as a gravitational sector without a derivation.
6. Test anisotropy on multidimensional extensions rather than infer isotropy from a one-dimensional chain.

No result for these extensions is claimed by this note.

## 7. Dimension selector and global topology: audit safeguards

The arithmetic statement
[
q+2=2^q
]
has the unique positive-integer solution (q=2). This is conditional arithmetic. The current `scripts/q2_dimension3_fixed_point_gate.py` explicitly rejects every (q\ne2); therefore it checks the frozen (q=2) recurrence and does not independently derive the selector.

Use distinct symbols until a bridge is proved:
[
n_{\rm route}=\text{number of route bits},\qquad
d_{\rm local}=\dim\mathcal H_{\rm local}.
]
Numerical equality of both to two does not establish identity of the concepts.

The selected 16-cell completion has exact kinematic gluing properties, including its dual (Q_4) structure and opposite shared-face flux cancellation. This establishes compatibility of the selected completion, not uniqueness under all global face pairings.

The required unbiased topology experiment must generate face pairings from the frozen local rules without seeding the 16-cell answer. Canonicalize complexes modulo relabelling/automorphisms; test all vertex/edge/face links, orientability, connectedness, homology and PL-manifold conditions; preserve rejected examples and rejection reasons. A finite random sample can find counterexamples to uniqueness but cannot prove uniqueness of an infinite class. A uniqueness claim requires an exhaustive finite classification with a completeness proof, or a separate mathematical theorem.

## 8. New gate registry (specifications, not yet implemented)

| Gate | Input | Required output | Current status |
|---|---|---|---|
| `Q_SELECTOR_FROM_MICRO_RULES` | Frozen microscopic route rules, no hard-coded q | Complete admissible-q derivation and rejection ledger | Not implemented |
| `UNBIASED_GLOBAL_TOPOLOGY_ENUMERATION` | Local face/edge rules only | Canonical gluing classes, link tests, completeness/bias diagnostics | Not implemented |
| `CLOCK_DEPARAMETRIZATION_VALIDITY` | Actual graph-changing constraint and clock candidate | Branch/domain/self-adjointness/equivalence report | Not implemented |
| `RIGGING_MAP_REFINEMENT_COMPATIBILITY` | Physical test spaces, rigging maps and refinement embeddings | Positivity, non-degeneracy, kernel and inter-level residuals | Not implemented |
| `PHYSICAL_RETARDED_RESPONSE` | Physical observables, state, inner product and interventions | Raw retarded kernels, operational checks, errors and front estimates | Not implemented |
| `COMMON_IR_CONE_UNIVERSALITY` | Same-carrier S/V/TT/matter responses at comparable levels | Universal/isotropy residuals and independent metric comparison | Not implemented |

These names are proposed acceptance gates. This document does not imply that code, fixtures, CI workflows, or passing artifacts for them already exist.

## 9. Relation to the existing production plan

This audit does not relax the existing production requirements. The unified Lorentzian E+S+R production carrier remains a prerequisite. All operator sectors must be generated from the same immutable carrier, background, basis and conventions. Preserve raw complex amplitudes and both forward and adjoint refinement actions. Complete the Lorentzian Schur audit, direct microscopic HDA residual, independent metric calibration and comparable physical L0-L3 levels before physical response extraction.

Do not assemble E, S and R from incompatible carriers. Do not treat E-only L1 shards, Euclidean Krylov precursors, Lorentzian V2/V3, or separate sparse-Fourier blocks as one physical package.

## 10. Acceptance rule

The full BQG causal-cone claim may be considered only after all of the following are independently satisfied:

1. A physical clock or boundary-history ordering is derived from the actual theory.
2. A positive, non-degenerate physical inner product/rigging map is constructed and controlled under refinement.
3. A physical relational generator or equivalent ordered-history amplitude is obtained.
4. Physical retarded response agrees with operational interventions and passes negative controls.
5. Fronts converge with regulator, volume, state, clock-resolution and estimator stability.
6. S/V/TT/matter sectors approach one isotropic cone on the same physical background.
7. The causal metric agrees with an independently reconstructed spectral metric.
8. The full Lorentzian HDA and production-carrier requirements remain satisfied.

Until then, the correct statement is: **the reduced causal cone is proved within its reduced model; the full BQG physical causal cone remains open.**

## 11. Provenance and execution disclaimer

This note was prepared from source-level inspection of the target branch at the audited base HEAD stated above and the existing relational-projector, metric-source, reduced-cone, dimension-selector, global-gluing and benchmark artifacts. The source files were inspected; this commit does not claim that their Python gates were executed in this session. No new numerical run, CI pass, production amplitude package, or physical result is asserted here.
