# BH-0: physical projector, Lorentzian history and metric gate

**Status:** open physical gate. The finite master-constraint theorem is exact, but the candidate BQG theory has not yet produced a refinement-stable Lorentzian physical-history amplitude or a physical spacetime metric.

This gate precedes BH-1. A zero of a coordinate expansion computed on an assumed Schwarzschild metric is only a benchmark until the metric and its causal structure are obtained from the same BQG physical dynamics.

## 1. Separate the three mathematical objects

The programme must not conflate:

1. a **kinematical constraint operator** and its finite-regulator zero-sector;
2. a **physical inner product / rigging map** obtained after regulator and boundary limits;
3. a **Lorentzian history amplitude** with relational observables and a derived causal interpretation.

At a finite regulator, for a declared family of constraints (C_A^{(\epsilon)}) and a strictly positive constraint-label metric (G^{AB}),

[
\mathbb M_\epsilon
=\sum_{A,B}(C_A^{(\epsilon)})^\dagger
G^{AB}C_B^{(\epsilon)}\geq0,
\qquad
\ker\mathbb M_\epsilon=\bigcap_A\ker C_A^{(\epsilon)}.
]

If zero is isolated, the exact finite projector is

[
P_\epsilon=\mathbf 1_{\{0\}}(\mathbb M_\epsilon).
]

The heat-kernel family (e^{-T\mathbb M_\epsilon}) may approximate this spectral projector as (T\to\infty). **The parameter (T) is not physical Lorentzian time.** Nor does this finite theorem prove that the zero sector survives refinement.

## 2. Physical-limit requirement

Freeze the regulator sequence, constraint ordering, domains, boundary states and normalization prescription. Then demonstrate convergence of physical matrix elements, not merely convergence of kinematical spectra:

[
\mathcal A_{\epsilon}[B_f,B_i;J]
=
\langle B_f|\eta_\epsilon\,
\mathcal T\exp\!\left(i\int J_I\mathcal O_I^{\rm rel}\right)
|B_i\rangle .
]

The limiting rigging/boundary operation (\eta=\lim_{\epsilon\to0}\eta_\epsilon), if it exists, must be shown to be finite on a nontrivial dense test domain and compatible with the constraints. The source insertions must be relational or otherwise gauge invariant.

Do not demand that an oscillatory Lorentzian path-integral integrand itself be a positive probability density. The correct requirements are:
- a positive physical inner product (or positive semidefinite physical Gram form);
- normalized Born probabilities for a specified measurement/decoherence prescription;
- a well-defined, consistently normalized Lorentzian amplitude;
- composition/constraint consistency and regulator stability.

These are distinct conditions. A Euclidean heat kernel is a projector-control device, not a substitute for Lorentzian dynamics.

## 3. Relational time cannot be assigned by relabelling a spectrum

There are two admissible routes, neither completed by the present finite controls.

### A. Derived clock

A clock degree of freedom must be present in the microscopic model and its constraint must yield a deparametrized or controlled relational evolution, schematically

[
C_{\rm total}=P_T+H_{\rm phys}\approx0,
\qquad
i\hbar\partial_T\Psi=H_{\rm phys}\Psi.
]

The clock's range, monotonicity, backreaction and validity domain must be quantified. The finite (q=2) choice (R=J) is a representation-theoretic positive control, not a derived gravitational Hamiltonian.

### B. Boundary/history evolution

Alternatively, a boundary-to-boundary amplitude may define relational transition observables without a preferred global clock. In that case, the composition law, causal orientation and observable time labels must be derived from the same microscopic amplitudes.

Clock-only averaging is not sufficient: the existing C8 audit shows that it can erase nontrivial history characters. A combined clock-plus-system constraint is the correct finite positive control, but the actual graph-changing BQG constraint must replace the toy step.

## 4. Lorentzian metric reconstruction

The metric used in BH-1 must be reconstructed from relational observables of the physical history. At minimum the construction must supply, in one common normalization:

- spatial metric (q_{ab}^{\rm phys});
- lapse (N^{\rm phys}) and shift (N^a_{\rm phys}), or an equivalent causal structure;
- extrinsic curvature (K_{ab}^{\rm phys}) or an equivalent first-history derivative;
- a nondegenerate metric with Lorentzian signature ((-+++)) in the claimed semiclassical regime.

The QFI resistance graph is a spatial information metric only. It cannot supply lapse, shift, causal cones or extrinsic curvature by itself.

The effective constraint algebra must reproduce the hypersurface-deformation relations in the same regime. In particular, with standard convention choices,

[
\{H[N],H[M]\}
=
D\!\left[q^{ab}(N\partial_bM-M\partial_bN)\right],
]

together with the spatial-diffeomorphism brackets. At the quantum level, the corresponding commutators must be anomaly-controlled on the declared physical domain. An arbitrary background lapse or a fitted causal cone is a failure of this gate.

## 5. Metric response and physical TT sector

Construct the source-deformed physical functional from the same rigging/history object:

[
Z_{\rm phys}[J_g]
=\frac{\mathcal A_{\rm phys}[J_g]}
{\mathcal A_{\rm phys}[0]},
\qquad
W_{\rm phys}=-i\hbar\log Z_{\rm phys}.
]

After relational/gauge reduction and a controlled Legendre transform, obtain the retarded metric kernel and its TT projection,

[
K_{TT}^{\rm ret}(\omega,\mathbf k)
=
\Pi_{TT}\,
\Gamma_{\rm phys}^{(2),\rm ret}(\omega,\mathbf k)
\,\Pi_{TT}.
]

A physical frequency is admitted only after the clock or boundary-history construction above. Constraint eigenvalues, C8 characters, Lanczos coordinates and the master-constraint heat parameter are not (\omega).

The low-energy kernel must have a stable, positive-residue massless TT pole with the correct leading causal cone before a Schwarzschild/Kerr correspondence test is meaningful. The existing reduced free-TT calculations are controls, not a substitute for this connected physical-history kernel.

## 6. Falsifiable BH-0 pass conditions

BH-0 passes only if all of the following are demonstrated on a frozen microscopic model:

1. **Constraint closure:** the regulated constraints and their domains are explicit; anomalies are controlled.
2. **Physical projector:** finite zero-sector projectors are computed and their normalized physical matrix elements have a nontrivial, refinement-stable limit.
3. **Physical positivity:** the induced physical inner product is positive semidefinite and nontrivial; Born probabilities are normalized under a declared measurement prescription.
4. **Lorentzian history:** a normalized causal amplitude and its composition/orientation rules are derived, without identifying Euclidean projector time with physical time.
5. **Metric:** relational (q_{ab},N,N^a,K_{ab}) (or equivalent observables) define a nondegenerate Lorentzian metric with a controlled semiclassical regime.
6. **Causal/TT response:** the same history yields a stable massless TT pole and leading Einstein cone after gauge reduction.
7. **Refinement:** the above statements survive at least two independent regulator/refinement paths with quantified errors.

Any failure blocks BH-1 through BH-7 as physical BQG claims. Reduced Schwarzschild, Regge, Planck-flux and random-unitary controls may still be reported as mathematical tests.

## 7. Repository inputs already available

The following main-branch documents provide exact finite pieces and explicit open boundaries:

- `MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md`: positive master constraint, exact finite kernel theorem and isolated-zero projector; continuum rigging limit remains open.
- `Q2_RELATIONAL_HISTORY_PROJECTOR.md`: exact combined C8 clock/system projector; its system step (R=J) is a positive control, not physical evolution.
- `Q2_RELATIONAL_METRIC_SOURCE_GENERATING_FUNCTIONAL.md`: exact finite source functional and metric tangent response; not a spacetime propagator.
- `HAMILTONIAN_CONSTRAINT_TO_EFFECTIVE_ACTION.md`: required projector → physical amplitude → connected (W) → gauge-reduced (\Gamma^{(2)}) chain.
- `TT_TO_REAL_PHYSICS_OBSERVABLES.md`: observable map from a future frozen physical six-component TT Wilson vector; it does not generate that vector.

The next implementation is therefore not another Schwarzschild toy model. It is a common-regulator physical-projector/history experiment that couples the actual graph-changing constraints to relational metric sources and measures convergence of physical inter-block matrix elements. Only after that succeeds should the horizon finder consume a BQG-derived metric.
