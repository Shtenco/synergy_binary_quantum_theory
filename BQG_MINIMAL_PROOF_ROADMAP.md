# BQG minimal proof roadmap: from W5 kinematics to falsifiable IR physics

Status: **active proof contract; no claim of complete physical theory**  
Snapshot audited: 2026-11-15  
Audited branch: `research/bqg-w5-production-audit-20261113`  
Snapshot HEAD: `16aeb223c0caca30abf5083deb8fe669fe5b86ce`  
Related PR: [#83](https://github.com/Shtenco/synergy_binary_quantum_theory/pull/83) (Draft; do not merge on this roadmap alone)

## 0. Purpose and rules of evidence

This document is the dependency-ordered proof contract for the shortest defensible route from the W5 SU(2) construction to an experimentally falsifiable low-energy theory.

Keep three evidence classes separate:

- **Mathematical theorem:** follows from stated definitions and assumptions, with a complete proof.
- **Finite computational control:** a reproducible calculation for a declared finite basis, cutoff, graph, seed and tolerance.
- **Physical result:** follows for the candidate BQG model after its actual constraint family, physical measure/history amplitude, relational observables and regulator limits have been constructed.

A literature citation, a green unit test, a finite positive control, a synthetic matrix example, or an exact identity in a truncated model cannot be silently promoted to a physical result.

No gate may be marked `proved` unless its hypotheses match the actual BQG objects. No gate may be marked `physical` or `experimental_test` merely because a finite proxy passes.

## 1. Dependency graph

```text
W5 representation and SU(2) covariance
  -> well-defined operators, domains, adjoints and linearity
  -> anomaly-controlled complete constraint family
  -> positive physical projector / rigging map and physical inner product
  -> derived relational clock and conditional dynamics
  -> causal retarded response from physical connected correlators
  -> continuum / IR TT kernel and common causal cone
  -> frozen Wilson coefficients and one common physical scale
  -> blind held-out observational predictions and falsification
```

A downstream gate cannot compensate for an upstream failure. In particular, an IR fit cannot repair an undefined physical inner product, and a causal-looking finite kernel cannot replace a retarded response derived from the physical history amplitude.

## 2. Gate A — SU(2) representation and operator correctness

### A1. Representation statement

For every declared W5 node and charged sector, specify the complete representation space, including total `(J,M)`, multiplicity/recoupling labels, endpoint representation factors and the gauge action on every open holonomy index.

For a **fully contracted rank-zero operator** `O^{(0)}`, the target selection rule is

$
\langle J'M';\alpha'|O^{(0)}|JM;\alphaangle
=\delta_{J'J}\delta_{M'M}O^J_{\alpha'\alpha}.
$

The operator may mix multiplicity labels. Do not apply this scalar rule to a partially contracted endpoint tensor.

### A2. Required checks

1. Prove the transformation law of the complete sine-ordered `H_E`, not just a label convention or one cached highest-weight output.
2. Independently evaluate all $M=-J,\ldots,J$ components and all five W5 source nodes for the frozen seed.
3. Check $[J_+,H_E]=[J_-,H_E]=[J_z,H_E]=0$ on the complete output, with scale-robust relative and absolute residuals.
4. Check forbidden-sector weight before any cleanup; preserve the raw leakage ledger and use the frozen acceptance contract.
5. Demonstrate cache-on/cache-off equivalence and a fully instrumented no-pruning reference. An incomplete shadow run is not evidence that pruning is harmless.
6. Prove operator linearity and complex homogeneity on a fixed input/output basis. Accumulate all coherent contributions before any sparsification. Any fixed mask must be input-independent and have a certified operator-norm or induced-norm error bound.
7. Specify a common invariant dense domain for unbounded operators; prove symmetry/adjoint relations there and establish essential self-adjointness or a declared self-adjoint extension where required.
8. Check the complete constraint commutator algebra, including the correct structure functions and graph-changing outputs. Vanishing of a projected commutator alone is not closure.

### A3. Current evidence and limits

- A production W5 artifact previously reported 172 output states and a wrong-sector fraction around $4.42\times10^{-21}$ for one finite seed/node/cutoff. It is evidence for that artifact only; the no-pruning contribution was not thereby excluded.
- The direct all-M gate now tests all three magnetic components for two distinct charged recoupling channels at each of the five source nodes, in both production and instrumented-shadow modes. It also checks the finite 6x6 compression of the sine-ordered operator for Hermiticity. This compressed adjointness check is necessary but not a proof of global self-adjointness. Shadow mode bypasses every threshold routed through PRUNE.keep; it is not a complete no-prune proof until the entire transitive import closure is statically inventoried. The cache-on/cache-off workflow is wired to rerun when this direct SU(2) gate changes. Exact-head CI results are still pending.
- Strict compensated complex accumulation has been integrated at selected V10/V11 boundaries. The repository's own remediation plan still identifies lower-level threshold sites and pruning coverage as incomplete.
- Finite HDA controls are not the complete physical anomaly-free constraint algebra.

**Exit condition A:** the actual W5 operator is SU(2)-covariant within a preregistered bound independently of cache/pruning, all operator/domain/adjoint conditions are stated and checked, and the full constraint algebra closes with the correct structure functions on the declared habitat.

## 3. Gate B — physical projector and positive measure

### B1. Exact finite theorem

For a specified finite family (C_A) and a positive-definite Hermitian constraint metric (G),

$
M_G=\sum_{A,B}C_A^\dagger G^{AB}C_B,
\qquad
\langle\psi|M_G|\psi\rangle
=\|G^{1/2}C\psi\|^2\ge0.
$

so

$
\ker M_G=\bigcap_A\ker C_A.
$

If zero is an isolated eigenvalue, $P_0=\mathbf1_{\{0\}}(M_G)$ is an orthogonal projector. This is an exact finite-dimensional theorem. It does not by itself define the continuum physical Hilbert space.

### B2. Theory-specific obligations

1. Freeze the **actual complete regulated BQG constraint family** $C_A^{(\epsilon,J)}$, including its graph-changing action, domains, adjoints and anomaly prescription.
2. Construct $M_G^{(\epsilon,J)}$ and verify positivity and common-kernel equivalence on that actual family, not only on randomly generated synthetic matrices.
3. State the physical boundary states and their normalization, gauge treatment and observable insertions.
4. Define the normalized spectral-window, heat-kernel or rigging-map sequence. If the gap closes or zero enters continuous spectrum, use the appropriate distributional/RAQ construction rather than assuming a finite-rank projector survives.
5. Prove refinement compatibility and convergence of physical matrix elements, independently of the arbitrary positive metric (G) wherever the theorem says the exact kernel must be invariant.
6. Construct the physical inner product
   $
   (\eta(\phi),\eta(\psi))_{phys}=\eta(\psi)[\phi]
   $
   (or a declared convention-equivalent form), prove Hermiticity and positivity, quotient its null space, and establish a nontrivial completed physical Hilbert space.
7. Prove that the physical observables preserve the rigging-map domain and are well-defined on the quotient.

### B3. Current evidence and limits

- `MASTER_CONSTRAINT_PHYSICAL_PROJECTOR.md` proves the positive-master-operator/common-kernel identity and the isolated-zero finite spectral projector.
- The original `master_constraint_physical_projector_gate.py` remains a synthetic operator-algebra self-test. A separate gate now assembles all 60 actual finite K5 Thiemann-type local constraint maps on the 32D fully-active domain and constructs the corresponding master operator and induced positive inner product. Its strict Gauss-projection-loss, spectral, and physical-Gram checks are running in CI; until that exact-head run passes, treat it as an implementation under verification, not a closed gate.
- `PHYSICAL_PROJECTOR_HISTORY` remains explicitly `open_physical` in `physicalization_gates.json`. The repository itself says the continuum/rigging-map limit and physical boundary amplitude remain open.

**Exit condition B:** a nontrivial, positive, regulator-compatible physical Hilbert space has been constructed from the actual BQG constraints, with physical observables and boundary amplitudes defined on it.

## 4. Gate C — relational dynamics and causality

### C1. Derived relational time

The clock must be part of the microscopic theory or a rigorously justified matter/reference sector. Do not add an external clock merely to create a frequency variable.

If a deparametrization is valid, establish a constraint of the form

$
P_T+H_{phys}=0,
\qquad
i\hbar\partial_T\Psi=H_{phys}\Psi,
$

including self-adjointness of `H_phys`, the clock's physical range, and the domain on which conditional states and probabilities are defined.

### C2. Required relational checks

1. Derive the clock observable and its coupling from the BQG graph-changing dynamics.
2. Construct conditional states/probabilities from the physical inner product, not from an externally imposed tensor factor unless that factorization is itself derived.
3. Prove nontrivial conditional evolution while the total state satisfies the constraints.
4. Compare the relational propagator with gauge-invariant boundary/projector amplitudes in their common regime.
5. Test clock changes and admissible clock reparametrizations; quantify any multiple-choice problem rather than hiding it.
6. Define retarded response from the physical connected generating functional. Establish the correct analytic prescription, spectral positivity/unitarity conditions and vanishing of response outside the allowed causal support.
7. Derive the effective causal cone from the same physical IR action for gravitational and matter/photon sectors. A finite Lieb–Robinson bound is a regulator-level control; its continuum scaling and common-cone limit must be shown.

### C3. Current evidence and limits

- `Q2_RELATIONAL_HISTORY_PROJECTOR.md` contains an exact finite C8 Page–Wootters-style positive control. Its clock is external and `R=J` is not the physical graph-changing gravitational evolution operator.
- The finite source/metric generating-functional controls are not a spacetime physical TT kernel.
- `PHYSICAL_PROJECTOR_HISTORY`, `CONNECTED_INTERBLOCK_HISTORY` and `PHYSICAL_TT_KERNEL` remain open physical gates.

**Exit condition C:** a clock derived from the candidate theory yields unitary/consistent conditional dynamics and a retarded physical response with a demonstrated causal support/cone, agreeing with the physical boundary formulation.

## 5. Gate D — continuum and IR limit

### D1. First required IR object

Derive the physical connected metric generating functional from Gate B, then its Legendre effective action and transverse-traceless two-point kernel:

$
K_{TT}(\omega,\mathbf{k})
=\Pi_{TT}\frac{\delta^2\Gamma[g]}{\delta g\,\delta g}\Pi_{TT}.
$

The frequency `omega` must be the physical relational frequency from Gate C. A constraint spectral parameter, finite character angle or Euclidean/source parameter is not interchangeable with it.

### D2. Required continuum properties

1. Jointly remove graph/spin and other regulators along declared refinement sequences; report each truncation, pruning and finite-size error separately.
2. Recover the leading massless Einstein/Fierz–Pauli TT pole with the correct residue/sign and exactly two physical tensor polarizations after constraints.
3. Demonstrate absence of negative-norm propagating modes and test stability/unitarity in the declared regime.
4. Establish long-wavelength isotropy and the Lorentzian dispersion relation, with direction-dependent spread tending to zero and the dynamical exponent approaching the target value.
5. Derive the six on-shell quartic TT Wilson coefficients from the physical kernel without post-hoc reduction or per-channel fitting.
6. Derive Maxwell/photon dynamics and the matter/gravity common cone from the same connected physical history functional.
7. Derive background/scalar cosmology and lensing/dynamics closure from that same physical measure; do not import a phenomenological component to close a mismatch.
8. Set one common physical scale once, or derive it internally. Hold it fixed for every remaining observable.

### D3. Current evidence and limits

The physicalization registry keeps `PHYSICAL_TT_KERNEL`, `IR_SIX_VECTOR`, `COMMON_SCALE_CALIBRATION`, `DYNAMICAL_MAXWELL_KERNEL`, `PHYSICAL_BACKGROUND_COSMOLOGY`, `PHYSICAL_SCALAR_COSMOLOGY` and `LENSING_DYNAMICS_CLOSURE` open. Existing finite cosmology, lensing and interference controls are references/positive controls, not physical BQG derivations.

**Exit condition D:** one common regulator-independent IR effective action reproduces the required low-energy degrees of freedom and fixes all declared dimensionless coefficients and one common scale before external comparison.

## 6. Gate E — experimental test and falsification

Before inspecting held-out values, preregister:

- observable and estimator definitions;
- complete frozen equations, parameters, priors and allowed inputs;
- one common unit/scale map;
- numerical solver, version, tolerances and uncertainty propagation;
- data cuts, calibration and nuisance treatment;
- accept/reject thresholds and explicit failure conditions;
- immutable prediction artifact/hash and held-out dataset identifier.

Use existing multimessenger and gravitational tests as constraints, not as proof that BQG is correct. A successful fit to data used in calibration is not a prediction.

A physical claim is rejected for the tested phase if, for example, the model cannot produce a stable nondegenerate continuum metric, the required two tensor polarizations, a common causal cone, the massless GR pole, regulator-independent coefficients, or a single scale consistent with held-out observables. Do not rescue a failed frozen prediction by post-hoc parameter or sector-specific scale changes.

**Exit condition E:** at least one independently held-out prediction passes the preregistered test, while all upstream mathematical and physical gates remain satisfied.

## 7. Current CI snapshot (2026-11-15)

The branch HEAD at audit time was `16aeb223c0caca30abf5083deb8fe669fe5b86ce`. PR #83 was open, Draft and reported `mergeable=false`.

At the time of the snapshot:
- W5 full shadow/no-prune run `43695810207` was pending on the audited HEAD.
- W5 cache-on/cache-off run `43695804737` was pending on the immediately preceding HEAD.
- A direct all-(M) SU(2) run `43695675137` was queued on an older HEAD.
- Several rapidly superseded workflow runs were cancelled. Cancellation is neither a pass nor a physical failure.
- A source-staging workflow completed successfully; staging success is not a numerical proof.

Always refresh the run status and verify the exact `head_sha` before treating a workflow artifact as evidence for the current source.

The PR base had diverged: compare against `research/physical-history-causality-20261112` reported 251 commits ahead and 38 behind, with merge-base `8f79f75c64d216de387b19bd731eaeb14046860a`. Do not merge or rewrite history as part of this proof programme; reconcile the base deliberately and rerun all gates afterward.

## 8. Immediate execution order

1. Finish and inspect current W5 all-(M), no-pruning and cache-equivalence artifacts at matching source HEADs. Repair the first genuine failure; do not infer results from cancelled jobs.
2. Complete transitive pruning inventory and linearity/adjointness checks; freeze operator-level sparsification bounds.
3. Audit the actual constraint family and HDA closure against the same declared basis/domain.
4. Extend the finite master-projector theorem to the actual BQG regulated constraints; then build the refinement/RAQ positive-inner-product limit.
5. Derive the clock and physical retarded kernel from that physical construction.
6. Only then calculate the physical TT IR kernel, six-vector, common scale and blind observables.

This ordering is intentionally strict: it is the shortest route that does not skip a logical dependency.
