# BQG mathematical remediation: literature-to-code plan
Date: 2026-11-13

This document maps peer-reviewed mathematical methods to the nine active BQG blockers. It distinguishes published theorems from proposed adaptations to the repository. A citation is not evidence that its hypotheses hold for the BQG operators.

## 1. W5 charged-sector SU(2): full representation, not a pointwise label assertion

### Literature

1. E. Alesci, I. Mäkinen, J. Yang, **Graphical Calculus of Spin Networks** (2023), arXiv:2304.00268. It develops SU(2) graphical/recoupling calculus and illustrates Hamiltonian-constraint matrix elements in the spin-network basis: https://arxiv.org/abs/2304.00268
2. A. Perez, **The Spin-Foam Approach to Quantum Gravity**, Living Reviews in Relativity 16, 3 (2013), DOI: 10.12942/lrr-2013-3. Spin-network states are built by contracting representation matrices at nodes with intertwiners: https://doi.org/10.12942/lrr-2013-3
3. Standard Wigner–Eckart theorem: a rank-k tensor maps (J) only to (J'\in |J-k|,...,J+k); rank zero implies (J'=J,M'=M). The theorem applies to the representation on which the operator actually transforms, not to an arbitrarily selected local factor.

### New production evidence

Artifact from Actions run `43446220696`, production matrix job, K seed 0, Jmax2=10, source node 0:

- Input charged label: `J2=0,M2=0`.
- Full sine-ordered (H_E) output support: 172.
- Correct-sector squared norm: 6.641081547822338.
- Wrong-sector squared norm: 2.9364927759944046e-20.
- Wrong-sector fraction: approximately 4.42e-21.
- Relative wrong-sector amplitude: approximately 6.65e-11.
- 780 individual sequence rows contain substantial (J=1) contributions; after the complete oriented-spec sum, the wrong-sector norm is at roundoff scale.

The guard currently throws on the mere presence of any wrong-sector key. This is stricter than the downstream acceptance condition `wrong_charge_fraction < 1e-18`. The production artifact therefore points to a numerical residual/acceptance-contract mismatch, not yet to a physical SU(2) anomaly. The separate shadow-no-prune job has not completed, so pruning is not yet excluded as a contributor.

### Repair

- Replace pointwise exact support rejection in `restore_scalar_charge_M` with a per-charged-node norm test, using the same preregistered fraction tolerance as the existing physical gate ((10^{-18}) in squared norm, i.e. (10^{-9}) relative amplitude).
- Preserve and report the unprojected wrong-sector norm and fraction before discarding any numerical residue.
- Reject if the fraction exceeds tolerance; do not relabel wrong-J states as correct-J states.
- Only after the bound passes, remove the residual wrong-sector rows and restore the original M label for the retained J,M block.
- Add unit tests for below-tolerance residual, above-tolerance leakage, wrong-M leakage, empty output, and multiple charged endpoints.
- Compare the cached highest-weight result against a direct all-M computation on a frozen finite W5 seed/node. Verify the covariance identity (U(g)H_E=H_EU(g)) at matrix level and report \\([J_a,H_E]\\) residuals for all three generators.
- Keep the scalar assertion active as a norm-bounded theorem check; do not disable it.

This is a proposed numerical contract correction. It does not prove the representation assignment until the direct all-M commutator test passes.

## 2. Strict complex accumulation and pruning coverage

### Literature and mathematical principle

A linear operator must satisfy additivity and complex homogeneity. State-dependent intermediate thresholding violates linearity because (P_\epsilon(x+y)\ne P_\epsilon(P_\epsilon x+P_\epsilon y)) in general. The existing L1/L2 expected failures already exhibit this.

A robust approach is a fixed matrix-level sparsity mask chosen independently of the input vector, or an explicit operator approximation with a certified norm bound. The strict sparse accumulator is a correct reference for coherent summation, but an unbounded dictionary is not itself a production solution.

### Implementation contract

1. Freeze the input/output basis and operator matrix.
2. Accumulate all complex contributions to each matrix element before sparsification.
3. If sparsification is required, use a basis/operator-level mask independent of the state being applied.
4. Store both the unpruned matrix-element residual and a certified operator-norm or induced-norm error bound.
5. Inventory every direct comparison against an amplitude threshold transitively, including helper functions and cache wrappers. The gate remains blocked while `lower_level_pruning_coverage_complete=false`.
6. Audit cache keys and replayed pruning traces; cached and uncached runs must produce identical outputs and loss ledgers.

## 3. V10/V11 W5 common carrier

The common carrier must be assembled only after V10 has a closed, finite, adjoint-consistent operator on a frozen basis. V11's common carrier is a shared domain/codomain construction, not evidence that the underlying operator is physical.

Required controls: full complex matrix comparison, common basis hash, direct-adjoint versus conjugate-transpose residual, cache on/off equivalence, cutoff convergence, zero omitted spin-wall branches, zero cap hits, and explicit operator ordering.

## 4. Physical projector and inner product

### Literature

- D. Giulini and D. Marolf, **On the generality of refined algebraic quantization**, Class. Quantum Grav. 16 (1999) 2479–2488, DOI: 10.1088/0264-9381/16/7/321.
- D. Giulini and D. Marolf, **A uniqueness theorem for constraint quantization**, Class. Quantum Grav. 16 (1999) 2489–2498, DOI: 10.1088/0264-9381/16/7/322.
- D. Marolf, **Group Averaging and Refined Algebraic Quantization: Where are we now?** (2000): https://escholarship.org/uc/item/1dw733w0
- K. Giesel and T. Thiemann, **Algebraic Quantum Gravity (AQG)** / Master Constraint Programme literature (see the repository's existing master-constraint references; select the exact operator implementation only after constraint algebra audit).

### Construction

For a positive master constraint (\hat M=\sum_I \hat C_I^\dagger\hat C_I\), if zero is an isolated eigenvalue, the physical projector is the spectral projector (P_{phys}=E_{\hat M}(\{0\})). Numerically, heat-kernel filters (e^{-t\hat M}) may approximate this projector only with a controlled spectral gap/limit study. If zero lies in continuous spectrum, the naive strong limit need not yield a nontrivial normalized projector; use a distributional rigging map/group averaging with a specified test space and convergence prescription.

The physical inner product must be induced by the rigging map, be positive semidefinite, and quotient its null space. A kinematical Gram matrix or heat-kernel parameter is not automatically the physical inner product/time.

## 5. Relational clock from graph dynamics

### Literature

- J. Tambornino, **Relational Observables in Gravity: a Review**, SIGMA 8 (2012) 017, DOI: 10.3842/SIGMA.2012.017, arXiv:1109.0740.
- P. A. Höhn, A. R. H. Smith, M. P. E. Lock, **Trinity of relational quantum dynamics**, Phys. Rev. D 104, 066001 (2021), DOI: 10.1103/PhysRevD.104.066001.
- Altaie, Hodgson, Beige, **Time and Quantum Clocks: A Review of Recent Developments**, Front. Phys. 10 (2022), DOI: 10.3389/fphy.2022.897305.

### Construction

Choose a graph observable (T_G) only if its classical/semiclassical flow is monotonic on the target sector and its quantum POVM has controlled resolution/backreaction. Construct complete observables (O_{|T_G=\tau}) or a quantum reduction map. Cross-check the clock-neutral Dirac, Page–Wootters conditional, and reduced relational-Heisenberg descriptions (the Höhn–Smith–Lock trinity) on the same finite model. Include multiple-clock patches where no single global monotonic clock exists.

## 6. Physical retarded kernel

### Literature

- E. R. Livine and D. Oriti, **Implementing causality in the spin foam quantum geometry**, Nucl. Phys. B 663 (2003) 231–279, DOI: 10.1016/S0550-3213(03)00378-X. Positive-lapse restriction defines causal transition amplitudes, while integration over both lapse signs implements the constraint projector.
- The relational-observable literature above supplies the gauge-invariant clock conditioning needed to make response times physical.

### Construction

Define the retarded response only after the physical inner product and relational clock exist. For complete physical observables (A_T(\tau)), (B_T(\tau')), define a response by differentiating the conditional expectation with respect to a physical source switched on at (\tau'). In a regime where linear response applies, compare with
[
G^R_{AB}(\tau,\tau')=-i\,\theta(\tau-\tau')\langle[A_T(\tau),B_T(\tau')]\rangle_{phys}.
]
The step function is ordered by the physical clock, not by an auxiliary master-constraint heat-kernel parameter. Check gauge invariance, positivity of the physical state, source normalization, reciprocity/adjoint relations, and stability under clock changes. A positive-lapse spin-foam amplitude is a causal transition amplitude, but it is not automatically the retarded commutator kernel.

## 7. Universal causal cone across sectors

### Literature

- H. Wilming and A. H. Werner, **Lieb-Robinson bounds imply locality of interactions**, Phys. Rev. B 105, 125101 (2022), DOI: 10.1103/PhysRevB.105.125101. For k-body systems, exponential Lieb–Robinson bounds are equivalent to exponential decay of interactions under the paper's assumptions.
- **Causality and quantum criticality in long-range lattice models** (2019), PNAS, discusses how response to local perturbations distinguishes linear from non-linear causal regions in long-range systems.

### Construction

Derive a graph-distance commutator/response bound for the regulated BQG operator using the actual interaction range and operator norms. Do not infer the cone from symmetric two-point functions, QMI, or a single reduced dispersion alone. Extract sector-specific front velocities from the same physical retarded kernel for TT, scalar, vector, and matter excitations. A universal cone requires convergence to one common velocity after one shared physical calibration, with regulator, graph-volume, clock-resolution, state, and boundary-condition convergence. If interaction tails are long-range, use the appropriate non-exponential bound rather than forcing a Lieb–Robinson linear cone.

## 8. Independent experimental falsification

No current finite BQG computation is an experimental confirmation. Published independent datasets already constrain particular Lorentz-violating propagation models, but those limits are not a direct confirmation or refutation of BQG until its physical retarded kernel predicts the same observables. Useful reference datasets include:

- GW170817/GRB 170817A multimessenger arrival-time comparison; the event was used to constrain the relative speed of gravity and light, with source-emission timing assumptions explicitly controlling the bound: https://doi.org/10.1103/PhysRevLett.119.161101 and https://doi.org/10.1103/PhysRevD.97.083013.
- LHAASO observation of GRB 221009A; the collaboration published stringent photon-dispersion/LIV constraints: https://doi.org/10.1103/PhysRevLett.133.071501.
- Vardanyan et al., **Revisiting tests of Lorentz invariance with gamma-ray bursts: Effects of intrinsic lags**, Phys. Rev. D 108, 123023 (2023), DOI: 10.1103/PhysRevD.108.123023. The intrinsic source lag can change limits by factors of a few: https://doi.org/10.1103/PhysRevD.108.123023.
- Addazi et al., **Quantum gravity phenomenology at the dawn of the multi-messenger era—A review**, Prog. Part. Nucl. Phys. 125 (2022) 103948, DOI: 10.1016/j.ppnp.2022.103948: https://doi.org/10.1016/j.ppnp.2022.103948.

Define falsifiable observables before fitting parameters:

- relative arrival-time dispersion/birefringence for photons from distant transients;
- gravitational-wave versus electromagnetic arrival and polarization constraints;
- energy-dependent propagation/dispersion of neutrinos or other matter probes;
- equivalence of front velocity among tensor, vector, scalar, and matter sectors.

Map the BQG prediction to an effective-field-theory parameterization with explicit assumptions about source emission lags, cosmological propagation, backgrounds, and detector systematics. Use independent datasets and preregistered likelihoods. Existing astrophysical bounds can exclude parameter regions, but cannot establish BQG or replace a derived physical kernel.

## Execution order

1. Fix the SU(2) norm-contract mismatch and complete direct all-M commutator test.
2. Finish strict accumulator integration and transitive pruning inventory.
3. Rerun V10 at cutoffs 10/12/14, cache on/off, production/shadow; require zero unresolved loss.
4. Build V11 common carrier and joint Gram only from a closed V10 matrix.
5. Construct RAQ/rigging map and physical inner product.
6. Construct and validate graph-derived relational clock.
7. Build physical retarded response.
8. Derive common causal bound and compare all sectors.
9. Perform independent astrophysical/experimental falsification.

## Current status

This is a literature-grounded remediation plan. It is not a claim that the remaining physical gates have closed. The production probe is a concrete numerical finding. The cached scalar restore guard has been changed to a per-node wrong-(J,M) squared-weight contract of `1e-18`, matching the existing downstream acceptance criterion, while retaining fail-closed behavior and recording the discarded norm. Seven unit tests were added; the first CI run found a missing `math` import, now fixed, and the rerun is pending. An independent direct all-M commutator gate now evaluates the actual sine-ordered `H_E` on a J=1 multiplet and checks `[J_+,H_E]=[J_-,H_E]=0` in production and instrumented-shadow modes. The no-prune diagnostic and direct all-M gate remain necessary before promoting the SU(2) repair.
