# BQG five-year physics, mathematics and quantum-theory literature audit (2021-11-14 to 2026-11-14)

**Purpose:** identify recent results that can shorten the route to closing the repository's open BQG gates. This is a literature-to-implementation audit, not a claim that published results prove this repository's operators correct.

**Evidence rule:** a published theorem is reusable only after its hypotheses are mapped to the actual BQG Hilbert spaces, operators, graph family, domains, cutoffs and limits. Numerical tests validate implementations; they do not replace those hypotheses.

## Executive findings

1. **A 2026 spin-foam continuum-limit theorem gives a precise warning and a constructive route:** overly strong convergence can force a topological continuum theory; a weaker distributional limit can instead yield a cylinder rigging map and physical Hilbert space. This directly informs the physical projector/inner-product work.
2. **A 2026 causal Lorentzian spin-foam vertex introduces Toller matrices and a Feynman i-epsilon prescription.** It is a promising causal transition-amplitude construction, but it is not automatically a retarded commutator response.
3. **A 2025 numerical construction implements graph-changing Thiemann Hamiltonian dynamics on 3- and 4-valent spin networks.** This provides a concrete benchmark against which any graph-preserving BQG approximation should be tested.
4. **2025-2026 LQG coarse-graining and covariant-LQG UV fixed-point work provide candidate RG architectures.** They are not drop-in proof of this project's continuum limit.
5. **Recent quantum-clock work shows finite-resource clocks can produce nonlocal-in-time and nonunitary conditional dynamics in one measurement prescription.** The clock gate must include clock quality, measurement prescription and backreaction.
6. **Recent Lieb-Robinson results sharpen locality tests, including long-range interactions and volume-filling operator tails.** The project must derive its cone from measured BQG couplings, not assume a linear cone.
7. **GRB/LHAASO analyses and hierarchical lag models provide a better falsification pipeline**, but intrinsic source lags and model assumptions dominate some current limits.

## 1. SU(2) charged sector and operator covariance

For a genuinely gauge-scalar operator at the relevant node,
\[
[J_a,H_E]=0,\qquad a\in\{x,y,z\},
\]
or equivalently
\[
U(g)H_EU(g)^{-1}=H_E\quad\forall g\in SU(2).
\]
In a coupled basis |J,M;alpha>, with multiplicity/recoupling label alpha, Schur/Wigner-Eckart gives
\[
\langle J'M';\alpha'|H_E|JM;\alpha\rangle
=\delta_{J'J}\delta_{M'M}\,h^J_{\alpha'\alpha}.
\]
The operator may mix multiplicity labels; it cannot mix total (J,M) if the **complete contracted operator** is rank zero under that node's gauge action.

For a rank-k tensor operator, allowed total-spin channels instead obey |J-k| <= J' <= J+k, with magnetic selection M'=M+q. Do not apply scalar selection rules to an open holonomy index before all endpoint indices and intertwiners have been contracted.

### Recent implementation resources

- Alesci, Mäkinen, Yang, **Graphical calculus of spin networks** (2023), arXiv:2304.00268: practical SU(2) recoupling and Hamiltonian matrix elements: https://arxiv.org/abs/2304.00268
- East, Martin-Dussaud, Van de Wetering, **Spin-networks in the ZX-calculus** (2021): diagrammatic proof/circuit representation of SU(2) identities and 3j/4j/6j symbols: https://arxiv.org/abs/2111.03114
- Yang, Zhang, Ma, **Relating spin-foam to canonical loop quantum gravity by graphical calculus** (2021; published just outside the strict five-year window): https://doi.org/10.1103/PhysRevD.104.044025

### Fastest valid BQG route

1. Freeze one small charged W5 basis and explicitly build representation matrices D_J(g) on the **whole local state space**, including external/open holonomy indices.
2. Build the full production matrix H_E without state-dependent pruning.
3. Compute all three residuals R_a=J_a H_E-H_E J_a, plus finite group checks D(g)H_ED(g)†-H_E for several noncommuting group elements.
4. Report absolute and scale-normalized Frobenius/operator-norm residuals, support and charged-sector leakage for each node and each M.
5. Only then expand to all seeds, all K5 nodes and larger cutoffs.

A ladder-only test on one J=1 seed is a useful regression, not a complete covariance proof. Keep the current direct all-M gate as a diagnostic subgate; add full matrix commutators and a representation-action test.

## 2. Strict linearity, sparse sums and pruning

The core algebraic requirements are
\[
A(x+y)=Ax+Ay,\qquad A(\alpha x)=\alpha Ax
\]
for all complex alpha. A state-dependent threshold map P_epsilon generally violates
\[
P_\epsilon(x+y)=P_\epsilon(x)+P_\epsilon(y).
\]
Summing paths/nodes after individually deleting small amplitudes is not an innocuous optimization: destructive interference may remove large intermediate terms and leave a small exact result, or constructive interference may make many individually small paths significant.

### Correct production pattern

- Freeze domain and codomain bases.
- Assemble each matrix element from **all** coherent contributions using stable complex summation (real and imaginary math.fsum, or an equivalent reproducible pairwise accumulator).
- Apply only a **fixed operator/basis mask** independent of the input vector, or publish a controlled approximation with a certified bound.
- For discarded matrix E=A-A_kept, report a bound relevant to the actual norm: e.g. ||E||_2, or a certified upper bound from induced norms. Entrywise tolerance alone is not an operator error bound.
- Test linearity with adversarial cancellation, complex scaling, permutations of summation order, and cache on/off.
- Inventory all threshold sites and all caches, but distinguish amplitude pruning from spectral-domain checks, acceptance thresholds and exact representation selection rules.

### Priority

The V10 node-level SparseLinearAccumulator integration is useful but does not close lower-level pruning in K/V builders or all V09/V10/V11 merges. Build a small **dense reference matrix** for a frozen basis and compare the production sparse operator against it. This provides a golden oracle for every optimization.

## 3. Graph-changing Hamiltonian dynamics: high-priority missing-physics check

Guedes, Mena Marugán, Müller and Vidotto, **Taming Thiemann's Hamiltonian constraint in canonical loop quantum gravity: Reversibility, eigenstates, and graph-change analysis**, Phys. Rev. D 112, 026024 (2025), DOI 10.1103/hjdk-kdhk:
https://doi.org/10.1103/hjdk-kdhk

They give a complete derivation of the graph-changing Hamiltonian action for 4-valent spin networks (and update the 3-valent case), including the volume operator, and implement it numerically. Their calculations report geometric-observable behavior differing from graph-preserving approximations and identify candidate constraint solutions.

### Implication for BQG

If the current BQG production operator is graph-preserving or truncates graph-changing branches, agreement across spin cutoffs is not enough. Add a graph-change benchmark:

- represent graph-changing states canonically (stable graph hash + edge/node labels + intertwiners);
- implement one small 3-/4-valent reference case;
- compare BQG's graph-preserving projection with the full graph-changing output;
- measure volume expectation, constraint residual, norm in omitted graph sectors and reversibility;
- only then classify the fixed-graph approximation as controlled or physically inadequate.

This is not a demand to replace BQG wholesale. It is a required systematic-error measurement.

## 4. Continuum limit, renormalization and coarse graining

### 4.1 Distributional continuum limit and rigging map (major 2026 result)

Bruno, Colafranceschi, Mele, Rovelli, **Structure of the continuum limit of spin foams**, Phys. Rev. D 114, 066005 (2026), DOI 10.1103/7493-9nb7:
https://doi.org/10.1103/7493-9nb7
Preprint: https://arxiv.org/abs/2603.16999

The paper develops an axiomatic framework assigning Hilbert spaces/amplitudes to triangulated manifolds. Under its stated natural assumptions, sufficiently strong notions of convergence force a topological theory (a no-go result). A weaker **distributional** continuum limit can make the cylinder amplitude into a rigging map and define a physical Hilbert space.

**Direct consequence:** do not demand naive pointwise/strong convergence of every spin-foam amplitude if the target is a non-topological gravitational theory. Specify the test space, topology of distributions, cylinder consistency, positivity and null quotient. Check the theorem's exact hypotheses against the BQG amplitude system before claiming it applies.

Suggested rigging-map contract:
\[
\eta:\Phi\to\Phi^*,\qquad
\langle\eta(\phi),\psi\rangle_{phys}:=\eta(\phi)[\psi],
\]
with positivity, reality/symmetry, constraint annihilation in the dual sense, and quotient by the null space. For a group-averaging realization (where the group/action and measure make sense),
\[
\eta(\phi)[\psi]=\int_G dg\,\langle\phi|U(g)|\psi\rangle_{kin}.
\]
The integral can be distributional and need not define a bounded operator on the kinematical Hilbert space.

**Highest-value implementation:** add a cylinder/rigging-map gate on a finite test space, then a distributional refinement sequence. Do not identify the finite master-constraint projector with the continuum rigging map.

### 4.2 Candidate UV fixed point in covariant LQG

Han, **Ultraviolet fixed point in covariant loop quantum gravity**, Phys. Rev. D 114, 044040 (2026), DOI 10.1103/d8s7-jqfl:
https://doi.org/10.1103/d8s7-jqfl

The work studies the complete amplitude summed over two-complexes, introduces spin-network/spinfoam stacks, and identifies a candidate UV fixed point with small-spin condensation. At leading order the fixed point is topological; triangulation ambiguities are reduced to a finite set of boundary coefficients.

**Use:** a candidate ansatz for organizing the sum over complexes and reducing infinitely many triangulation ambiguities to finitely many boundary couplings. Do not treat the candidate fixed point as an established universal fixed point of all LQG models or proof that BQG is at it.

### 4.3 Canonical LQG coarse-graining and RG flow

Assanioussi and Zeiß, **Coarse-grained models for loop quantum gravity and their renormalization** (2026 preprint), arXiv:2608.31152:
https://arxiv.org/abs/2608.31152

It constructs coarse spin-network Hilbert spaces, effective Hamiltonians and explicit nonperturbative renormalization flow equations.

**Use:** adapt the coarse-map/effective-operator architecture to BQG:
\[
\mathcal H_{fine}\xrightarrow{C}\mathcal H_{coarse},\qquad
H_{coarse}\approx C H_{fine} C^\dagger
\]
(with the actual isometry/embedding and inner-product conventions made explicit; this schematic formula is not a universal RG identity). Track generated couplings and test whether the operator family closes under coarse graining.

### 4.4 Spin-foam refinement and tensor networks

Asante, Dittrich, Steinhaus, **Spin Foams, Refinement Limit, and Renormalization** (2022/2023 chapter), arXiv:2211.09578:
https://arxiv.org/abs/2211.09578

It surveys consistent boundary formulations, tensor-network coarse graining, restricted spin foams and effective spin foams, including perturbative refinement limits.

**Use:** choose one concrete blocking map and define refinement consistency before adding larger graph sizes. Compare (i) direct fine amplitude, (ii) blocked fine amplitude, and (iii) coarse effective amplitude on the same boundary data.

## 5. Physical projector, master constraint and physical inner product

The finite-regulator identity already used in the repository is
\[
\mathbb M=\sum_{A,B}\hat C_A^\dagger G^{AB}\hat C_B,\qquad G>0,
\]
which gives
\[
\ker\mathbb M=\bigcap_A\ker\hat C_A
\]
when the domains and positivity assumptions are satisfied.

If zero is an isolated eigenvalue with gap Delta>0, then for positive self-adjoint M,
\[
P_0=E_{\mathbb M}(\{0\}),\qquad
\|e^{-t\mathbb M}-P_0\|\le e^{-t\Delta}.
\]
This is a finite spectral-filter estimate, not a continuum physical inner product. If zero is continuous spectrum, the normalized strong-limit projector may fail or vanish; a distributional rigging map is then a separate construction.

### Required gates

1. Prove symmetry and self-adjointness (or specify a self-adjoint extension) on a dense invariant domain.
2. Verify the complete constraint family, not a hand-selected subset.
3. Check anomaly closure/commutators, including graph-changing sectors and regulator dependence.
4. Measure the near-zero spectral density as cutoff/refinement grows; distinguish isolated zero eigenvalue from accumulating continuum spectrum.
5. Construct the rigging map on a specified test space and prove positivity/null quotient.
6. Compare finite master-projector matrix elements against the rigging-map limit on common test states.

The 2026 continuum-limit result makes step 5 central, not optional.

## 6. Relational physical clocks and quantum reference frames

Hausmann, Schmidhuber, Castro-Ruiz, **Measurement events relative to temporal quantum reference frames**, Quantum 9, 1616 (2025), DOI 10.22331/q-2025-01-30-1616:
https://doi.org/10.22331/q-2025-01-30-1616

The paper compares twirled-observable and purified-measurement Page-Wootters approaches. They agree for ideal clocks, but finite-resource clocks can yield time-nonlocal equations and, in one prescription, nonunitary evolution and ambiguity in temporal order; discreteness can restore unitarity/definite order under their assumptions.

### BQG clock contract

Do not define physical time as a graph coordinate or master-constraint heat-kernel parameter. Define a clock POVM E_T(tau), conditional states/observables, and clock quality. For each candidate graph clock, measure:

- monotonicity on the target semiclassical sector;
- clock resolution and energy/resource spread;
- clock-system correlations and backreaction;
- conditional norm conservation and complete positivity where applicable;
- temporal-order distinguishability;
- agreement of clock-neutral Dirac, Page-Wootters and reduced relational-Heisenberg descriptions in an ideal-clock control;
- stability across clock choices and overlapping clock patches.

Höhn, Smith and Lock's **Trinity of relational quantum dynamics** (PRD 104, 066001, 2021-09-01) is foundational but just outside the strict five-year interval. It establishes equivalence of the three descriptions under its quantum-reduction-map assumptions: https://doi.org/10.1103/PhysRevD.104.066001

## 7. Causal spin foams and retarded response

### 7.1 New causal Lorentzian vertex

Bianchi, Chen, Gamonal, **Causal spinfoam vertex for 4D Lorentzian quantum gravity**, Phys. Rev. D 113, 126020 (2026), DOI 10.1103/fwql-t4yr:
https://doi.org/10.1103/fwql-t4yr
Preprint: https://arxiv.org/abs/2601.23162

The construction encodes causal data in Toller matrices with
\[
T^{(+)}+T^{(-)}=D,
\]
and a Feynman i-epsilon representation. In the large-spin limit it selects Lorentzian Regge geometries compatible with its causal data and a single Regge-phase exponential.

Bianchi, Chen, Gamonal, **Toller matrices and the Feynman i-epsilon in spinfoams**, Phys. Rev. D 114, 046014 (2026), DOI 10.1103/v3kc-4n3n:
https://doi.org/10.1103/v3kc-4n3n

It gives analytic, integral-over-boost-spectrum and hypergeometric representations of the Toller matrices.

**Use:** prototype a causal orientation/lapse sector in the BQG spin-foam transition amplitude and compare the Feynman prescription against existing orientation and proper-time conventions. This can sharpen the causal-amplitude gate.

**Do not conflate:**
- Feynman propagator: time-ordered boundary condition;
- positive-lapse/causal transition amplitude: oriented histories;
- retarded Green function: response vanishes before the source in physical time;
- retarded commutator: physical-state expectation of a commutator multiplied by a physical-clock step function.

A Feynman i-epsilon prescription alone does not establish the retarded response.

### 7.2 Retarded physical kernel

After the physical inner product and clock are established, define complete observables A_T(tau), B_T(tau') and test
\[
G^R_{AB}(\tau,\tau')=
-i\,\theta(\tau-\tau')\,
\langle[A_T(\tau),B_T(\tau')]\rangle_{phys}.
\]
For a source-coupled system, compare this commutator expression with the derivative of the conditional expectation with respect to a source switched on at tau'. Check gauge invariance, state positivity, source normalization, adjoint relations and clock-choice stability.

## 8. Emergent causal cone and long-range graph couplings

### 8.1 Exponential interactions

Wilming and Werner, **Lieb-Robinson bounds imply locality of interactions**, Phys. Rev. B 105, 125101 (2022), DOI 10.1103/PhysRevB.105.125101:
https://doi.org/10.1103/PhysRevB.105.125101

For k-body lattice systems under their assumptions, exponential Lieb-Robinson bounds are equivalent to exponential spatial decay of interactions. Infer the admissible interaction-decay class from measured BQG commutator response instead of assuming locality.

A schematic exponential bound is
\[
\|[A_X(t),B_Y]\|\le C\|A_X\|\|B_Y\|
e^{-\mu(d(X,Y)-v_{LR}|t|)}.
\]
The constants and metric must be derived for the actual graph interaction and local Hilbert-space assumptions.

### 8.2 Long-range interactions

Sigal and Zhang, **On propagation of information in quantum many-body systems**, Annals of Physics 480 (2025) 170103, DOI 10.1016/j.aop.2025.170103:
https://doi.org/10.1016/j.aop.2025.170103

It proves maximal-velocity/light-cone approximation and Lieb-Robinson-type bounds for general lattice many-body systems, including long-range interactions.

McDonough, Yin, Lucas, Zhang, **Lieb-Robinson Bounds with Exponential-in-Volume Tails**, PRX Quantum 6, 040322 (2025), DOI 10.1103/n242-m1l2:
https://doi.org/10.1103/n242-m1l2

It refines the tail for operators spreading over volume-filling supports, with suppression of the form
\[
\exp\left[-\frac{(r-vt)^d}{(vt)^{d-1}}\right]
\]
in the stated regime and assumptions.

Chen, Lucas and Yin, **Speed limits and locality in many-body quantum dynamics**, Rep. Prog. Phys. 86 (2023), DOI 10.1088/1361-6633/acfaae:
https://doi.org/10.1088/1361-6633/acfaae

A technical review covering power-law interactions, bosons and information-theoretic speed limits.

### 8.3 BQG test

1. Derive effective graph couplings J_XY from the actual regulated operator, not from graph adjacency alone.
2. Measure decay by graph distance and interaction support size.
3. Check bounded local operator norms and summability assumptions.
4. Numerically compute commutator fronts for several operators/states/cutoffs.
5. Fit competing exponential, power-law and stretched-exponential envelopes with confidence intervals.
6. Apply only the theorem whose decay and boundedness hypotheses are satisfied.
7. Compare TT, scalar, vector and matter sectors using the same physical clock and physical retarded kernel.

A shared graph metric does not by itself imply a universal physical light cone.

## 9. Independent experimental falsification

### 9.1 GRB 221009A / LHAASO

LHAASO Collaboration, **Stringent Tests of Lorentz Invariance Violation from LHAASO Observations of GRB 221009A**, PRL 133, 071501 (2024), DOI 10.1103/PhysRevLett.133.071501:
https://doi.org/10.1103/PhysRevLett.133.071501

The paper reports 95% lower limits E_QG,1 > 10 E_Pl for linear dispersion and E_QG,2 > 6 x 10^-8 E_Pl for quadratic dispersion, within its assumptions and data treatment. Piran and Ofengeim, PRD 109, L081501 (2024), provide an independent time-of-flight analysis: https://doi.org/10.1103/PhysRevD.109.L081501

### 9.2 Intrinsic source lags

Vardanyan et al., **Revisiting tests of Lorentz invariance with gamma-ray bursts: Effects of intrinsic lags**, PRD 108, 123023 (2023), DOI 10.1103/PhysRevD.108.123023:
https://doi.org/10.1103/PhysRevD.108.123023

Their 49-burst Swift analysis shows source-lag modeling can shift inferred LIV limits by factors of a few. Bartlett et al. (2021) provide source-by-source Monte Carlo forward modeling and nuisance-noise marginalization: https://doi.org/10.1103/PhysRevD.104.103516

A hierarchical Bayesian GRB analysis published in 2026 also emphasizes intrinsic-lag modeling as the dominant systematic and reports no statistically significant LIV evidence under its spline-based model: https://doi.org/10.3847/1538-4357/ae3374

### 9.3 Formula and pipeline

For a commonly used phenomenological dispersion ansatz
\[
v_g(E,z)\simeq c\left[1-s\frac{n+1}{2}
\left(\frac{E(1+z)}{E_{QG,n}}\right)^n\right],
\]
the propagation delay between two observed energies is schematically
\[
\Delta t_{LIV}\simeq
s\frac{n+1}{2H_0}
\frac{E_h^n-E_l^n}{E_{QG,n}^n}
\int_0^{z_s}\frac{(1+z)^n}{h(z)}\,dz,
\qquad h(z)=H(z)/H_0.
\]
This formula depends on the chosen dispersion relation, sign convention and cosmological model; it is not a universal BQG prediction.

The observed delay must be modeled as
\[
\Delta t_{obs}=\Delta t_{source}+\Delta t_{propagation}
+\Delta t_{instrument}+\epsilon.
\]
Use source-level random effects/hierarchical inference, redshift uncertainties, energy-band response, detector systematics, and both linear and quadratic dispersion. Pre-register the likelihood and nuisance model before fitting BQG parameters. Include GW170817/GRB 170817A as a separate multimessenger constraint with source-emission-time uncertainty.

## 10. Prioritized implementation plan for the repository

### P0 — mathematical correctness before scale

1. Build a dense exact reference matrix for a frozen W5 charged basis.
2. Validate SU(2) representation action and all three commutators for every charged node.
3. Remove all state-dependent pruning from coherent operator assembly; quantify remaining fixed-mask error in operator norm.
4. Make pruning inventory transitive across imports, dynamic dispatch, caches and all state-merge paths; pair static inventory with runtime loss ledger.
5. Compare cache-on/off outputs and pruning ledgers byte-for-byte on deterministic frozen fixtures.

### P1 — actual dynamics and continuum

6. Benchmark graph-preserving versus graph-changing action on at least one exact 4-valent reference case.
7. Rerun W5/V10 cutoffs with a fixed basis hash, full outputs, no cancelled jobs and preserved artifacts.
8. Construct coarse maps and effective operators; measure RG flow of couplings and boundary amplitudes.
9. Test the 2026 distributional continuum/rigging-map architecture on a finite test space, including positivity and null quotient.
10. Revisit the V10/V11 common carrier only after the fine operator and coarse-map error budget are controlled.

### P2 — physical time and causal response

11. Define a graph-derived clock POVM and clock-resource/backreaction metrics.
12. Cross-check Dirac, Page-Wootters and reduced relational dynamics.
13. Add causal spin-foam orientation/Toller-matrix prototype as a separate amplitude experiment.
14. Derive the retarded commutator from physical observables and clock; do not relabel a Feynman amplitude.
15. Derive sector-resolved graph locality bounds using measured coupling tails.

### P3 — falsification

16. Generate dispersion predictions from the derived BQG retarded kernel.
17. Map to a declared EFT/LIV parameterization with explicit approximations.
18. Fit GRB 221009A and multi-GRB samples with intrinsic-lag nuisance models.
19. Include GW/EM multimessenger constraints and publish posterior predictive checks.
20. Report excluded parameter regions and null results; never describe a bound as a detection.

## 11. What this literature does NOT solve automatically

- No recent paper proves the BQG W5 operator is SU(2)-invariant.
- No paper makes state-dependent thresholding linear.
- The 2026 rigging-map construction is not automatically the rigging map for this repository's constraint family.
- A finite master-constraint projector is not a continuum physical Hilbert space.
- Toller/Feynman causal amplitudes are not automatically retarded Green functions.
- Lieb-Robinson bounds for bounded lattice systems cannot be imported until BQG local norm and interaction-decay hypotheses are checked.
- Published GRB bounds do not imply BQG predicts LIV; a dispersion law must first be derived from the physical kernel.
- A candidate RG fixed point is not proof that the BQG trajectory reaches it.

## Bibliography: core papers

1. Wilming & Werner (2022), PRB 105, 125101, https://doi.org/10.1103/PhysRevB.105.125101
2. Asante, Dittrich & Steinhaus (2022/2023), arXiv:2211.09578, https://arxiv.org/abs/2211.09578
3. Alesci, Mäkinen & Yang (2023), arXiv:2304.00268, https://arxiv.org/abs/2304.00268
4. Chen, Lucas & Yin (2023), Rep. Prog. Phys. 86, 116001, https://doi.org/10.1088/1361-6633/acfaae
5. Vardanyan et al. (2023), PRD 108, 123023, https://doi.org/10.1103/PhysRevD.108.123023
6. Bianchi & Martin-Dussaud (2024), Universe 10(4), 181, https://doi.org/10.3390/universe10040181
7. LHAASO Collaboration (2024), PRL 133, 071501, https://doi.org/10.1103/PhysRevLett.133.071501
8. Hausmann, Schmidhuber & Castro-Ruiz (2025), Quantum 9, 1616, https://doi.org/10.22331/q-2025-01-30-1616
9. Guedes et al. (2025), PRD 112, 026024, https://doi.org/10.1103/hjdk-kdhk
10. Sigal & Zhang (2025), Ann. Phys. 480, 170103, https://doi.org/10.1016/j.aop.2025.170103
11. McDonough et al. (2025), PRX Quantum 6, 040322, https://doi.org/10.1103/n242-m1l2
12. Bianchi, Chen & Gamonal (2026), PRD 113, 126020, https://doi.org/10.1103/fwql-t4yr
13. Han (2026), PRD 114, 044040, https://doi.org/10.1103/d8s7-jqfl
14. Bianchi, Chen & Gamonal (2026), PRD 114, 046014, https://doi.org/10.1103/v3kc-4n3n
15. Bruno et al. (2026), PRD 114, 066005, https://doi.org/10.1103/7493-9nb7
16. Assanioussi & Zeiß (2026 preprint), arXiv:2608.31152, https://arxiv.org/abs/2608.31152
17. Höhn, Smith & Lock (2021; just outside strict interval), PRD 104, 066001, https://doi.org/10.1103/PhysRevD.104.066001

**Final assessment:** the fastest defensible path is not to import a single grand new formula. Use the 2025 graph-changing implementation as a benchmark, the 2026 distributional continuum theorem to structure the physical Hilbert-space construction, the 2026 causal vertex as a separate causal-amplitude prototype, and modern Lieb-Robinson results only after deriving the actual BQG interaction tails. These methods reduce ambiguity and prevent false closure claims, but every BQG-specific gate still needs its own proof and reproducible calculation.
