# BQG physical environment, Born measure and causal-response frontier

**Status: finite reference gates added; all theory-specific physical gates remain `open_physical`.**

This note separates exact conditional mathematics from the missing construction of physical inputs by the graph-changing BQG constraints.

## 1. Physical environment overlap and positivity

Let \\(\mathcal H_{phys}\\) be the physical state space after constraint imposition and quotienting null states. Let its physical sesquilinear form be represented at a regulator by a Hermitian positive-semidefinite matrix \\(\eta_{phys}\\). For physical environment records \\(|E_a\rangle\\), define

\\[
G_{ab}=\langle E_a|E_b\rangle_{phys}.
\\]

For every complex vector \\(c\\),

\\[
c^\dagger Gc=
\left\|\sum_a c_a|E_a\rangle\right\|_{phys}^{2}\geq0.
\\]

Thus the overlap matrix is Hermitian positive semidefinite, and

\\[
|G_{ab}|^2\leq G_{aa}G_{bb}.
\\]

This is an exact Gram-matrix theorem conditional on a positive physical inner product. Null vectors must be quotiented; zero-norm representatives cannot count as distinct physical records. The new engine computes the Gram matrix from a supplied physical metric, tests positivity/Cauchy-Schwarz, and rejects indefinite forms.

**What is not yet done:** no actual BQG physical environment states or theory-specific rigging-map metric have been produced. The finite tests use a supplied positive metric. A kinematical tensor product \\(H_S\otimes H_E\\) and a naive partial trace are not automatically physical in a constrained theory. The environment algebra/factorization must be gauge invariant (or justified through boundary/edge modes), and overlaps must be evaluated after the physical quotient.

## 2. Born rule: conditional theorem versus BQG axiom

There are two logically distinct claims.

### Conditional Hilbert-space result

If the theory supplies a complex physical Hilbert space with positive inner product, normalized states, a physical PVM/appropriate POVM, and a normalized noncontextual additive probability measure on orthogonal alternatives (with the dimensional/domain assumptions of Gleason's theorem or a suitable POVM extension), then the probability measure has trace form \\(p(P)=\operatorname{Tr}(\rho P)\\). For a pure state and rank-one event this is the squared physical amplitude.

Likewise, once physical class operators \\(C_\alpha\\), a physical state \\(\rho_{phys}\\), and the physical inner product are independently constructed,

\\[
D(\alpha,\beta)=\operatorname{Tr}_{phys}
(C_\alpha\rho_{phys}C_\beta^\dagger)
\\]

is a strongly-positive decoherence functional, and \\(\mu(A)=D(A,A)\\) is quadratic. Ordinary additive probabilities for a history family require decoherence (vanishing or controlled off-diagonal terms) on that family.

### BQG status

The candidate BQG has not yet derived the physical history class operators, physical inner product/rigging map, complete physical measurement algebra, or the required probability assignment. Therefore the honest verdict is:

- Quadratic probabilities follow conditionally from the stated Hilbert-space measurement assumptions.
- A quadratic history measure follows once the physical decoherence functional is independently defined.
- **The Born rule has not yet been derived from full BQG.** At the current BQG level it remains an additional probability postulate/interpretation until the required physical structure and assumptions are derived.

Sorkin \\(I_3=0\\) is a consistency consequence of a bilinear decoherence functional and quadratic measure, not an independent Born-rule proof.

## 3. Relational clock and physical projector

At finite regulator the existing master-constraint theorem is

\\[
M_\epsilon=\sum_{A,B}C_A^{(\epsilon)\dagger}
G^{AB}C_B^{(\epsilon)}\geq0,
\qquad
P_\epsilon=\mathbf 1_{\{0\}}(M_\epsilon),
\\]

when zero is an isolated spectral sector. This selects the common kernel; it does **not** create physical time. The heat-kernel parameter in \\(e^{-sM_\epsilon}\\) is a spectral regulator parameter, not a clock.

The existing C8 Page-Wootters gate proves that a combined clock+system constraint can preserve conditional evolution while clock-only averaging erases it. But that clock is externally introduced and its system step \\(R=J\\) is a positive-control operator, not the physical graph-changing gravitational evolution.

Two legitimate routes remain:

1. **Derived deparametrized clock:** derive a matter/boundary clock \\(T\\), prove monotonicity in the regime and suitable self-adjointness/domain conditions, and derive \\(C=P_T+H_{rel}\approx0\\). Only then may one obtain \\(i\hbar\partial_T\Psi=H_{rel}\Psi\\).
2. **Relational-observable/boundary-history route:** if no global deparametrization exists, construct complete observables or boundary-history amplitudes. A conditional physical state is schematically \\( |\psi;\tau\rangle_{phys}\propto E_T(\tau)P_{phys}|\psi_{kin}\rangle\\), normalized with the physical inner product, with a proof that the observables preserve the constraints.

In either route the normalized physical matrix elements must converge under graph refinement. No arbitrary clock may be chosen just to manufacture a desired propagation cone.

## 4. Physical retarded response

Only after clock and physical observables are derived can one define

\\[
G^R_{AB}(\tau)=\frac{i}{\hbar}\theta(\tau)
\langle[\mathcal O_A(\tau),\mathcal O_B(0)]\rangle_{phys}.
\\]

The expectation must use the physical state/inner product. Operators must be gauge-invariant relational or properly dressed boundary observables. An unprojected kinematical commutator, symmetric correlator, mutual information, master-constraint heat-kernel parameter, or constraint resolvent spectral variable is not a substitute.

## 5. One cone across photon, matter and gravity

At every graph refinement \\(\ell\\), compute channel-resolved physical retarded functions with the same clock prescription, physical distance map and boundary conditions:

\\[
G^R_{a,\ell}(x,y;\tau),\qquad
a\in\{\gamma,\mathrm{matter},TT,scalar,vector\}.
\\]

For each level:
- determine each response frontier with stated threshold, precision and uncertainty;
- test exact zero support where analytically possible;
- otherwise derive a rigorous Lieb-Robinson-type bound, without calling it exact microcausality;
- extract channel speeds separately;
- compare all channels without averaging away disagreements;
- extrapolate each channel as the graph scale is removed.

The target is

\\[
\lim_{\ell\to\infty}v_{a,\ell}=v_*
\quad\text{for every physical channel }a.
\\]

A finite sampled threshold frontier is only a diagnostic. It cannot prove exact microcausality, physical distance calibration, or a continuum common cone. The new engine labels the refinement test `FINITE_REFINEMENT_DIAGNOSTIC_ONLY`.

## 6. New fail-closed physical gates

All four are added to `physicalization_gates.json` and must remain open until theory-specific evidence exists:

- `PHYSICAL_ENVIRONMENT_OVERLAP`: derive the physical environment algebra/factorization and evaluate its Gram matrix after the physical quotient.
- `PHYSICAL_BORN_MEASURE`: derive the physical history measure and establish the assumptions under which quadratic probabilities follow; otherwise state the Born map as an explicit additional axiom.
- `PHYSICAL_RELATIONAL_CLOCK_RESPONSE`: derive the clock, rigging map/projector, relational observables and physical retarded functions from the actual graph-changing constraints.
- `PHYSICAL_COMMON_CAUSAL_CONE`: establish common photon/matter/gravity fronts with refinement stability and exact-support or rigorous bounded-tail semantics.

Passing the new reference engine cannot freeze any of these gates.

## 7. Actual finite tests

The new CI reference engine tests:
- physical Gram positivity and Cauchy-Schwarz for a supplied positive metric;
- rejection of an indefinite purported physical metric;
- two nonzero history branches reaching the same system output, with orthogonal environment records suppressing the cross term;
- survival of the cross term without a distinguishable environment record;
- support depth versus coherent depth after path cancellation;
- photon/matter/tensor front comparison at several synthetic refinement levels.

These are finite reference controls. The synthetic refinement data are not physical BQG predictions.

## 8. Remaining production blockers

1. Actual all-node W5 K-K-V branches must feed physical history class operators.
2. Actual BQG physical environment states and the rigging-map inner product are absent.
3. The C8 clock remains an external positive control.
4. No theory-specific physical retarded photon/matter/TT response has been computed.
5. No physical response data across graph refinements exist.

The next scientific closure step is to connect W5 operator outputs to physical history class operators and a derived rigging map. Until then, all four new physical gates stay `open_physical`.
