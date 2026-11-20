# Black holes in BQG: horizon, entropy, evaporation and information recovery

**Branch:** \`research/bqg-infoton-black-holes-20261120\`  
**Status:** research specification and falsifiable reduced-model programme; no physical BQG black hole is claimed.

## 1. Non-negotiable consistency with BQG

A black hole is not introduced as a new microscopic constituent. It must be a state/phase of the same binary, SU(2), tetrahedral and relational-history degrees of freedom used elsewhere in BQG.

The dependency chain is

\[
\text{microscopic BQG amplitudes}
\to \text{physical constraints/history}
\to \text{Lorentzian metric sector}
\to \text{trapped region}
\to \text{horizon cut algebra}
\to \text{unitary radiation channel}.
\]

The current repository has an exact finite master-constraint kernel theorem, a combined C8 clock/system projector positive control, a finite metric-source functional, a reduced free TT propagator and several finite information/entanglement results. But the refinement-stable physical rigging map, normalized Lorentzian history amplitudes, physical inner product and causal metric remain open. The detailed stop conditions are in [BH-0: physical projector, Lorentzian history and metric](BH0_PHYSICAL_METRIC_GATE.md). Consequently, the programme distinguishes (A) exact finite-network statements, (B) conditional effective-theory checks and (C) physical BQG derivations still required.

## 2. Geometry and the horizon are derived, not inserted

### 2.1 Emergent metric

Use the existing local information-resistance construction only as a spatial graph metric:

\[
g_e^{QFI}=\tfrac14F_Q(\rho_e,K_e),\qquad
\ell_e=\ell_*\frac{g_*}{g_e^{QFI}},\qquad
d(i,j)=\min_{\gamma:i\to j}\sum_{e\in\gamma}\ell_e.
\]

This is not yet a Lorentzian metric. A physical lapse, shift, causal cone and extrinsic curvature must come from the same constraint algebra/history amplitudes. Do not identify the QFI graph distance alone with a spacetime metric.

### 2.2 Trapping horizon

Once a Lorentzian physical metric \(g_{\mu\nu}^{phys}\) exists, define a closed codimension-two surface \(S\) with outgoing/ingoing null normals \(l^\mu,n^\mu\), normalized by \(l\cdot n=-1\). Its expansions are

\[
\theta_{(l)}=q^{\mu\nu}\nabla_\mu l_\nu,\qquad
\theta_{(n)}=q^{\mu\nu}\nabla_\mu n_\nu.
\]

A future outer trapping horizon requires \(\theta_{(l)}=0\), \(\theta_{(n)}<0\), and \(\mathcal L_n\theta_{(l)}<0\). These are operational conditions on the emergent metric, not labels assigned to a high-density graph cut.

The event horizon is global and depends on the full future history; it cannot be defined from one instantaneous slice alone. The programme therefore starts with trapping horizons and only later asks whether a global event horizon exists.

### 2.3 Required Schwarzschild/Kerr correspondence

In a controlled semiclassical regime, the BQG effective metric must reproduce the appropriate vacuum exterior and its invariant observables. For a nonrotating benchmark only, the comparison target is

\[
ds^2=-\left(1-\frac{r_s}{r}\right)c^2dt^2
+\left(1-\frac{r_s}{r}\right)^{-1}dr^2+r^2d\Omega^2,
\quad r_s=2GM/c^2.
\]

This is a correspondence test, not an assumption about the microscopic BQG solution. The rotating sector must be tested separately against Kerr; spherical symmetry cannot be used to claim Kerr closure.

## 3. Horizon area and entropy

The physical area must be computed from the emergent metric or its validated discrete area operator:

\[
A_H=\int_S\sqrt{\det q}\,d^2x
\quad\text{or}\quad
\widehat A_H=\sum_{f\cap S}\widehat a_f.
\]

The horizon entropy must be computed from the physical state and the horizon/boundary algebra:

\[
S_H=-\operatorname{Tr}(\rho_H\ln\rho_H),
\]

with gauge constraints and edge-mode/boundary contributions treated consistently. A tensor-network cut bound

\[
S(\rho_L)\leq\sum_{e\in C}\ln\chi_e
\]

is only an upper bound. It is not the Bekenstein–Hawking law.

The target correspondence is

\[
S_H\longrightarrow \frac{k_B A_H}{4\ell_P^2}
\]

in the semiclassical limit. BQG must derive both the area scaling and the coefficient \(1/4\) from the microscopic measure/physical state. The coefficient may not be inserted by choosing a virtual bond dimension or fitting a boundary degeneracy after the fact.

Required checks: area scaling across graph refinements; independence from arbitrary cut placement; gauge-invariant boundary algebra; regulator independence; logarithmic/subleading corrections; comparison of microcanonical and entanglement entropy.

## 4. Hawking channel and evaporation

The outgoing flux must be derived from the physical Lorentzian state and its asymptotic mode decomposition. In the semiclassical correspondence regime it must reproduce the Hawking temperature

\[
T_H=\frac{\hbar c^3}{8\pi G M k_B}
\]

for Schwarzschild, including the greybody-filtered spectrum at infinity. A near-horizon Bogoliubov calculation can be used as a controlled bridge only after the BQG causal structure and mode normalization are established.

The mass-loss law is then a consequence of the computed flux, not a separately imposed rule:

\[
\frac{dM}{dt}=-\frac{P_{\rm out}(M)}{c^2}.
\]

The leading dimensional benchmark is \(P\propto M^{-2}\), hence \(t_{\rm evap}\propto M^3\); the numerical coefficient depends on particle content and greybody factors and must be computed, not copied as a BQG result.

## 5. Information conservation and the Page curve

The total physical evolution must be unitary on the physical Hilbert space. At retarded time \(u\), define the radiation subsystem \(R(u)\) using asymptotic outgoing modes and compute

\[
S_R(u)=-\operatorname{Tr}\rho_R(u)\ln\rho_R(u).
\]

The required qualitative pattern for an initially pure collapsing state is:

1. early radiation entropy increases;
2. it reaches a maximum near the Page time;
3. it decreases as correlations transfer to the radiation;
4. the final radiation state is pure (up to the chosen asymptotic sector and any explicitly retained remnant).

The Page curve is not obtained by plotting the minimum of two guessed curves. It must be calculated from the same microscopic BQG amplitudes that generate the horizon and outgoing flux. A unitary random-circuit model is a useful numerical control, but not a derivation from BQG.

The central test is simultaneous consistency of:
- total-state norm and constraint preservation;
- positive reduced density matrices;
- radiation entropy and mutual information;
- energy balance between black hole, radiation and any remnant;
- regulator/refinement stability;
- late-time purity without discarding a hidden sector.

## 6. Information-paradox resolution criteria

The project must state which mechanism actually resolves the paradox. Candidate mechanisms are not interchangeable:

* **Unitary evaporation:** information is encoded in radiation correlations.
* **Remnant:** a residual physical sector stores information; its state count and stability must be derived.
* **Non-unitary evolution:** requires explicit modification of the physical evolution law and quantitative consistency tests.
* **Holographic reconstruction:** requires a proved boundary algebra/code and reconstruction map, not just a cut entropy inequality.

No mechanism is selected by this document. The BQG calculation must decide.

## 7. Staged gates and stop conditions

| Gate | Required output | Pass condition |
|---|---|---|
| BH-0 | Physical projector, Lorentzian history and metric | anomaly-controlled constraints; nontrivial positive physical inner product; normalized causal amplitudes; derived Lorentzian metric and TT pole; refinement stability |
| BH-1 | Trapping horizon | stable zero of outgoing expansion; convergent under refinement |
| BH-2 | Area operator | gauge-invariant spectrum and continuum scaling |
| BH-3 | Horizon entropy | \(S/A\to k_B/(4\ell_P^2)\) without fitted coefficient |
| BH-4 | Hawking channel | normalized outgoing spectrum and flux |
| BH-5 | Evaporation | energy-conserving mass loss from computed flux |
| BH-6 | Page curve | entropy from full radiation density matrix, late-time purity |
| BH-7 | Robustness | same results across regulators, refinements and initial data |

A failure at an earlier gate blocks claims at every later gate. In particular, no Page curve is called a BQG result before BH-0 through BH-5 pass.

## 8. Relationship to existing results

* \`research/holography_20261109/README.md\`: finite cut-entropy bound and parity code; explicitly not a physical BH entropy derivation.
* \`BQG_LOCAL_INFORMATION_RESISTANCE_METRIC.md\`: local QFI resistance graph metric; spatial construction only.
* \`TT_PROPAGATOR_FIRST_PASS.md\`: exact reduced free TT transfer; not full physical-history propagator.
* \`TT_VACUUM_TWO_POINT_RESULT.md\`: reduced Gaussian TT vacuum has \(P(k)\sim k^{-1}\); observer smoothing exponent is not the vacuum spectrum.
* \`docs/information_infoton_review/oldbridge.md\`: finite \(j=2\) and conditional resonance bridge; no physical infoton established.

## 9. Completion definition

The black-hole sector is complete only when one frozen microscopic BQG model, without switching models between stages, passes BH-0 through BH-7 and reproduces the semiclassical exterior, area law, Hawking flux, energy balance and unitary Page curve in a common continuum regime.

Until then this document is a precise construction plan and falsification contract, not a claim that BQG has solved black-hole physics.
