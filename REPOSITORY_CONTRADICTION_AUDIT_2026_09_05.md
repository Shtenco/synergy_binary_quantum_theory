# Repository contradiction audit — 2026-09-05

Status: **repo-wide scientific consistency audit. This document does not close any physical gate and does not claim dark matter or dark energy.**

This audit compares the current canonical `main` status surface with the active physical-scalar research branch and the newly completed regulator-safe Euclidean Peter–Weyl diagnostics.

## 1. Executive result

No fatal algebraic contradiction has been identified in the current exact/finite q=2 -> geometry -> Peter–Weyl -> GR/HDA-control chain.

However, the repository contains several **active status/semantics contradictions** left by earlier scientific eras. The most serious are not numerical: they concern the word `closed` and whether an oracle/reference reconstruction is being described as a dynamically generated end-to-end theory.

The correct current scientific statement is:

```text
structural dictionaries / finite controls          = extensive and largely internally consistent
microscopic end-to-end generative dynamics         = NOT CLOSED
full graph-changing physical projector/history     = OPEN PHYSICAL
connected interblock physical history              = OPEN PHYSICAL
physical TT kernel                                  = OPEN PHYSICAL
physical Maxwell kernel                             = OPEN PHYSICAL
physical background/scalar cosmology                = OPEN PHYSICAL
BQG dark matter                                     = NOT DERIVED
BQG dark energy                                     = NOT DERIVED
experimental confirmation                           = NO
```

A physical projector is **not** a procedure for projecting the theory toward DM/DE. Its role is only to impose the quantum constraints / physical-history quotient. A dark-sector interpretation is allowed only after an independently derived gauge-reduced physical effective action has the required infrared behavior.

---

## 2. Newly completed numerical cross-checks

### 2.1 Five-node Euclidean boundary master

Preregistered definition:

\[
M_B[i,j]=\sum_{v=0}^4\langle H_v^E b_i\mid H_v^E b_j\rangle.
\]

At `Jmax=5/2`, on the frozen 32-dimensional all-`j=1/2` q=2 boundary:

```text
dimension                         = 32
rank                              = 32
nullity                           = 0
lambda_min                        = 9.651811183254074
lambda_max                        = 14.48385071910081
condition number                  = 1.500635522608476
Hermiticity error                 = 0
first-action boundary return      = 0 by parity
node-trace relative spread        = 4.03e-14
quadratic-form regression error   <= 7.01e-16
```

Scientific meaning:

\[
\boxed{
\bigcap_{v=0}^{4}\ker H_v^E\cap\mathcal H_{q=2,32D}=\{0\}.
}
\]

This strengthens the older two-node/summed-operator 32D no-go. It does **not** prove that the enlarged graph-changing Peter–Weyl habitat has no physical states.

### 2.2 First outgoing Peter–Weyl layer

For

\[
g_{v i}=H_v^E b_i
\]

the full labelled Gram matrix has

```text
labelled columns                  = 160
rank(G1)                          = 160
nullity(G1)                       = 0
generated odd dimension           = 160
total K1 dimension by parity      = 32 + 160 = 192
G1 eigenvalue min                 = 0.4633234351993575
G1 eigenvalue max                 = 8.652092493235065
cross-node Gram weight fraction   = 0.41313794379358565
boundary projection max norm      = 0
pair reconstruction error         <= 2.78e-16
```

Important interpretation correction:

\[
\boxed{\lambda_{min}(G_1)=0.4633\text{ is NOT a constraint/master gap.}}
\]

It is an eigenvalue of the overlap/Gram geometry of the generated basis. It must not be compared with the boundary-master `9.6518` or interpreted as a near-zero physical mode.

The actual enlarged master gap requires the separate K1 Ritz calculation

\[
D c=\lambda G_1 c,
\qquad
D_{\alpha\beta}=\sum_w\langle H_wg_\alpha|H_wg_\beta\rangle,
\]

which was still computing at the time of this audit.

### 2.3 Consistency with the old 32D baseline

The retained earlier operator

\[
K_{32}=(H_{E,0}+H_{E,1})^\dagger(H_{E,0}+H_{E,1})
\]

also had

```text
rank     = 32
nullity  = 0
lambda_min = 4.306075987001585
```

There is no contradiction. The operators are different: the old object uses one summed constraint with node cross terms; the new master sums five positive node norms without allowing cancellation between constraint labels. Both agree on the only invariant conclusion relevant here: **the bare 32D q=2 carrier has no Euclidean common normal-constraint zero.**

---

## 3. Direct contradiction C1 — stale physicalization closure language

`PHYSICALIZATION_MASTER_CLOSURE_2026_08_17.md` states:

```text
The architecture between microscopic dynamics and experiment is now closed.
No undefined “RG magic” remains between these steps.
```

This is incompatible with the newer canonical truth surface, which keeps open:

```text
PHYSICAL_PROJECTOR_HISTORY
CONNECTED_INTERBLOCK_HISTORY
PHYSICAL_TT_KERNEL
IR_SIX_VECTOR
COMMON_SCALE_CALIBRATION
DYNAMICAL_MAXWELL_KERNEL
PHYSICAL_BACKGROUND_COSMOLOGY
PHYSICAL_SCALAR_COSMOLOGY
LENSING_DYNAMICS_CLOSURE
```

Classification: **DIRECT STATUS CONTRADICTION / SUPERSEDED ACTIVE DOCUMENT.**

Required fix: move the 2026-08-17 document to `docs/archive/` or add an unambiguous top banner:

```text
SUPERSEDED HISTORICAL FRONTIER — NOT A CURRENT CLOSURE CLAIM
```

and point to `THEORY_STATUS.md` + `physicalization_gates.json`.

---

## 4. Direct contradiction C2 — `C6_PHYSICAL_KERNEL_CLOSURE.md` is named too strongly

The six-dimensional parity-even quartic TT quotient and six-observable extractor are exact algebraic results. That is a genuine closure of the **observable dictionary**.

But the same file is titled `Physical C6 closure` and describes a `physical` six-Wilson momentum sector while the current truth ledger says the interacting physical TT kernel and physical six-vector are still open.

Classification: **STATUS/INTERPRETATION CONTRADICTION, not an algebraic error.**

Correct wording:

```text
C6 TT OBSERVABLE-DICTIONARY CLOSURE
```

not

```text
PHYSICAL KERNEL CLOSURE.
```

The exact six-dimensional quotient survives unchanged.

---

## 5. Direct contradiction C3 — `core_theory_closed_declared=true` versus open generative arrows

`theory_gates.json` currently declares

```text
core_theory_closed_declared = true
```

and `CANONICAL_THEORY_PACKAGE.md` says `CORE THEORY PACKAGE: CLOSED IN DECLARED SCOPE`.

But `BIT_TO_SPACETIME_CENTRAL_EQUATION.md` and the deprecated final certificate explicitly retain as open:

```text
microscopic dynamics -> required geometric qubit/two-form state
bare graph rule -> dynamically selected global topology
full graph-changing off-shell HDA
uniform regulator/refinement removal
physical scale
blind physical prediction
```

The Plebanski/Urbantke/S4 reconstruction is a valid oracle/reference reconstruction but does not prove that the frozen binary microscopic dynamics dynamically generates the required B/two-form state.

Classification: **SEMANTIC/GENERATIVE CLOSURE CONTRADICTION.**

Recommended replacement status schema:

```text
structural_control_package_complete = true
generative_micro_to_geometry_closed = false
full_constraint_dynamics_closed = false
physical_theory_closed = false
experimentally_confirmed = false
```

This preserves all finite theorems while removing the misleading implication that the complete theory has been dynamically derived end-to-end.

---

## 6. Scope drift C4 — HDA `PASS` versus HDA `OPEN`

There is no mathematical contradiction if scopes are kept explicit:

```text
finite three-node graph-changing HDA control = TESTED_FINITE
arbitrary-graph/full Lorentzian physical HDA = OPEN
```

Some active documents still use only `multi-node HDA OPEN`, while others say `HDA structural closure`.

Classification: **DOCUMENTATION SCOPE DRIFT.**

Required fix: every status surface should carry the two separate flags above.

---

## 7. Foundational issue C5 — starting complex qubit versus emergent complex structure

The canonical theory begins with

\[
|\psi\rangle=\alpha|0\rangle+\beta|1\rangle,
\qquad \alpha,\beta\in\mathbb C,
\]

and later derives a real history/orientation operator

\[
J^2=-I,
\qquad
W(\theta)=e^{-\theta J}.
\]

These are compatible if BQG is explicitly **quantum gravity from qubits**.

They become circular only if the same result is advertised as a derivation of complex quantum mechanics itself from pre-quantum classical bits.

Classification: **FOUNDATIONAL CLAIM FORK / POTENTIAL CIRCULARITY.**

Recommended current canonical programme:

```text
quantum two-level carrier is assumed;
C4/C8 derives an internal geometric/history realization of complex phase;
Born rule is not derived;
physical photon/Maxwell dynamics remains open.
```

A stronger `classical bits -> full QM` programme must be a separate research branch and may not use complex qubit/SU(2) machinery upstream of the derivation.

---

## 8. Scope issue C6 — selected topology versus dynamically forced topology

The 16-cell boundary is an exact/stable selected PL `S3` completion of the local q=2 shell in the tested refinement family.

The bare microscopic rule has not uniquely forced that global completion.

Classification: **NO MATHEMATICAL CONTRADICTION; STRONGER UNIQUENESS CLAIM REMAINS OPEN.**

Any arrow written as

```text
q=2 -> unique global S3
```

must be replaced by

```text
q=2 local shell -> selected canonical global S3 completion with tested stability.
```

---

## 9. Scope issue C7 — smoothing is not a physical scalar propagator

The `b^-2`, `b^-3`, `b^-4` smoothing hierarchy depends on weak-correlation/self-averaging assumptions.

It is a coarse-reconstruction law, not a physical connected two-point function. Long-range correlations can alter it.

Classification: **CONDITIONAL CONTROL, not contradiction.**

Forbidden inference:

```text
b^-2 smoothing -> 1/k^2 Poisson propagator -> dark matter.
```

---

## 10. DM/DE clarification C8 — what the projector actually does

The physical/master projector is not a dark-sector generator.

Its logical role is

\[
\mathcal H_{kin}
\xrightarrow{\text{constraints / rigging}}
\mathcal H_{phys}.
\]

Only after a physical history measure exists may one build a source response / effective action.

A dark-matter-like outcome would require, after gauge reduction, either

- a modified long-range scalar constraint/source response that passes dynamics+lensing+growth+cluster tests; or
- an additional stable physical pole with correct degree count, residue, sound speed, clustering and abundance.

A dark-energy-like outcome requires a homogeneous background contribution to the same physical history effective action, preserving relative absolute amplitudes between geometries and deriving `rho_hist(a), p_hist(a)` before `w(a)`.

If neither structure appears, the result is simply:

```text
BQG does not explain DM/DE in the tested physicalization branch.
```

No projector depth, spectral window, source normalization or basis choice may be tuned to force the desired phenomenology.

---

## 11. Gauge-order correction C9

Lapse and shift are constraint multipliers in ADM, not ordinary propagating operators on the final physical Hilbert space.

Two legal routes remain:

1. background/off-shell `Z[g,T] -> Gamma[g]` followed by Ward/Dirac/Schur constraint reduction;
2. fully gauge-invariant relational observables after projection.

Forbidden route:

```text
P_phys -> ordinary lapse operator -> invert correlator -> call pole dark scalar.
```

Classification: **CORRECTION CONSISTENT WITH THE NEW CANONICAL BRIDGE.**

---

## 12. Dark-energy normalization correction C10

For fixed-background connected correlators, `Z_corr[0]=1` is legitimate.

For vacuum/background gravity, separately imposing `Z[g;0]=1` for every geometry erases the geometry-dependent zero-source free energy and can destroy the very cosmological-volume term being tested.

A safe object is instead

\[
\Delta W[g;g_{ref}]
=-\log\frac{\mathcal Z_{phys}[g]}{\mathcal Z_{phys}[g_{ref}]}
\]

with one geometry-independent measure normalization convention.

Classification: **IMPORTANT NORMALIZATION FIREWALL; consistent with `CONSTANTS_ZERO_FIT_LEDGER.md`, which already states that the cosmological constant is not derived.**

---

## 13. Lorentzian status C11

The genuine full ordered `epsilon^{abc}` Lorentzian amplitude remains a computational no-result at the time of this audit; the active 8-shard rerun is still inside the amplitude step.

Therefore the following are forbidden at present:

```text
full Lorentzian physical kernel closed
Lorentzian zero/nonzero logical return measured
physical scalar pole from H_L
physical Maxwell/gravitational common cone derived from H_L
```

Classification: **OPEN / COMPUTING, not zero and not nonzero.**

---

## 14. Research-branch merge hygiene C12

The active research branch adds important physicalization corrections, but the canonical truth ledgers (`THEORY_STATUS.md`, `theory_gates.json`, `physicalization_gates.json`, `q2_scalar_frontier.json`) have not yet been synchronized to cite all of them.

This is acceptable while the work remains a research branch.

It becomes a contradiction risk if the branch is merged without simultaneously updating the canonical status hierarchy.

Recommended precedence after cleanup:

```text
1. physicalization_gates.json          physical truth
2. q2_scalar_frontier.json             scalar child truth
3. theory_gates.json                   structural truth, after closure rename
4. THEORY_STATUS.md                    human-readable summary
5. README.md                           canonical narrative
6. dated/superseded documents          historical only
7. docs/archive/*                      archaeology only
```

---

## 15. Final audit verdict

### No fatal contradiction found in these retained results

```text
q=2 selector in declared family
Walsh tetrahedral frame
rank-two X/Z metric tangent and Y orientation no-go
selected/stable 16-cell completion
exact d*=3 causal-volume fixed point
finite Regge/HDA/Peter-Weyl controls in their declared scopes
six-dimensional parity-even S4 quartic TT quotient
six-observable extractor
finite master/projector theorems
current Euclidean 32D no-go and K1 span results
```

### Real contradictions requiring canonical cleanup

```text
A. stale 2026-08-17 “microscopic dynamics -> experiment closed” claim
B. `physical C6 closure` wording versus open physical TT kernel
C. `core_theory_closed_declared=true` versus explicitly open generative arrows
D. HDA status wording that fails to distinguish finite three-node PASS from full arbitrary-graph OPEN
```

### Foundational / interpretational blockers, not algebraic failures

```text
E. complex-qubit start versus any claim to derive full QM from classical bits
F. selected topology versus unique dynamical topology
G. smoothing exponent versus physical scalar propagator
H. Hamiltonian-only low mode versus full Dirac physical mode
I. normalized local source Hessian versus vacuum/background action
J. constraint Gram eigenvalue versus physical mass/frequency/dark mode
```

The strongest current physics conclusion is therefore deliberately modest:

\[
\boxed{
\text{BQG has a substantial internally checked structural candidate architecture,}
}
\]

but

\[
\boxed{
\text{the full generative/physical theory is not yet closed and DM/DE are not derived.}
}
\]

That conclusion is stronger scientifically than forcing a dark-sector interpretation, because the current Euclidean calculations already eliminate the simplest bare-q=2 shortcut and require any future effect to survive genuine graph/representation/history dressing, full constraint reduction, physical source response and refinement.