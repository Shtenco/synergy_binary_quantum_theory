# BQG closure hierarchy and current RG frontier

Date: 2026-11-10
Branch: `main`

## Purpose

This ledger consolidates the repository's already established self-consistency and closure results. It prevents a local open numerical test from being mistaken for a failure of the whole construction, and prevents a finite zero or algebraic identity from being promoted to a continuum physical theorem.

The central rule is to keep four logically different statements separate:

1. an exact algebraic/representation identity;
2. a finite-regulator numerical null or rank certificate;
3. a theorem that converts a controlled residual into physical-projector control;
4. an actual sequence of BQG residuals that tends to zero under refinement.

The first three are present in several branches of the project. The fourth is still the decisive RG task.

## 1. Already closed structural and self-consistency chain

### 1.1 Binary selector and three-dimensional causal-volume fixed point

Source: `Q2_DIMENSION3_FIXED_POINT_CLOSURE.md`.

The local selector

```math
q+2=2^q
```

selects the integer solution (q=2). With the frozen route rewrite, (B=2^q=4), active causal edges grow by (2B=8), and causal distance doubles. The exact vertex count is

```math
N_g=(4\,8^g+10)/7.
```

The finite-step exponent increases monotonically to

```math
d_g\nearrow d_*=\log_2 8=3.
```

This is an exact asymptotic result for the declared route rewrite, not a universal theorem for every graph ensemble.

### 1.2 Canonical geometry and gravity structure

Source: `PHYSICALIZATION_MASTER_CLOSURE_2026_08_17.md`.

The ledger records independent topology/scaling routes to the 3+1 interpretation, the DeWitt/HDA selection (c=1/2), the onsite (S_4) decomposition, the complete six-dimensional parity-even quartic TT Wilson basis, and the exact-rational blind extractor. These are structural/algebraic closures in their stated domains.

The six microscopic interblock transfer amplitudes and their refinement/locality controls remain numerical tasks. The algebraic map to observables is closed; the physical six-vector is not yet computed and frozen.

### 1.3 Canonical representation refinement

Source: `BQG_CANONICAL_REFINEMENT_RESIDUAL_DECOMPOSITION.md`.

The first (j=1/2\to j=1) carrier map is fixed by the unique (S_4[2,2]) intertwiner

```math
W=\begin{pmatrix}0&2/3\\1&0\\0&-\sqrt5/3\end{pmatrix},
\qquad W^\dagger W=I_2.
```

This embedding is not a fitted rotation. Its symmetry and isometry properties are fixed.

For arbitrary fine/coarse operators, the exact residual splits into orthogonal pieces:

```math
R=M_cW-WM_f=L+WD,
```

```math
L=(I-WW^\dagger)M_cW,
\qquad D=W^\dagger M_cW-M_f,
```

and

```math
R^\dagger R=L^\dagger L+D^\dagger D.
```

Thus exact intertwining requires both (L=0) and (D=0). If the fine operator is defined as the Galerkin pullback (W^\dagger M_cW), then (D=0) by construction and the remaining defect is precisely carrier leakage (L). This special case must not be confused with an independently generated fine BQG operator.

## 2. Already proved finite operator theorems

### 2.1 Master residual controls physical-projector residual

Source: `BQG_REFINEMENT_PROJECTOR_STABILITY_THEOREM.md`.

For positive finite master operators with isolated zero sectors and gaps (Delta_f,Delta_c),

```math
\|P_c\iota-\iota P_f\|
\le
\frac{\|M_c\iota-\iota M_f\|}
{\min(\Delta_f,\Delta_c)}.
```

This is an exact finite-dimensional theorem. It does not assert that the actual BQG residual decays.

### 2.2 Scale-free comparison

Source: `BQG_SCALE_FREE_REFINEMENT_RESIDUAL_THEOREM.md`.

Independent positive rescalings of the two master operators are removed by

```math
\epsilon^{sf}
=
\inf_{\rho>0}
\frac{\|M_c\iota-\rho\,\iota M_f\|}
{\min(\Delta_c,\rho\Delta_f)}.
```

This is the correct scale-independent physical-projector control. A post-hoc Frobenius scalar fit is only a diagnostic and is not automatically the minimizing (ho) in this gap-normalized operator-norm expression.

### 2.3 Graph-changing constraints require two-sided habitat refinement

Source: `BQG_TWO_SIDED_HABITAT_REFINEMENT_THEOREM.md`.

For one-hit constraint maps and embeddings into the graph-changing image habitat,

```math
E_a=C_{c,a}\iota_0-\iota_1C_{f,a},
\qquad
F_a=C_{c,a}^\dagger\iota_1-\iota_0C_{f,a}^\dagger,
```

the full master residual has the exact identity

```math
R=\sum_a(C_{c,a}^\dagger E_a+F_aC_{f,a}).
```

The forward residuals (E_a) control the pullback/internal mismatch; the return residuals (F_a) are also needed to control the full off-carrier master leakage. This is the correct route for graph-changing BQG, rather than comparing only base-carrier return kernels.

## 3. Finite zero/rank/closure certificates already present

### 3.1 Depth-6 finite master sector

Source: `BQG_DEPTH6_2111_MASTER_CLOSED_2026-10-01.json`.

The declared (K_5), depth-6 sectors ([1^5]) and ([2,1,1,1]) are closed with zero common kernel in the tested representation sectors. For ([2,1,1,1]), the recorded dimension and rank are both 103318. The exact 16-dimensional (H_0)-zero subspace is separately killed injectively by (H_1).

Scope is finite depth-6 and the named irreps. The other depth-6 irreps and the full all-depth/continuum theorem are not thereby closed.

### 3.2 Serial multiplicity master

Source: `BQG_SERIAL_MULTIPLICITY_MASTER_RESULT.md`.

The same (S_4)-twirled local Euclidean master has a finite PASS on (j=2,5/2,3,7/2,4), with nonzero lowest-channel isolation gaps and symmetry defects around (10^{-14}). The gaps are nonmonotonic. This establishes finite channel selection, not physical-projector convergence.

### 3.3 W5 scalar/Gauss-charge null checks

Sources:
- `docs/numerical_continuation/reverse_inverse_charge_checks.json`
- `docs/numerical_continuation/HE_full_sum_charge_check.json`
- `docs/numerical_continuation/HE_projection_residual.json`
- `BQG_W5_SCALAR_SYMMETRY_BLOCKER.md`

The saved reverse/inverse checks contain 25 entries with forbidden (J,M) squared norm reported as zero. The complete all-Gauss-source Euclidean sum reports forbidden Gauss norm squared (1.4116580\times10^{-30}), fraction (4.4748360\times10^{-31}), with 24 words and 98,304 branches processed. The projection-residual audit reports outside-spectator norm squared (1.71104795\times10^{-31}), but explicitly says the primitive branches were thresholded and the full operator error is not certified.

For the newer exact-key W5 failure-context calculation ((K=0), node 0, (J_{max}^{(2)}=10)), the unpruned independent diagnostic gives wrong-sector weight (4.01671023\times10^{-30}), whereas production pruning gives (5.39051261\times10^{-20}). This is a strong numerical warning about pruning/summation order. It is not yet a proof of which individual pruning stage causes the difference, nor a multi-seed production fix.

### 3.4 K-holonomy inverse identity

Source: `docs/numerical_continuation/K_all_holonomy_closure.json`.

The four-edge inverse identity residuals are approximately (1.3\times10^{-15}) to (1.6\times10^{-15}), with inverse-basis leakage (4.44\times10^{-16}). The file explicitly limits its scope to the computed (K(\phi)) and does not certify global unitarity.

## 4. What the current j=1/2 -> j=1 numerical run actually tests

Source: `scripts/bqg_jhalf_to_j1_return_rg_gate.py` and `BQG_JHALF_TO_J1_RETURN_RG.md`.

The current executable constructs the 32-dimensional logical return kernels

```math
K_f=A_f^\dagger A_f,
\qquad K_c=A_c^\dagger A_c,
\qquad A=H_{E,0}^{sine}+H_{E,1}^{sine},
```

by applying the one-hit map to the 32 embedded logical basis states and taking their Gram matrices.

It therefore measures the pullback (W_{global}^\dagger M_cW_{global}) on the selected logical carrier. It does **not** construct the full coarse-domain (M_cW_{global}) image, so it cannot measure the off-carrier leakage (L=(I-WW^\dagger)M_cW). It also does not assemble the full graph-changing master or its two-sided (E_a,F_a) residuals.

The earlier finite run found a large raw return-kernel mismatch (operator-relative (0.5232), Frobenius-relative (0.4529)). Its post-hoc Frobenius scale/shape decomposition left a (0.0719) relative shape residual after a scalar fit. That fit is not a physical renormalization and does not establish the scale-free (epsilon^{sf}) criterion.

## 5. Correct next calculation order

1. Inspect the one-pass W5 tolerance-ablation run (workflow 42820210450). It generates each branch once and evaluates projection/accumulator tolerance variants in parallel. Compare full maps, forbidden-sector weights, allowed-sector weight, and grouping/order effects against the saved zero/null controls. Do not alter the production guard before this is resolved.
2. For the (j=1/2\to1) return-kernel result, retain the current 32x32 pullback comparison as a diagnostic, but do not label it a full master refinement test.
3. Inspect the running coarse-domain carrier-leakage calculation (workflow 42797344820). It targets the same finite positive return operator and reports (L) and (D) separately without constructing an unnecessarily dense full coarse matrix. If its Hermiticity/pullback checks pass, extend the method to the full graph-changing habitats.
4. Only then proceed to the actual graph-changing constraint family: freeze (iota_0,iota_1), measure (E_a,F_a), assemble/bound (R), and compute (epsilon^{sf}) using independently computed gaps.
5. Repeat at a second refinement pair and a declared regulator variation. A single (j=1/2\to1) datum cannot establish asymptotic decay or summability.
6. Keep HDA, Lorentzian adjoint closure, graph-changing physical projector, and TT transfer extraction as separate gates. Passing one does not imply the others.

## 6. Current evidence boundary

The repository already contains many exact finite identities, symmetry-protected zeros, finite master-rank closures, and theorems linking residuals to physical projectors. These are substantive results and remain valid within their stated scopes.

The unresolved continuum claim is narrower and more precise:

```math
\epsilon_j^{sf}\to0
\quad\text{and preferably}\quad
\sum_j\epsilon_j^{sf}<\infty,
```

for the **actual graph-changing BQG master**, with a dynamically selected microscopic refinement lineage and compatible forward/return habitat embeddings.

That is the next theorem-to-data bridge. No existing finite zero should be discarded, and no finite zero should be promoted beyond its declared sector.


## 7. New implementation added after this audit

The script `scripts/bqg_jhalf_to_j1_return_carrier_leakage.py` has now been added. It is designed to compute a narrower but genuinely new datum: the carrier leakage (L) and internal mismatch (D) for the compressed 243-dimensional coarse (j=1) local-singlet return operator.

It obtains the coarse action on the embedded 32D carrier by applying the declared Hermitian sine-ordered two-node operator twice, then independently checks the pullback against the Gram matrix ((AW)^\dagger(AW)). The script reports the discrepancy rather than silently assuming that cutoff/pruning preserves Hermiticity.

Scope remains deliberately limited: the 243D product of local (j=1) singlet multiplicity spaces. Spin-changed graph-changing habitats and the full master are excluded. The numerical run is in progress in workflow 42797344820; no (L,D) values are claimed until its artifact is available and its Hermiticity/pullback checks are inspected.


### 1.4 Refinement ancestry is not necessarily the lowest local master branch

Sources: `BQG_LOWEST_CHANNEL_REFINEMENT_NOGO.md`, `BQG_MASTER_BRANCH_TRACKING.md`, and `BQG_RG_ANCESTRY_PATH.md`.

A finite counterexample already rejects the shortcut “follow the lowest local master eigenchannel at every scale.” The recorded lowest-channel overlaps oscillate strongly (including values near (0.00108) and (0.9920)); the all-channel resolution is required because the microscopic binary refinement may switch between multiplicity branches.

The correct finite ancestry diagnostic is the full overlap matrix

```math
\mathcal F_j^{ab}
=\frac12\operatorname{Tr}
(P^{block}_{j,a}P^{master}_{j+1/2,b}),
\qquad \sum_b\mathcal F_j^{ab}=1
```

for a normalized rank-two source copy, followed by a no-fit maximum-product path through consecutive scales. This selects a candidate representation-RG lineage; it is not a physical-time probability and does not replace the graph-changing habitat residuals.

The global correlated-ancilla master-selection note separately identifies a two-dimensional (S_5)-trivial ancilla sector and proposes the actual symmetric Euclidean master as the no-fit selector. Its final physical-master interpretation remains finite-control only.


## 8. Five-step full ancestry-matrix result (finite representation-RG certificate)

The full-flow workflow `42797583866` completed all five consecutive
`s_2=3,4,5,6,7` shards, covering
`j=3/2 -> 2 -> 5/2 -> 3 -> 7/2 -> 4`. The exact matrices are persisted in
`docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_{3,4,5,6,7}.json`.

The one-step matrices are:

```math
F_{3/2\to2}=(0.0540421326,\ 0.9459578674)
```

```math
F_{2\to5/2}=
\begin{pmatrix}
0.9656128797&0.0343871203\\
0.1311507029&0.8688492971
\end{pmatrix}
```

```math
F_{5/2\to3}=
\begin{pmatrix}
0.8397696212&0.1602303788\\
0.1173115629&0.8826884371
\end{pmatrix}
```

```math
F_{3\to7/2}=
\begin{pmatrix}
0.0010798725&0.9126526548&0.0862674727\\
0.0118507795&0.2045240738&0.7836251467
\end{pmatrix}
```

```math
F_{7/2\to4}=
\begin{pmatrix}
0.9920365891&0.0045659110&0.0033975000\\
0.0019135244&0.8233707897&0.1747156859\\
0.0060241961&0.1270873565&0.8668884474
\end{pmatrix}
```

Across these five finite matrices, the maximum row-stochasticity defect is
(1.11\times10^{-15}), the maximum target paired-eigenvalue spread is
(2.13\times10^{-14}), the maximum generated-pair spread is
(2.22\times10^{-16}), and Gauss-projection leakage ranges from
(1.67\times10^{-31}) to (2.89\times10^{-30}). Thus the five individual
channel resolutions pass their stated finite numerical checks.

The maximum-product dynamic-programming path is

```math
a_0^{(3/2)}\to a_1^{(2)}\to a_1^{(5/2)}
\to a_1^{(3)}\to a_2^{(7/2)}\to a_2^{(4)}.
```

Its five edge overlaps are
(0.9459578674471627, 0.8688492971205302, 0.8826884371011211,0.7836251466890334, 0.866888447374884), with product
(0.49282787562801145). The all-low-index shortcut has product
(4.694570240506031\times10^{-5}), so the maximum-product local path score
is about (1.0498\times10^4) times larger.

**Important refinement of the claim:** this is a maximum-product *local*
ancestry candidate. The product of scalar one-step matrices (F_jF_{j+1})
is equivalent to an intermediate channel measure-and-prepare map; it discards
off-diagonal channel coherence and traceless structure inside each
multiplicity block. The maximum-product path is a different statistic again:
it selects one intermediate channel rather than summing over channels. Neither
is automatically the overlap of the composed microscopic refinement.

The direct two-step comparison is specified in
`BQG_MULTISCALE_CHANNEL_COMPOSITION_GATE.md` and implemented by
`scripts/bqg_two_step_channel_composition.py`. It compares the actual
generated rank-2 carrier propagated through two adjacent microscopic
refinement tensors against (F_jF_{j+1}), while recording intermediate
off-diagonal/internal-anisotropy blocks and the maximum-product local path
score. Its workflow is pending; no factorization outcome is claimed yet.

The full five-step data and local path are also recorded in
`BQG_FULL_CHANNEL_FLOW_MATRIX.md`,
`BQG_RG_ANCESTRY_PATH.md`, and
`docs/numerical_continuation/BQG_RG_ANCESTRY_PATH_5STEP.json`.

**Claim boundary:** these are finite representation-RG projector-overlap
results. They are not physical-time probabilities, physical histories,
continuum fixed points, or a proof of physical-projector convergence. The
full graph-changing habitat residuals (E_a,F_a), off-carrier leakage, and
independently defined master gaps remain separate gates.
 

## 9. HDA/gravity closure family: what is closed and what is conditional

The wider HDA branch contains several independent results which must remain distinct:

- `ADM_HDA_PARAMETER_SELECTION.md` gives the classical local-ADM identity and finite spectral regression selecting (c=1/2) and (AB=1), leaving the familiar (G) and (Lambda) freedoms. This is not a full quantum K5 HDA proof.
- `BINARY_HDA_UNIQUENESS_CHAIN.md` links qubit full isotropy to (D=3), (c_{DW}=1/2), (AB=1), (z=1), and two tensor configuration modes, but explicitly marks the microscopic dynamical hypotheses as conditional/unproved.
- `DENSITIZED_QUANTUM_HDA_TARGET.md` freezes a preferred inverse-volume-free graph-changing quantum identity and a directional/null-channel test. The complete graph-changing dual response still has to pass.
- `GRAPH_CHANGING_HDA_TARGET.md` corrects the interpretation of fixed-sector leakage: leaving the original spin sector is not automatically a graph change or an HDA anomaly. The safe first HH column has (44.4033\%) of its norm squared in (j=0)-link sectors which reduce to a different abstract graph; the remaining (55.5967\%) stays on the same K5 graph with nonzero changed spins.
- `OFF_SHELL_HDA_HABITAT_TARGET.md` identifies the decisive Level-C test: a nontrivial off-shell vertex-smooth habitat response matching the diffeomorphism structure function. Kinematical covariance and group-averaged/on-shell vanishing alone are insufficient.
- `DUAL_K5_HDA_COCHAIN_TARGET.md` closes the exact combinatorial support: the discrete lapse bracket (N_vM_w-N_wM_v) is supported on the single dual edge ((vw)). The geometry-dependent sharp map from this cochain to a tangential deformation remains OPEN.

This branch therefore has strong algebraic and classical closure results, plus a sharply specified quantum falsification target. It does not yet have a full off-shell graph-changing quantum HDA pass.

## 10. Constraint resolvent is not yet a physical propagator

`FESHBACH_INTERBLOCK_EFFECTIVE_KERNEL.md` supplies an exact Feshbach/block-Lanczos definition of the coarse constraint resolvent. It explicitly keeps the constraint spectral parameter distinct from physical frequency. The physical history graviton propagator still requires the separate Hamiltonian-to-effective-action bridge, a specified physical state/measure, and the full interacting history kernel. This is another example where the operator definition is closed but the physical interpretation remains a separate gate.


## 11. Dimension and relational-history branches: conditional selectors versus exact finite controls

These branches contain further important distinctions:

- `BINARY_ADJOINT_DIMENSION_THEOREM.md` is a conditional structural theorem: if local geometry is the full traceless-adjoint observable space and local frame freedom is the full unitary conjugation group, then the binary carrier selects (D_{spatial}=3) and (d_{spacetime}=4). The assumptions are physical hypotheses, not empirical facts.
- `TWO_FORM_DIMENSION_PRINCIPLE.md` gives a conditional selector for (d=4) from a connection/curvature two-form pairing. Its first dynamical arrow remains open.
- `HODGE_DIMENSION_SELECTOR.md` is explicitly a falsification hypothesis: the required microscopic Hodge-duality map must be generated by a frozen rule and survive coarse graining.
- `MANIFOLD_DIMENSION_GATE.md` is a coordinate-free necessary topology test; it is not yet a completed manifold-emergence result.
- `BQG_QFI_SPECTRAL_DIMENSION_RESULT.md` provides exact hypercubic calibration and a branching-topology discriminator, but does not prove that irregular dynamical BQG graphs flow to a stable (d_s=3) fixed point.
- `Q2_RELATIONAL_HISTORY_PROJECTOR.md` resolves the clock-averaging issue exactly in a finite control: clock-only averaging kills the nontrivial phase sector, while the combined clock-system projector preserves relational correlations. The chosen (R=J) is explicitly a positive control, not the physical graph-changing evolution. The genuine physical history/rigging map remains open.

Thus the project already has several conditional dimension-selection theorems and exact finite projector controls. They are complementary, not interchangeable, and their dynamical hypotheses must not be silently promoted to proved facts.
