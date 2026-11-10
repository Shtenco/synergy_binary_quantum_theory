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

1. Finish and inspect the W5 tolerance-ablation run. Determine separately the effect of projection tolerance, accumulator tolerance, grouping, and summation order. Validate against the already saved zero/null controls. Do not alter the production guard before this is resolved.
2. For the (j=1/2\to1) return-kernel result, retain the current 32x32 pullback comparison as a diagnostic, but do not label it a full master refinement test.
3. Implement a coarse-domain carrier-leakage calculation for the same finite positive return operator, or derive a sparse adjoint/Gram route that obtains (M_cW_{global}) without constructing an unnecessarily dense full coarse matrix. Report (L) and (D) separately.
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

Scope remains deliberately limited: the 243D product of local (j=1) singlet multiplicity spaces. Spin-changed graph-changing habitats and the full master are excluded. The script has been committed but its numerical run has not yet been executed; no (L,D) values are claimed here.


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
