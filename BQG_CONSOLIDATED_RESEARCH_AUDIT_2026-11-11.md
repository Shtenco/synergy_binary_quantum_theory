# BQG consolidated research audit and result register — 2026-11-11

> **Status:** research consolidation / audit memorandum. This file records the results and code-level findings available in the research record as of 2026-11-11. It is not a claim that the full physical theory, continuum limit, or experimental predictions have been established.
>
> **Repository:** `Shtenco/synergy_binary_quantum_theory`  
> **Base commit inspected for this consolidation:** `23669209c83fe0aa2741612e5b2b119562ee6d7d` (`main`)  
> **Consolidation branch:** `research/consolidated-bqg-audit-20261111`

## 1. Executive register

The work has produced several nontrivial finite-dimensional mathematical and computational results:

1. An exact support/rank argument for the collective L1 E-sector at q4, with reported projected rank 24 and a matching 24-dimensional upper bound.
2. A strict-interior boundary Gram structure with six symmetry sectors and rank 6.
3. An S4-equivariant metric compression of the 24-dimensional L1 sector to six coarse edges, with a resolved (A_1\oplus E\oplus T_2) decomposition.
4. A local SU(2)/S4 bilinear refinement theorem: the unique equivariant bilinear map is a partial isometry on the selected product sector.
5. A global K5 ancilla orientation-covariance result, alongside an important non-uniqueness result.
6. A Peter–Weyl ancilla refinement construction that solves a local representation-theoretic obstruction after adding one binary spin, while leaving the many-edge Gauss-compatible habitat map open.
7. Exact finite-dimensional two-sided habitat and scale-free residual identities.
8. A canonical refinement residual decomposition with a useful orthogonal split into leakage and in-habitat mismatch.
9. A K5 graph/HDA diagnostic programme that has exposed important regulator, shell-compression, pruning, and gate-contract limitations.

These results should be preserved as **finite mathematical theorems, finite numerical certificates, or implementation audits according to their actual evidence class**. None alone establishes the continuum BQG dynamics or empirical validity.

## 2. Evidence discipline

Every result in this register belongs to one of these categories:

- **Exact algebraic identity/theorem:** follows from stated finite-dimensional definitions and hypotheses.
- **Finite numerical certificate:** a computed result reported by a committed script/workflow and its result document. It is not assumed to have been rerun during this consolidation.
- **Code audit finding:** follows from inspection of implementation logic; it is not a runtime test.
- **Open conjecture/target:** proposed next construction, not yet established.

No test or workflow is claimed to have been freshly executed in this consolidation. Historical workflow identifiers and reported numerical values are provenance, not fresh validation.

## 3. L1 collective E-sector: q4 rank-24 result

Primary theorem:
- [COLLECTIVE_L1_E_Q4_RANK_THEOREM.md](COLLECTIVE_L1_E_Q4_RANK_THEOREM.md)
- [worker](scripts/collective_l1_block_e_q4_worker.py)
- [collector](scripts/collective_l1_block_e_q4_collect.py)

Recorded provenance:
- Commit: `919f64856c2e2b232c94ffbd48593f1c4d0c2d6b`
- Historical workflow run: `31965359681`
- Historical shard accounting: 24 workers plus collector, 25/25 successful.

Finite setup and reported certificate:
- First barycentric subdivision of the 16-cell.
- 24 fine chamber columns, 384 fine tetrahedra, and 768 dual links.
- Physical sine Euclidean operator and exact (P_4) projection onto states with exactly four changed spin labels.
- Reported projected rank: 24.
- Reported Gram eigenvalue range: approximately 0.8436999771 to 1.8553221449.
- Reported minimum singular value: approximately 0.9185314241.
- Reported union support: 768.

The support argument is that a q4 plaquette touches four edges once. Starting from doubled spin 1, the changed labels are 1→0 or 1→2. q6/q8 terms have at least six changed plaquette edges and therefore cannot contribute to the exact four-changed-edge support. The finite rank squeeze is
[
24=\operatorname{rank}(P_4E\Omega)
\leq \operatorname{rank}(E\Omega)\leq24,
]
hence (operatorname{rank}(E\Omega)=24), provided the stated operator/support construction is the one implemented.

**Audit qualification:** the worker aggregates/prunes sequentially in `add()` at tolerance (10^{-10}), making the result potentially order-dependent. The reported Gram margin is encouraging but does not itself bound the total discarded amplitude. A rigorous stability certificate should bound the aggregate discarded vector/norm and establish a perturbation bound such as
[
\|\Delta G\|_2<\lambda_{\min}(G).
]
A tolerance sweep and order-invariance check should be added.

## 4. Strict-interior boundary Gram: rank 6

Primary theorem:
- [COLLECTIVE_L1_STRICT_INTERIOR_BOUNDARY_RANK_THEOREM.md](COLLECTIVE_L1_STRICT_INTERIOR_BOUNDARY_RANK_THEOREM.md)
- Research branch recorded in the theorem provenance: `research/l1-amplitude-package-20261110`.

Recorded finite construction:
- 24 chamber vectors.
- 20 microscopic states per chamber; 480 state occurrences.
- 192 distinct internal constrained states.
- 36 internal dual links and 24 open boundary magnetic indices.
- Double-layer tensor-network Gram construction, avoiding explicit enumeration of a (2^{24}) boundary tensor.

After parity rephasing and grouping by the unordered pair of first two permutation entries, the 24 chambers form six classes of four. The reported Gram has the structure
[
G=d\big[(1-r)\operatorname{diag}(J_4,\ldots,J_4)+rJ_{24}\big],
]
with
[
d=3.137686859428218\times10^{-10},\qquad
r=0.8523308467411363.
]
Reported relative structure defect: (4.084649355260909\times10^{-16}).

The spectrum is
[
\lambda_0=4d(1+5r),\qquad
\lambda_{\rm shape}=4d(1-r)
]
with shape eigenvalue multiplicity 5; the reported rank is 6.

**Scope:** this is a strict-interior q4/vacuum-exterior bridge. It is not yet the complete production (W_{\rm block}) including crossing sectors.

**Branch hygiene:** the `audit/l1-amplitude-package-20261110` and `research/l1-amplitude-package-20261110` branches have incompatible E-only helper signatures. The audit branch's `stable_sum_states` returns two values and its validator takes three arguments; the research branch returns three values (including discarded absolute/relative norms) and its validator takes four arguments including `shard_root`. Do not mix these implementations or transplant call signatures between them.

The research branch has stronger shard semantic reconstruction, but the E-only depth-two package remains a regression precursor, not the full E+S+R physical operator. Its backend equivalence check covers the unrefined 16-cell first/second hit and representative L1 positions 0 and 7, not all 24 positions. Its second-hit check reuses the same reference first-hit state for both backends rather than exhausting RR/RL/LR/LL combinations. The depth-two S4 residual also needs a complex-valued residual gate; taking only the real part while merely reporting imaginary components is insufficient.

## 5. L1 S4 metric compression

Primary result:
- [L1_Q4_S4_METRIC_COMPRESSION_RESULT.md](L1_Q4_S4_METRIC_COMPRESSION_RESULT.md)
- Script: [scripts/collective_l1_q4_s4_metric_compression.py](scripts/collective_l1_q4_s4_metric_compression.py)
- Script SHA recorded in the result: `6ac401c86ae8902d330002b227aca5a431cbb637`.

Reported finite result:
- Full complex 24×24 Gram matrix, not only the E-only depth-two package.
- All 24 left-regular S4 commutators reported at approximately (6.9\times10^{-16}).
- Compression map from 24 chamber coordinates to six coarse edges, with (B^\dagger B=I_6).
- Six-dimensional representation decomposes as (A_1\oplus E\oplus T_2).
- Reported eigenvalues: (A_1=0.843699977122487), (E=1.111191787558474), (T_2=1.022027850746478).
- The E/T2 splitting is approximately 0.0891639368119954 (about 8.36% relative to the T2 value).

**Implementation audit item:** the script forms `K=(B.T @ G @ B).real`, dropping the imaginary part. Exact S4 invariance and the real multiplicity-free decomposition motivate real scalar blocks, but the implementation must explicitly measure and gate
[
\frac{\|\operatorname{Im}(B^\dagger GB)\|}{\|B^\dagger GB\|}.
]
Do not silently discard a nonzero imaginary component.

## 6. Local microscopic bilinear refinement

Primary theorem:
- [BQG_MICROSCOPIC_BILINEAR_REFINEMENT_THEOREM.md](BQG_MICROSCOPIC_BILINEAR_REFINEMENT_THEOREM.md)
- Gate: [scripts/bqg_microscopic_bilinear_refinement_gate.py](scripts/bqg_microscopic_bilinear_refinement_gate.py)

Representation-theoretic result:
[
\mathcal A=\operatorname{Inv}_{SU(2)}[(1/2)^{\otimes4}]
\simeq[2,2]_{S_4},\qquad \mathcal A^{S_4}=0.
]
The tensor product decomposes as
[
[2,2]\otimes[2,2]=[4]\oplus[2,2]\oplus[1^4].
]
There is a unique (S_4)-equivariant bilinear refinement map (F), up to normalization, satisfying
[
2F^\dagger F=P_{\rm prod},\qquad
2FF^\dagger=P_{22}.
]
Thus (sqrt2F) is a partial isometry on the selected product sector, not on the entire four-dimensional product space.

**Open physical question:** what physical dynamics prepares/populates (P_{\rm prod})? The algebraic uniqueness theorem does not by itself establish dynamical selection.

## 7. Global K5 ancilla symmetry and master-selection caveat

Primary sources:
- [BQG_GLOBAL_K5_ANCILLA_SYMMETRY.md](BQG_GLOBAL_K5_ANCILLA_SYMMETRY.md)
- [scripts/bqg_global_k5_ancilla_symmetry_gate.py](scripts/bqg_global_k5_ancilla_symmetry_gate.py)
- [BQG_GLOBAL_ANCILLA_MASTER_SELECTION.md](BQG_GLOBAL_ANCILLA_MASTER_SELECTION.md)
- [scripts/bqg_global_ancilla_master_selection_gate.py](scripts/bqg_global_ancilla_master_selection_gate.py)

Recorded symmetry result:
- Logical tensor (V_5\in([2,2]_{S_4})^{\otimes5}), with S5 vertex relabeling and local S4 recoupling.
- Bare representation acts on the selected tensor with sign character:
  (R(g)|V_5\rangle=\operatorname{sgn}(g)|V_5\rangle).
- Orientation-corrected representation gives
  (R_{\rm oriented}(g)|V_5\rangle=|V_5\rangle).
- Reported maximum covariance error: (2.43\times10^{-15}).
- Both trivial and sign projectors have rank 2. Symmetry compatibility therefore does not uniquely select (V_5).

Recorded script/document SHAs:
- Symmetry document: `c353970eb623a2d0bc708f03466021b3550fec0d`.
- Symmetry script: `5ac00ee4400c445f2a689639aec7f5a8d84d0fe9`.
- Master-selection document: `6bf40471b343a5f6a281c113368aa2e9da4663e8`.
- Master-selection script: `dc3f68c13f64f9919e68a46d952d4fd7ef3b2ac9`.

**Gate-contract defect to resolve:** the symmetry gate appends candidates only when projector rank equals one and sets `passed=bool(candidates)`. If the mathematically reported projector rank is 2, the gate can fail despite the symmetry compatibility result being valid. Separate at least three statuses: (i) covariance compatibility, (ii) multiplicity/non-uniqueness, and (iii) unique selection.

**Master-selection limitation:** the master is
[
M_{\rm anc}=\sum_{v=0}^4(H_{E,v}^{\rm sine})^\dagger H_{E,v}^{\rm sine}
]
on the 32D all-(j=1/2) logical K5 carrier. The gate checks commutator compatibility, invariant weight, a gap, and overlap normalization. It does not require (V_5) to be the lowest eigenline, overlap with the lowest eigenline to be approximately one, or the minimum eigenvalue to vanish. It can therefore pass without dynamically selecting (V_5) or establishing a physical master kernel. The workflow `.github/workflows/bqg-global-ancilla-master.yml` is push-path/manual triggered; no fresh workflow result is asserted here.

A proposed global map (T_{V5}=\widetilde F^{\otimes5}(I_{\rm old}\otimes|V5\rangle)), with (widetilde F=\sqrt2F), should be treated as a target until the full 32×32 Gram, spectrum, rank, condition number, S5 covariance, and shared-edge/Gauss/orientation constraints are computed reproducibly.

**Reproducibility warning:** numerical global K5 correlated-ancilla Gram/rank/channel figures mentioned in earlier discussion are not promoted to established repository results in this register because the corresponding reproducible script/artifact was not available in the audit record. They must be independently regenerated before citation as certificates.

## 8. Peter–Weyl ancilla refinement

Primary theorem:
- [BQG_ANCILLA_PETER_WEYL_REFINEMENT_THEOREM.md](BQG_ANCILLA_PETER_WEYL_REFINEMENT_THEOREM.md)

The theorem establishes the local obstruction: there is no nonzero SU(2) intertwiner (V_j\to V_{j+1/2}) without enlarging the representation carrier. After adjoining one binary spin-(1/2), the symmetric (j+1/2) channel admits a unique SU(2)-intertwining isometry.

In the occupation basis the recorded map is
[
J_n|k\rangle_{n+1}
=\sqrt{\frac{n+1-k}{n+1}}|k\rangle_n|0\rangle
+\sqrt{\frac{k}{n+1}}|k-1\rangle_n|1\rangle.
]

**Open:** a many-edge, Gauss-compatible habitat map with consistent shared-edge ancilla bookkeeping. The local representation-theory result is not yet that global construction.

## 9. Two-sided habitat refinement and scale-free residuals

Primary sources:
- [BQG_TWO_SIDED_HABITAT_REFINEMENT_THEOREM.md](BQG_TWO_SIDED_HABITAT_REFINEMENT_THEOREM.md)
- [scripts/bqg_two_sided_habitat_refinement_gate.py](scripts/bqg_two_sided_habitat_refinement_gate.py)
- [BQG_SCALE_FREE_REFINEMENT_RESIDUAL_THEOREM.md](BQG_SCALE_FREE_REFINEMENT_RESIDUAL_THEOREM.md)
- [BQG_CANONICAL_REFINEMENT_RESIDUAL_DECOMPOSITION.md](BQG_CANONICAL_REFINEMENT_RESIDUAL_DECOMPOSITION.md)

The finite-dimensional two-sided identity is
[
E_a=C_{c,a}\iota_0-\iota_1C_{f,a},\qquad
F_a=C_{c,a}^\dagger\iota_1-\iota_0C_{f,a}^\dagger,
]
and
[
R=M_c\iota_0-\iota_0M_f
=\sum_a(C_{c,a}^\dagger E_a+F_aC_{f,a}).
]
This is an exact algebraic identity under the stated operator definitions. The current gate tests random finite matrices; it does not yet supply physical BQG (E_a,F_a) matrices.

The scale-free residual is
[
\epsilon^{sf}=\inf_{\rho>0}
\frac{\|M_c\iota-\rho\iota M_f\|}
{\min(\Delta_c,\rho\Delta_f)}.
]
It is invariant under independent positive rescalings. The existing gate uses random matrices and a coarse grid, so it is a regression check, not a continuum result or necessarily a certified global minimization.

The canonical isometry is
[
W=\begin{pmatrix}0&2/3\\1&0\\0&-\sqrt5/3\end{pmatrix},
\qquad W^\dagger W=I_2.
]
With (P=WW^\dagger), the residual decomposes as
[
R=M_cW-WM_f=L+WD,
quad L=(I-P)M_cW,
quad D=W^\dagger M_cW-M_f,
]
and
[
R^\dagger R=L^\dagger L+D^\dagger D.
]
This separates leakage out of the refined habitat from mismatch inside the habitat. The current gate uses random Hermitian matrices, not actual BQG master data.

## 10. K5 HDA: reachability, finite diagnostics, and regulator wall

Implementation-level follow-up: [K5_SAFE_HDA_COLUMN_IMPLEMENTATION_AUDIT_2026-11-11.md](K5_SAFE_HDA_COLUMN_IMPLEMENTATION_AUDIT_2026-11-11.md) records the threshold cascade, adjointness and volume checks, cutoff scope, and acceptance criteria for the safe Peter–Weyl column.

Primary sources:
- [K5_ORIENTED_QUANTUM_HDA_RESULT.md](K5_ORIENTED_QUANTUM_HDA_RESULT.md)
- [DUAL_K5_HDA_COCHAIN_TARGET.md](DUAL_K5_HDA_COCHAIN_TARGET.md)
- [scripts/k5_hda_reachable_basis_gate.py](scripts/k5_hda_reachable_basis_gate.py)
- [scripts/k5_quantum_hda_span_gate.py](scripts/k5_quantum_hda_span_gate.py)
- [scripts/bqg_k5_graph_sector_master_shell.py](scripts/bqg_k5_graph_sector_master_shell.py)
- [scripts/k5_oriented_quantum_hda_gate.py](scripts/k5_oriented_quantum_hda_gate.py)
- [scripts/k5_thiemann_constraint_gate.py](scripts/k5_thiemann_constraint_gate.py)

### 10.1 Reachability accounting

The combinatorial support enumerator starts from all (j=1/2) on K5. Its recorded local move touches triangle links once and the radial link twice. It reports:
- after one (H): 120 Gauss-admissible spin assignments, dimension upper bound 816;
- after an (HH) pair: 4193 assignments, dimension upper bound 24364;
- maximum spin (j=5/2) after two hits.

This is combinatorial reachability and a dimension upper bound; it does not include amplitude cancellations. It is not a proof that all these states have nonzero amplitudes.

### 10.2 Fixed-sector HDA diagnostic

The quantum HDA span gate uses the old 140D Gauss Hilbert space and 32D fully-active sector, and examines raw/Hermitian/anti-Hermitian node orderings. Its HDA pass criterion checks only a median full-sector leakage below 1% for the Hermitian completion. It does not test the (1/(3V)) coefficient or the actual structure function (D(k,l)).

The oriented quantum HDA result reports:
- orientation signs restore S5 covariance across all ten pairs;
- Hermitian completion commutator action norm: 3906.6220702801;
- graph leakage: (sqrt{37/69}=0.7322785563);
- fixed-sector rank 12, with (32/69) norm-squared inside and (37/69) outside;
- (Q_{ij}=-iP[H_i,H_j]P) spectrum (-768) (multiplicity 6), zero (multiplicity 20), (+768) (multiplicity 6);
- the projected operator is a pure weight-three Pauli operator, with weight-(leq2) residual 1 and weight-(leq3) residual below (1.3\times10^{-14});
- the projected SO(5) Lie skeleton is exact, while the full operator bracket leakage is (sqrt{296/297}=0.998315).

Interpretation: this is a strong finite symmetry/algebra diagnostic and a clear failure of fixed-sector closure for the full bracket. It is not a successful physical HDA closure result.

### 10.3 Regulator wall

The vector-5 (J_{\max}=1/2) diagnostic is not a safe regulator for repeated (HH) action. The documented safety target is
[
J_{\max}\geq j_{\rm in}+2,
]
so for (j_{\rm in}=1/2), use at least (J_{\max}=5/2) for the two-hit diagnostic. Higher shell depths require a new cutoff bound; do not reuse (5/2) automatically.

The old K5 Thiemann-like finite constraint has a 140D Gauss space, 32D active sector, 60 local kernels, and common kernel (V_5). Its scope is a finite Euclidean/canonical regularization, not Lorentzian dynamics, HDA closure, or continuum physics.

### 10.4 Dual cochain target

For lapse functions (N,M), the cochain
[
\omega_{vw}=N_vM_w-N_wM_v
]
has a basis-lapse-pair support property: the pair ((\delta_i,\delta_j)) selects dual edge ((ij)). The geometry-dependent musical map (\sharp_{E,q}) remains open.

## 11. Critical graph-sector shell audit: exact leakage correction

Script:
- [scripts/bqg_k5_graph_sector_master_shell.py](scripts/bqg_k5_graph_sector_master_shell.py)

The current shell code forms compressed node operators (C_v=P_SH_vP_S) and a compressed master
[
M_S=\sum_v C_v^\dagger C_v.
]
It also records only a scalar total outside-shell norm per vertex, not the full leakage operator.

Let (P_S) project onto the finite shell and (Q_S=I-P_S). Define
[
C_v=P_SH_vP_S,\qquad L_v=Q_SH_vP_S.
]
If the inside/outside decomposition is orthogonal in the actual Hilbert-space inner product and the full action (H_v P_S) is represented without truncation, then the exact pullback identity is
[
P_SH_v^\dagger H_vP_S
=C_v^\dagger C_v+L_v^\dagger L_v.
]
Therefore
[
M_{\mathrm{pullback},S}
=\sum_v(C_v^\dagger C_v+L_v^\dagger L_v)
=M_S+\sum_vL_v^\dagger L_v.
]
Consequently,
[
\ker M_{\mathrm{pullback},S}\subseteq\ker M_S.
]
The compressed-only master can have spurious zero modes because it omits positive semidefinite leakage contributions.

**Important precision:** the leakage identity itself does not require (H_v) to be self-adjoint. It follows from the orthogonal decomposition (H_v P_S=P_S H_v P_S+(I-P_S)H_v P_S), and is an identity for (P_S H_v^dagger H_v P_S). Self-adjointness remains a separate requirement if the intended constraint is specifically a symmetric/Hermitian Hamiltonian. The identity does require the correct Hilbert inner product (orthogonal projection, or the corresponding Gram-weighted formulation), complete representation of graph-changing outputs, and no unaccounted pruning/truncation.

### Required implementation correction

For every input shell basis vector (e_j) and node (v):
1. Compute the full sparse output (H_ve_j), without silently discarding small amplitudes.
2. Split it into inside coefficients (C_v[:,j]) and outside coefficients (L_v[:,j]).
3. Key outside amplitudes by the full canonical spin-network basis key.
4. Accumulate the coherent Gram
   [
   (G_v^{\rm leak})_{jk}
   =\langle L_ve_j,L_ve_k\rangle
   =\sum_{s\notin S}\overline{L_v(s,j)}L_v(s,k).
   ]
   **Do not replace this with only per-column squared norms.** Off-diagonal terms encode interference between distinct shell columns that reach the same outside basis state.
5. Build
   [
   M_{\rm corrected}=\sum_v(C_v^\dagger C_v+G_v^{\rm leak}).
   ]
6. Independently check (|C_v-C_v^\dagger|), and check full-action adjointness before compression.
7. Report mode-resolved leakage: for eigenvector (u_i), (ell_i^2=u_i^\dagger(\sum_vG_v^{\rm leak})u_i).
8. Repeat under shell-depth, (J_{\max}), pruning tolerance, and deterministic ordering sweeps.

### Additional code defects identified

- **Tautological Hermiticity diagnostic:** `master_matrix` symmetrizes the matrix before `run` computes `herm_err=norm(M-M†)`. That reported error cannot diagnose pre-symmetrization non-Hermiticity. Check every (C_v), the raw master before symmetrization, and the unsymmetrized full operator action.
- **Leakage under-reporting:** a scalar sum of squared outside amplitudes per node does not construct (L_v^\dagger L_v), and cannot reveal coherent off-diagonal leakage.
- **Scale-invariance defect:** the projector routine uses `scale=max(max(abs(ev)),1.0)`. For a matrix whose natural norm is much smaller than one, this introduces an absolute floor and can classify eigenvalues differently under positive rescaling. Use a threshold relative to the actual matrix scale, with a separately justified floating-point/backward-error floor and explicit scale sweep.
- **Pruning hides support:** shell generation calls `K5.prune_state(...,1e-10)`. This can remove weak transitions before reachability/leakage accounting. Aggregate deterministically first, then record discarded absolute and relative norm budgets.
- **Weak structural pass contract:** the current `passed_structural_checks` checks small post-symmetrization Hermiticity error, minimum eigenvalue above a loose bound, and at least one group. It is not a physical kernel, gap, or convergence certificate.

### Correct comparison programme

Compute and compare:
- (M_A=\sum_v C_v^\dagger C_v), current compressed-only operator;
- (M_B=\sum_v(C_v^\dagger C_v+L_v^\dagger L_v)), exact full-master pullback to the current shell, subject to complete action and the correct orthogonal/Gram-weighted decomposition;
- (M_C), the newly compressed operator on an expanded shell.

Compare nullity, positive gap, projectors, mode-resolved leakage, shell-depth convergence, and cross-graph blocks. (M_B) and (M_C) act on different spaces and need not have identical spectra; the convergence comparison must use explicit embeddings/restrictions and controlled residuals.

The physical projector is not yet established by a finite shell calculation. It requires shell/cutoff/pruning convergence, stable kernel and gap, graph transitions, and the correct physical inner product.

## 11A. Coherent leakage-Gram implementation kernel

New utility:
- [BQG_SPARSE_LEAKAGE_GRAM_DESIGN_AND_TESTS_2026-11-11.md](BQG_SPARSE_LEAKAGE_GRAM_DESIGN_AND_TESTS_2026-11-11.md)
- [scripts/bqg_sparse_leakage_gram.py](scripts/bqg_sparse_leakage_gram.py)

The helper assembles the full coherent outside-shell Gram by matching canonical output keys across input columns and adding complex outer products. It also provides a scale-relative PSD/nullspace classifier without a unit absolute scale floor. Its input contract requires complete, coherently aggregated sparse columns; it performs no pruning itself.

The embedded tests cover a shared outside key with a nonzero imaginary off-diagonal, agreement with a dense full-output Gram on a toy problem, positive-rescaling invariance, and rejection of significant negative eigenvalues.

**Execution status:** the committed repository script was not run in this turn. A local scratch prototype of the same core algebra passed, but the committed file still requires direct execution in the target repository environment. The helper is not yet integrated into the K5 producer, whose current nested thresholds remain a blocker.

## 12. Ordered next research programme

Priority 0 — establish trustworthy operator numerics:
1. Inspect and verify the actual (H_v) adjoint construction in `k5_peter_weyl_safe_hda_column.py`, especially `apply_H_cached_state`, `apply_T`, orientation signs, local volume, and basis normalization.
2. Implement coherent sparse leakage Gram matrices and raw adjointness tests.
3. Remove absolute scale floor from zero-mode classification; add positive-rescaling invariance tests.
4. Make pruning budgets explicit and run deterministic ordering/tolerance sweeps.
5. Split the K5 gates into compatibility, finite-sector closure, regulator safety, and physical convergence statuses.

Priority 1 — resolve the habitat map:
1. Construct the one-hit graph-changed habitat, using the 816-dimensional one-hit upper-bound sector only as a starting combinatorial estimate.
2. Build edgewise binary ancilla embeddings with shared-edge consistency, Gauss projection, orientations, spins, and intertwiner multiplicities.
3. Verify (\iota_1^\dagger\iota_1=I) in the actual physical inner product.
4. Instantiate the exact two-sided (E_a,F_a) identities using real BQG operators, not random matrices.
5. Test the scale-free residual and the (L,D) canonical decomposition on actual master data.

Priority 2 — complete the collective L1 bridge:
1. Reconcile the audit/research E-only branch interfaces without mixing incompatible signatures.
2. Extend backend equivalence across all 24 L1 positions and all second-hit orderings.
3. Use complex S4 residuals and gate the imaginary component.
4. Bound all discarded amplitudes and prove Gram stability.
5. Add S/R sectors and crossing sectors to reach the full production (W_{\rm block}).

Priority 3 — global ancilla selection:
1. Separate S5 covariance from uniqueness.
2. Require the master gate to test the actual lowest eigenline, target overlap, kernel condition, and gap.
3. Compute and publish the full global refinement Gram and its singular spectrum.
4. Verify edgewise blocking exactly once per shared edge, and test Gauss/orientation compatibility.

Priority 4 — physical HDA:
1. Raise (J_{\max}) safely for each shell depth and derive cutoff growth rather than assuming it.
2. Measure full graph-changing bracket action, not only fixed-sector projection.
3. Test the expected volume coefficient and the actual (D(k,l)) structure function.
4. Establish convergence across shell depth, regulator, graph, and inner product before making physical-closure claims.

## 13. Bottom line

The strongest current achievements are precise finite-dimensional results: L1 q4 rank 24; strict-interior boundary rank 6; S4 metric compression and representation splitting; local bilinear refinement; K5 orientation-covariance; and exact finite-dimensional refinement-residual identities.

The most important newly clarified blocker is equally precise: **the K5 compressed master currently omits the coherent positive leakage Gram**. Until that is included and the underlying node operators are verified self-adjoint in the actual graph-changing Hilbert space, its compressed kernel cannot be treated as the physical master kernel.

This is a productive research position: the finite algebraic structures are becoming clearer, while the remaining work is now stated as concrete operator, convergence, and reproducibility obligations rather than vague requests for “more testing”.


## 11. K5 pruning and numerical-integrity audit — 2026-11-11

New source-level audit:
- [K5_PRUNING_AND_NUMERICAL_INTEGRITY_AUDIT_2026-11-11.md](K5_PRUNING_AND_NUMERICAL_INTEGRITY_AUDIT_2026-11-11.md)
- Inspected implementation: `scripts/k5_peter_weyl_safe_hda_column.py`, blob SHA `e7e1e0705c4b833dcee167a2036ef5f8929a4388`.

The exact source inspection identified a multi-stage hard-coded pruning cascade: (10^{-13}) in dictionary accumulation, (10^{-12}) in local channel/path coefficients, (10^{-11}) at primitive output, (10^{-10}) in cached T outputs, and (10^{-8}) in state/HH composition and top-level H regression. The graph-shell wrapper adds another (10^{-10}) threshold. These are distinct truncation layers; the current HH result must therefore be described as a legacy-thresholded finite regression, not a pruning-certified complete operator action.

The audit derives the valid telescoping error bound
[
\|x_n-\widetilde x_n\|_2
\leq\sum_s\|A_{n-1}\cdots A_{s+1}\|_2\,\|\delta_s\|_2,
]
where each (delta_s) is the coherently aggregated discarded vector at that stage. It explicitly rejects a naive sum of thresholds or dropped-key counts as a global error bound.

Additional code-level findings:
- `volume123_matrix` symmetrizes raw (Q) before measuring its anti-Hermitian defect, concealing a needed diagnostic.
- The reverse-sequence adjoint construction is structurally plausible but lacks independent bilinear action tests in the inspected file.
- Hard-coded regression constants test regression stability, not pruning completeness, adjointness, volume integrity, or HDA closure.
- Scalar outside-shell norms cannot replace the coherent off-diagonal leakage Gram (L_v^\dagger L_v).

The prescribed next milestone is a policy threaded through every operator layer (`legacy`, `no_prune`, `aggregate_then_prune`), deterministic coherent aggregation, a tiny dense oracle, bilinear adjoint tests, raw volume integrity checks, and only then a recomputed shell with full leakage Gram.

**Execution status:** this is a source audit and mathematical design record. No fresh execution of the expensive K5 HH calculation is claimed.


## 12. Deterministic coherent sparse-pruning primitive — 2026-11-11

New implementation and design record:
- [BQG_DETERMINISTIC_SPARSE_PRUNING_POLICY_2026-11-11.md](BQG_DETERMINISTIC_SPARSE_PRUNING_POLICY_2026-11-11.md)
- [scripts/bqg_sparse_pruning.py](scripts/bqg_sparse_pruning.py)
- [scripts/test_bqg_sparse_pruning.py](scripts/test_bqg_sparse_pruning.py)

The new independent helper separates coherent aggregation from pruning:
- `coherent_aggregate` groups all complex contributions by canonical key, deterministically sorts them, and uses compensated `math.fsum` separately for real and imaginary parts.
- `apply_pruning` only accepts an already coherently aggregated mapping and supports `no_prune` and `aggregate_then_prune`.
- The ledger records stage, threshold, key counts, maximum discarded amplitude, discarded L1/L2 norms, and a SHA-256 digest of the sorted discarded vector.
- A verifier checks the full-minus-pruned identity against the ledger.
- A telescoping error-bound helper implements the propagated stage-error inequality.

Twelve committed unit tests cover permutation invariance, cancellation before thresholding, constructive accumulation, a dense-vector oracle, ledger identity, no-prune retention, telescoping-bound verification, and rejection of pre-aggregation pruning.

**Verification boundary:** both source files were fetched back from GitHub and inspected. A separate local Python harness passed core checks for deterministic aggregation, cancellation/constructive accumulation, ledger norms, dense-vector comparison, and the telescoping bound. The committed unittest file itself has not yet been executed against a local checkout. Crucially, the K5 producer is not yet wired to this helper; existing HH outputs remain legacy-thresholded and are not promoted to pruning-certified status. The integration requirements and acceptance gate are recorded in the new design document.


## 13. K5 no-prune integration contract — 2026-11-11

A new source-grounded integration contract is preserved in [K5_NO_PRUNE_INTEGRATION_CONTRACT_AND_TEST_PLAN_2026-11-11.md](K5_NO_PRUNE_INTEGRATION_CONTRACT_AND_TEST_PLAN_2026-11-11.md). It enumerates every threshold in the inspected K5 producer, requires policy identity in all dependent caches, separates legacy/no-prune/aggregate-then-prune semantics, and defines five acceptance gates: primitive coherent aggregation, actual K5 tiny dense oracle, bilinear adjointness, raw volume integrity, and threshold-ladder/error propagation. It also connects the no-prune prerequisite to the full coherent graph-shell leakage Gram.

**Newly explicit source boundary:** the K5 safe-column producer is not automatically the active V10 KKV production path discussed in the two commits currently on `main`. The V10 L1/L2/L3 linearity audit is methodologically relevant, but its call-site findings are not attributed to the K5 file without a confirmed caller. At this snapshot, `main` contains two commits absent from this research branch (`8d55ff4` and `254f19f`); this branch remains isolated and no merge or history rewrite was performed.

**Verification status:** this update adds a source-audit/design contract only. It does not claim policy integration, direct execution of the committed unit tests, a K5 dense-oracle PASS, adjoint certification, raw-volume integrity, a no-prune HH run, or physical HDA closure. The existing HH result remains a legacy-thresholded finite regression.


## 14. Reproducible CI entry point for sparse-pruning helper — 2026-11-11

Added [`.github/workflows/bqg-sparse-pruning.yml`](.github/workflows/bqg-sparse-pruning.yml). It installs NumPy and runs the committed test file directly as `python scripts/test_bqg_sparse_pruning.py` (the direct invocation is intentional because that test currently imports `bqg_sparse_pruning` from the script directory). The workflow is scoped to the helper, its tests, and itself, and supports manual dispatch.

**Status boundary:** workflow definition is committed and fetched back from GitHub; no CI run result is claimed in this audit update. A green CI result must be confirmed from the actual run and job logs before the eight tests are called CI-verified. This workflow validates only the standalone pruning primitive, not the K5 producer integration or any physical BQG claim.


## 15. Sparse-pruning helper input-integrity hardening — 2026-11-11

The standalone helper received additional fail-closed validation: non-finite aggregated amplitudes are rejected by `apply_pruning`; non-finite scales are rejected even for empty mappings; ledger verification rejects negative/non-finite `atol`; and non-finite overflow in the full-minus-pruned difference is rejected. Four negative-path regression tests were added, bringing the committed suite to twelve tests. The helper and test files were fetched back from GitHub after the updates.

The GitHub Actions workflow is configured to run the committed test file directly. The repository's combined commit status query for the latest test update returned no statuses at the time of this record; therefore no CI PASS is claimed. These tests still cover only the standalone helper, not K5 producer integration.


## 16. Initial K5 producer no-prune integration — 2026-11-11

The K5 producer source scripts/k5_peter_weyl_safe_hda_column.py has been updated with explicit --mode legacy|no_prune. The default remains legacy for historical regression continuity. The strict mode bypasses the identified local coefficient, primitive output, cached T, state/HH composition, top-level H/commutator, and V5-input amplitude cutoffs; coherent aggregation is used at the main branch-combination boundaries. mode is part of the T_cached key. A MemoryError produces an explicit resource_limit result and nonzero exit rather than silent fallback. The volume routine records the raw pre-symmetrization anti-Hermitian defect by encountered spin triple.

Added six focused K5 wiring tests in scripts/test_k5_sparse_policy_wiring.py and extended the GitHub Actions workflow to install SymPy and run both the 12 sparse-pruning tests and the K5 wiring tests. Workflow paths include the producer and both test files for push and pull-request events.

Evidence boundary: source files and workflow were fetched back from GitHub after updates. No CI PASS is asserted yet; no full no-prune HH run, K5 dense oracle, production bilinear adjoint test, threshold ladder, or certified leakage Gram is claimed. Raw volume defect is currently recorded but not yet used as a fail gate. This is initial producer integration, not completion of the K5 acceptance plan.


## 17. Exact K5 integration snapshot

Latest producer blob: a87cc00c3be07b4a66c249c16ec840cb8cb5629d. Source inspection finds all previously catalogued amplitude thresholds either inside the explicit legacy branch or after a no-prune early return. The no-prune mode now coherently aggregates at primitive, T-sequence, input-state, node-H, HH-composition, commutator and V5-superposition boundaries; cache identity includes mode; MemoryError returns resource_limit with no fallback. The raw pre-symmetrization volume defect is recorded by spin triple.

The committed verification suite consists of 12 sparse-pruning helper tests and 6 K5 wiring tests. Latest test blobs: a1c2cbbb772b814355a520194e4bb32b19aa2383 and 476255880c123289b7f35d13b4bc81d03f1a4f50. Workflow blob: ffd43126b2d6003a7e551a5f51e2d5d59bde60de. The current GitHub status query contains no status entries, so CI execution remains unverified. No full no-prune HH calculation or K5 dense-oracle / adjoint certification is claimed.


## 18. Expanded K5 policy-path tests

The K5 producer wiring test suite has grown from six to nine tests. New controlled tests verify retention of a 1e-13 local primitive coefficient under no_prune (and its removal under legacy), coherent aggregation of two HH input-state contributions before the legacy 1e-8 state cutoff, and distinct T_cached outputs for identical operator arguments under different policy modes. Together with the twelve standalone helper tests, the committed suite now contains 21 tests.

These are focused software-wiring tests, not full K5 numerical certification. Latest K5 test blob: b83733713aca6bcfa03fdf9212ec7b7f97613e5c. The latest combined commit status returned no status entries; CI PASS is not asserted. Attempting a local clone to run the tests failed because the environment could not resolve github.com. Therefore the new tests have been committed and reviewed from GitHub, but not executed in this session.


## 19. K5 acceptance semantics corrected and tested

A control-flow audit found and fixed an invalid acceptance rule: the historical regression constants had been applied to no-prune output, even though the two modes intentionally implement different numerical conventions. This could make a successfully completed no-prune computation exit nonzero merely because it differed from the legacy-thresholded reference.

A dedicated apply_acceptance_gate now applies frozen regression targets only to legacy. For no_prune, it reports calculation_completed=true, passed_regression=null, regression_errors=null, and an explicit not-physically-certified acceptance scope. CLI exit code zero in no_prune mode now means execution completed; it does not mean physical validation passed. Two tests cover the distinction.

K5 wiring tests: 11. Standalone sparse-pruning tests: 12. Combined committed tests: 23. Latest producer blob: 7db5bae2bb22535cc99f9d86f5dda0f024db7fbc; latest K5 test blob: 0d43625f3d82f079f3f70c8b547fb4ac88039b8b. Files were fetched back from GitHub. CI status remains empty and no local test run was possible in this session; no test PASS is claimed.


## 20. Synthetic dense oracle and complex adjoint reference harness

Added a deterministic 5-dimensional synthetic complex sparse operator test in scripts/test_bqg_sparse_dense_oracle.py. The dense reference is assembled independently of the sparse helper. The four tests cover sparse/dense matrix equality, sparse/dense matrix-vector action, the complex adjoint inner-product identity, and tiny-amplitude retention versus exact coherent cancellation.

This is a meaningful independent algebraic control, but it is deliberately not represented as a production K5 dense oracle or as verification of the actual volume/primitive adjoints. Those production checks remain outstanding. The workflow now invokes the new test file. Committed suite inventory: 12 sparse-helper tests + 11 K5 policy-wiring tests + 4 synthetic dense-oracle/adjoint tests = 27 test methods. No local or CI execution has been confirmed; no PASS is claimed.


## 21. Production raw K5 volume operator gate

Exposed the unsymmetrized raw triple-product operator as volume123_raw_operator and refactored volume123_matrix to consume it, preserving the prior symmetrization and sqrt(abs(Q)) behavior. Added three small-sector tests on four doubled-spin triples: raw Hermiticity before symmetrization, Hermiticity/positive semidefiniteness of the derived volume matrix, and consistency of the pre-symmetrization diagnostic ledger. This directly tests the production spin/volume implementation in finite sectors, but not the full HH column or continuum physics.

The workflow now runs the new volume test file. Total committed test methods: 30 across four scripts. Execution remains unconfirmed; no PASS claim.


## 22. Local production link-hit adjoint gate

Added a finite-sector test of actual production CG link-hit matrices. For all four binary endpoint channels in the doubled-spin transition 1 -> 2, reversing edge orientation and swapping endpoint labels reproduces the conjugate-transpose matrices and preserves the real normalization. Path and sequence adjunction are also tested as involutions.

This verifies the local link-hit convention only. It does not yet prove full primitive/T/H adjointness in the coupled K5 state space. The committed suite now contains 33 test methods across five scripts; CI execution is not confirmed.


A final workflow configuration audit found and corrected incomplete pull_request path filters. All five test scripts are now covered by both push and pull_request filters and each has an explicit execution step. Workflow blob SHA: 7caec171b996b3c662b4445681764ff24f1d2a1f. No workflow result is claimed.


## 14. Production K5 primitive matrix-level adjoint gate — 2026-11-11

The K5 verification ladder now advances beyond local link-hit identities. New file: [`scripts/test_k5_primitive_dense_adjoint.py`](scripts/test_k5_primitive_dense_adjoint.py).

It assembles the actual production `apply_primitive_oriented` map on a finite, explicitly enumerated recoupling basis for the closed triangle `(0,1,2,0)), then separately assembles the reversed production map. The comparison is matrix-level:
[
B_{\mathrm{reverse}}=A_{\mathrm{forward}}^\dagger.
]
The four binary endpoint-channel pairs are covered. The spin wall is doubled-spin `Jmax2=2`; the input has all ten links at doubled spin 1, while the output basis includes all eight assignments of doubled spins 0/2 on the three triangle edges. This closed path modifies two incident links at each affected node and preserves the parity required for invariant intertwiner channels.

This is a materially stronger check than testing `adjoint_path` as an involution or testing the individual CG matrices: it exercises the actual production primitive, its branch expansion, node tensor transformations, recoupling projection and coherent sparse output. It remains one finite-sector primitive test, not a proof for all paths or the assembled node Hamiltonian.

**Test inventory:** 34 committed test methods across six scripts (12 sparse-helper, 11 K5 policy-wiring, 4 synthetic dense-oracle, 3 raw-volume, 3 local hit-adjoint, 1 production primitive dense-adjoint). The workflow explicitly invokes all six scripts and includes the new test in both push and pull-request path filters.

**Verification status:** the new source, workflow and audit record are committed and fetched back from GitHub. No test execution or CI PASS is confirmed. The next acceptance step is to run this test and inspect any failure at the first differing matrix element; do not relax tolerances or symmetrize the matrices to force agreement. After it passes, extend the same independent matrix construction to the remaining primitive path families, individual T terms, and the assembled node-H map, with explicit finite-wall boundary accounting.
