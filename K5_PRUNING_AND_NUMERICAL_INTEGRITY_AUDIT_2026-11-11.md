# K5 safe-HDA pruning and numerical-integrity audit — 2026-11-11

> **Evidence class:** source-code audit of the exact branch file fetched from GitHub. No fresh execution of the expensive SymPy/CG computation is claimed here.

- Branch: `research/consolidated-bqg-audit-20261111`
- Source: [`scripts/k5_peter_weyl_safe_hda_column.py`](https://github.com/Shtenco/synergy_binary_quantum_theory/blob/research/consolidated-bqg-audit-20261111/scripts/k5_peter_weyl_safe_hda_column.py)
- Inspected blob SHA: `e7e1e0705c4b833dcee167a2036ef5f8929a4388`
- Source length: 389 lines.

## 1. Concrete pruning cascade

The inspected implementation contains multiple hard-coded truncation boundaries:

| Stage | Lines | Current behavior |
|---|---:|---|
| `add_dict` | 191–195 | after each individual update, retains only if `abs(v)>1e-13`; partial sums can be deleted and later recreated |
| Local intertwiner extraction | 254 | discards channel coefficient unless `abs(c)>1e-12` |
| Primitive path aggregation | 261–264 | screens path products at (10^{-12}), then final output keys at (10^{-11}) |
| Cached T action | 277–282 | aggregates input-state contributions, then removes outputs at (10^{-10}) |
| Generic state prune | 292–293 | default state threshold (10^{-8}) |
| HH composition | 326–331 | calls generic prune at (10^{-8}) |
| Top-level H regression | 337–338 | prunes H outputs again at (10^{-8}) |

The graph-shell consumer `scripts/bqg_k5_graph_sector_master_shell.py` wraps `apply_H_cached_state` in another (10^{-10}) prune. Its scalar outside-shell norm is therefore computed from thresholded amplitudes.

**Conclusion:** the reported HH column is a legacy-thresholded finite calculation. It is not yet a pruning-certified evaluation of the complete finite operator.

## 2. Why thresholding must follow coherent aggregation

For each canonical output key (k), the actual amplitude is
[
a_k=\sum_{h\to k}a_{h,k}.
]
Deleting individual path terms before all histories reaching (k) have been coherently summed is not equivalent to pruning the final (a_k). Small terms may add constructively; larger terms may cancel. The repeated deletion/recreation in `add_dict` can make the result depend on contribution order.

A correct policy therefore has one rule: **aggregate all contributions to a canonical key first; only then apply a declared threshold at a named stage.**

## 3. Mathematically valid cumulative error budget

For an exact linear pipeline (x_{s+1}=A_sx_s), let a numerical stage discard the coherently aggregated vector (delta_s):
[
\widetilde x_{s+1}=A_s\widetilde x_s-\delta_s.
]
Then
[
\|x_n-\widetilde x_n\|_2
\leq\sum_{s=0}^{n-1}
\|A_{n-1}\cdots A_{s+1}\|_2\,\|\delta_s\|_2.
]
This follows by telescoping the propagated stage errors. It is useful only if the downstream operator-norm factors are bounded or computed. A raw sum of threshold values, dropped-key counts, or local norms is not a certified global bound.

For branch-local deletion before coherent aggregation, (delta_s) must represent the complete coherent omitted contribution per output key. Otherwise this bound does not describe the actual truncation error.

## 4. Required policy architecture

Thread an explicit policy through every primitive, T-cache, H-cache, composition, commutator and shell wrapper:

- `legacy`: preserve exact historical thresholds for regression reproduction only.
- `no_prune`: never delete amplitudes; fail loudly on resource exhaustion.
- `aggregate_then_prune`: coherently aggregate by canonical key, then prune once at a declared boundary and record the removed vector.

Each pruning-ledger entry must record stage, threshold, pre/post key counts, maximum removed amplitude, discarded (\ell_1) and (\ell_2) norms, sorted-key digest, and whether aggregation was complete. Keep the full discarded vector for small reference runs.

Tests required:
1. disabled-pruning aggregation invariant under deterministic permutations;
2. a cancellation/constructive-addition example proving early pruning differs;
3. sparse-vs-dense oracle on tiny spaces;
4. ledger identity: (delta=output_{unpruned}-output_{pruned}) key by key;
5. telescoping-bound verification on small explicit linear maps;
6. threshold ladder for output vectors, Gram spectra, singular values, null projectors and ranks.

## 5. Other numerical-integrity findings

### 5.1 Volume Hermiticity is hidden by early symmetrization

At line 170, `volume123_matrix` replaces raw (Q) by ((Q+Q^\dagger)/2) before diagonalization. Thus the code cannot subsequently report the Hermiticity defect of the assembled operator. Preserve (Q_{raw}), measure
[
\epsilon_H=
\frac{\|Q_{raw}-Q_{raw}^\dagger\|_F}
{\max(\|Q_{raw}\|_F,\epsilon_{machine})},
]
then apply a declared symmetrization policy. Also report eigensolver and square-root reconstruction residuals.

### 5.2 Adjoint construction is plausible but not bilinearly certified

Reversing the operation sequence and path endpoints is structurally plausible for the formal adjoint. The source does not itself establish
[
\langle x,Ay\rangle=\langle A^\dagger x,y\rangle
]
through independent action-level tests over a finite basis. Structural plausibility must not be reported as tested adjointness.

### 5.3 Regression pass is narrower than numerical-integrity acceptance

Hard-coded regression targets are useful for detecting accidental changes. They do not validate pruning stability, adjointness, raw volume Hermiticity, or the discrete HDA structure function. Keep regression status separate from mathematical/numerical acceptance status.

## 6. Graph-shell consequence: scalar leakage is insufficient

For (C_v=P_SH_vP_S) and (L_v=(I-P_S)H_vP_S), the exact identity is
[
P_SH_v^\dagger H_vP_S=C_v^\dagger C_v+L_v^\dagger L_v.
]
The shell master built only from (C_v^\dagger C_v) omits the positive leakage Gram. Moreover, if the (L_v) columns were produced by thresholded calls, the leakage Gram is itself incomplete.

The required matrix is
[
(G_{leak})_{ij}
=\sum_{\alpha\notin S}
\overline{(L_v)_{\alpha i}}(L_v)_{\alpha j}.
]
A scalar outside-shell norm/fraction cannot reconstruct off-diagonal interference between input columns sharing an outside output state.

Two separate fixes are mandatory:
1. coherent leakage Gram, not only scalar outside-shell norms;
2. complete/no-prune shell columns or a valid propagated error bound.

## 7. Priority order

- **P0:** label existing outputs as legacy-thresholded finite diagnostics; do not promote ranks/gaps as pruning-certified.
- **P1:** introduce policy-threaded coherent sparse accumulation.
- **P2:** add tiny dense oracle and bilinear adjoint tests.
- **P3:** measure raw volume anti-Hermitian defect before symmetrization.
- **P4:** recompute safe HH column with no-prune or a full threshold/error-budget ladder.
- **P5:** rebuild shell using complete columns and full (G_{leak}=L^\dagger L); compare compressed, pullback and enlarged-shell masters.

## 8. Final assessment

The source is a substantive finite Peter–Weyl K5 implementation, but its current numerical output is a **legacy-thresholded finite HH-column regression**, not yet a pruning-certified operator evaluation. This is an actionable implementation limitation. The next defensible milestone is end-to-end amplitude-policy control plus independent dense-oracle tests, before promoting rank, nullity, spectral-gap or graph-shell conclusions.


## 8. Raw volume operator inspection

The K5 producer now exposes its raw unsymmetrized triple-product Q through volume123_raw_operator. The existing volume routine uses this shared raw object, records its relative anti-Hermitian defect, and then performs the historical Hermitian completion. A dedicated test file checks four low-spin triples before symmetrization, verifies the derived sqrt(abs(Q)) matrix is Hermitian/positive semidefinite, and checks diagnostic consistency. This tests finite-sector implementation details only; full production HH and continuum claims remain open.


## 9. Local link-hit adjoint convention

A dedicated test now checks the actual production Clebsch-Gordan link-hit matrices for the doubled-spin transition 1 -> 2 across all four endpoint channel pairs. Orientation reversal with endpoint-channel swap must equal the conjugate transpose and preserve normalization. Path and operator-sequence adjunction are tested as involutions. This is a local convention check, not full primitive/T/H adjoint certification.
