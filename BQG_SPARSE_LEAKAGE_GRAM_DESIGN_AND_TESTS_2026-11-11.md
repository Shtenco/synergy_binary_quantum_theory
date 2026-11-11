# Coherent sparse leakage Gram: construction and acceptance tests — 2026-11-11

## Purpose

This note specifies a small reusable numerical kernel for the most important K5 shell correction identified in the consolidated audit. It does not claim that the current K5 Hamiltonian output pipeline is already complete or unpruned.

Implementation:
- [scripts/bqg_sparse_leakage_gram.py](scripts/bqg_sparse_leakage_gram.py)

The module is deliberately decoupled from the current K5 operator code. Its input contract is strict: for every shell basis column (e_j), the caller supplies the **complete, coherently aggregated** sparse expansion of (H e_j), keyed by canonical output basis-state keys. The module itself performs no pruning.

## 1. Algebraic construction

Let (P) be the orthogonal projection onto an (n)-dimensional shell, and (Q=I-P). For the full action (H P), define
[
C=PHP,qquad L=QHP.
]
Then (HP=C+L), with (PC=C), (QL=L), and orthogonal ranges. Hence
[
P H^\dagger H P=(HP)^\dagger(HP)
=C^\dagger C+L^\dagger L.
]
This identity does not require (H=H^\dagger). It does require the correct Hilbert inner product, a complete action, and an orthogonal projection (or a Gram-weighted equivalent).

For shell columns (j,k), with outside basis keys (s\notin S),
[
(G_{\rm leak})_{jk}
=\langle Le_j,Le_k\rangle
=\sum_{s\notin S}\overline{L(s,j)}L(s,k).
]
The implementation groups amplitudes by the same outside key and accumulates the coherent outer product. It therefore preserves off-diagonal interference between distinct input columns that reach a common outside state.

The corrected pullback is
[
M_{\rm pullback}=C^\dagger C+G_{\rm leak}.
]

## 2. Why column-wise leakage norms are insufficient

A diagnostic containing only
[
\ell_j^2=\|Le_j\|^2
]
stores only the diagonal entries of (G_{\rm leak}). It loses the cross terms
[
\langle Le_j,Le_k\rangle,quad j\ne k.
]
Those terms affect eigenvectors, nullity, and the corrected master. They cannot be reconstructed from per-column norms.

The module's embedded two-column test uses
[
Le_0=|x\rangle,\qquad Le_1=i|x\rangle,
]
so
[
G_{\rm leak}=\begin{pmatrix}1&i\\-i&1\end{pmatrix}.
]
A diagonal-only implementation would return the identity and fail to reproduce the exact full-output Gram.

## 3. Scale-invariant spectral classification

The previous K5 projector code used an eigenvalue scale with a unit floor. This is not invariant under (M\mapsto\alpha M) when (\|M\|\ll1).

The new helper classifies the spectrum using
[
s=\max_i|\lambda_i|,\qquad
\tau=\max(\tau_{\rm rel},c\,\epsilon_{\rm mach}n)s.
]
It rejects significant negative eigenvalues instead of treating all values below a positive cutoff as zero. An identically zero matrix is handled explicitly. The zero projector should remain invariant under positive rescaling, up to eigensolver roundoff.

This helper is a diagnostic policy, not a proof that a particular physical master is positive. Positivity must follow from the actual operator construction or be independently verified.

## 4. Embedded tests

The module's `_self_test()` checks:
1. A shared outside key produces the expected complex off-diagonal leakage terms.
2. (C^\dagger C+G_{\rm leak}) agrees with an explicitly assembled full-output matrix on a small test problem.
3. Nullity, positive rank, and the zero projector are invariant under positive rescalings (10^{-24},10^{-12},10^{12},10^{24}).
4. A matrix with a significant negative eigenvalue fails the PSD diagnostic.

The self-test is included in the committed module and can be run with:
```bash
python scripts/bqg_sparse_leakage_gram.py
```

**Execution status:** this consolidation turn did not execute the committed repository copy of the module. A local scratch prototype of the same coherent-outer-product and rescaling test logic passed, but the committed file still needs a direct run in the repository environment before it is marked runtime-verified.

## 5. Integration contract for K5

Before connecting this helper to `bqg_k5_graph_sector_master_shell.py`, the producer must satisfy all of the following:

- Return the full sparse action for each input key, with no hidden thresholding in `apply_primitive_oriented`, `apply_T_cached_state`, `apply_H_cached_state`, or `prune_state`.
- Aggregate all histories reaching the same canonical output key before computing any norm or Gram.
- Establish that the canonical spin-network keys form an orthonormal basis for the inner product used. If not, provide the basis Gram and use its induced metric in the leakage contraction.
- Distinguish exact zero coefficients from small but nonzero coefficients; do not silently discard the latter.
- Track any deliberate truncation with a certified residual vector/norm and keep it separate from physical leakage.
- Verify that all output keys from different input columns use identical canonical conventions, orientation, and cylindrical normalization.

Until these producer obligations are met, the helper returns an exact Gram only **relative to the supplied columns**, not necessarily relative to the mathematical BQG operator.

## 6. Integration sequence

1. Run the standalone helper self-test in the target environment.
2. Add action-level tests for the K5 producer with two input columns sharing at least one outside state.
3. Compare the helper's result against direct dense (A^\dagger A) on tiny spin/cutoff sectors.
4. Add a no-pruning mode to the K5 producer and record output support and discarded-norm budgets.
5. Assemble all five (G_v^{\rm leak}) matrices and the corrected pullback master.
6. Compare compressed-only and leakage-corrected spectra/projectors, then repeat under cutoff and shell-depth growth.

**Do not** change the physical status of the K5 shell result until this integration and the independent adjointness/volume checks have passed.
