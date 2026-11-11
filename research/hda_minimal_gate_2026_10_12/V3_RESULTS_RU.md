# BQG HDA V3 — further exact calculations (2026-10-12)

**Status:** HDA on full Lorentzian microscopic K5 **OPEN/BLOCKED**. Positive results are restricted to theorems explicitly stated below and synthetic controls. Do not mark physicalization gates closed.

## 1. Ten pair norms are already sufficient

For any common-source residual vectors
`r_ij=([H_i,H_j]/(i hbar)-D_ij)psi`, Gram `G_ab=<r_a,r_b>` satisfies

```text
max_{i<j} ||r_ij||^2 <= lambda_max(G) <= tr(G) <= 10 max_{i<j} ||r_ij||^2.
```

Therefore `lambda_max(G_n)->0` **iff** the maximum of ten squared residual norms vanishes (same positive normalization at every n). Full Gram retains value for multi-lapse comparisons, but diagonalization is not necessary just for convergence. Per-source-vector result; physical operator closure requires uniformity over a common habitat and refinement.

## 2. Strict nodal aliasing counterexample

On five fixed distinct points in R^3, take `x0=0`, `x1=(1,0,0)`, `x2=(0,1,0)`, `x3=(0,0,1)`, `x4=(1/3,1/4,1/5)`, and

```text
N(x)=x^1 product_(j=1..4) (|x-x_j|^2 / |x0-x_j|^2), M(x)=1.
```

Then ALL `N(x_i)=0` exactly; `grad N(x0)=(1,0,0)`, `grad N(x_j)=0` for j>0. If `H[N]=sum_i N(x_i)H_i` uses only these fixed samples globally, then `[H[N],H[1]]=0`. But at Euclidean q=I, independent `beta=N grad M-M grad N=-grad N` acts on smooth probe `f=x0^1` as `D[beta]f=-1`. Thus this **fixed-nodal-only scheme cannot implement continuum HDA on all smooth lapse functions**. This is NOT a no-go for finite-element restricted lapses, graph-changing evaluation at newly generated vertices, or a continuum limit for fixed smooth lapses. Symbolic computations in `v3_core_gate.py`.

## 3. Numeric spectral obstruction from real stored SU(2) data

From `research/audit_2026_10_10/results.json` Euclidean one source column and pair 01, `||Cpsi||=2.879453814704954`, `||P32 Cpsi||=0.5503455636886702`, `||Q32 Cpsi||=2.826371177242493` (hbar=1 convention).

If D01 psi stays within original 32D code, then the associated Gram diagonal is at least `||Q32 Cpsi||^2=7.988374031547116`. Hence **lambda_max(G)>=7.988374031547116** on this single Euclidean source and finite cutoff, conditional on D's scope. Relative squared residual to that commutator >= `0.9634699130467539`. No claim for full H_E+lambda_L H_L, different source or n->infinity.

## 4. Lorentzian grading block, CONDITIONAL on D parity

According to `PETER_WEYL_LORENTZIAN_PARITY_FRONTIER.md`, `Pi E Pi=-E` and `Pi L Pi=L`. If independently derived D is EVEN and psi a parity eigenvector, define `A=[E,E]/i-D`, `B=([E,L]+[L,E])/i`, `C=[L,L]/i`. The vector residual has orthogonal even `(A+lambda^2 C)psi` and odd `lambda B psi` parts.

Hence `||r||^2=||(A+lambda^2 C)psi||^2+lambda^2||Bpsi||^2`, and **Gram(lambda)=G0+lambda^2 G2+lambda^4 G4**. Synthetic exact-parity Hermitian model: max Gram discrepancy 1.7758741495849052e-16. Not an actual computed BQG Lorentzian mixed amplitude; if D flips parity, this reduction cannot be used.

## 5. Nontrivial scale trap

Rescale `H_i^(n)=2^{-n}H_i`, `D_ij^(n)=4^{-n}D_ij`. Then Gram(lambda_max) shrinks as `16^{-n}`, reaching `5.960464477539063e-8` of its value by n=6, with NO improved closure. A dimensionless diagnostic is `Q_n=lambda_max(G_n)/sum_ij (||C_ij psi||^2+||D_ij psi||^2)`, with required positive denominator and a separately nonzero off-shell D probe. Synthetic n=0..6 gives **Q_n=0.3882771742294772** unchanged while the raw Gram shrinks seven orders. This is why bare Gram->0 is not a physical criterion by itself.

## 6. Sparse memory-efficient Gram

Locally reproducible script `streaming_gram.py` streams sorted final basis keys and ten complex C,D amplitudes. Contributions to each (key,pair) are SUMMED before Gram outer product. It needs O(10^2) accumulators rather than a dense 11-million-dimensional matrix. On synthetic 17-key/10-pair amplitudes: relative Gram error 2.2159223497801594e-16. **No physical C,D amplitude streams presently supplied.** Independent D provenance and a common graph-changing state basis are mandatory.

## Sources and status

- [DUAL_K5_HDA_COCHAIN_TARGET.md](../../DUAL_K5_HDA_COCHAIN_TARGET.md): exact graph cochain, not metric sharp.
- [OFF_SHELL_HDA_HABITAT_TARGET.md](../../OFF_SHELL_HDA_HABITAT_TARGET.md): nonconstant vertex probes and independent D required.
- [LORENTZIAN_HH_REACHABLE_SPACE.md](../../LORENTZIAN_HH_REACHABLE_SPACE.md): conservative 11,314,085 support, not actual amplitudes.
- [Euclidean two-hit result](../../research/audit_2026_10_10/results.json): one pair, one column.
- [Parity frontier](../../PETER_WEYL_LORENTZIAN_PARITY_FRONTIER.md): conditional operator grading.

**Open:** full 24-epsilon Lorentzian H on common domain, independent geometric D, all 10 actual full residuals, meaningful renormalized refinement, physical projector/history and interacting positive TT pole. No main branch physicalization status changed.
