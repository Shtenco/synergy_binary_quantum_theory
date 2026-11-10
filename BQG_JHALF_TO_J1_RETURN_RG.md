# First BQG symmetry-matched dynamical RG datum: j=1/2 -> j=1

Status: **FINITE RETURN-KERNEL TEST PASSED; DYNAMICAL RG MATCHING FAILED AT THE TESTED FINITE REGULATOR.** This is not a full-master or continuum RG verdict.

## Exact logical embedding

The fine K5 carrier has five local logical bits with face spin j=1/2 and local labels K in {0,2}. The coarse carrier has face spin j=1 and the exact symmetry-selected S4-[2,2] embedding

\[
|0_L\rangle=|K=2\rangle,\qquad
|1_L\rangle=\frac23|K=0\rangle-\frac{\sqrt5}{3}|K=4\rangle.
\]

Both logical spaces have dimension 2^5=32. No fitted basis rotation is used.

## Operator compared

The committed gate applies the same regulator-safe Euclidean sine-Hermitian operator

\[
H_{01}=H_{E,0}^{\rm sine}+H_{E,1}^{\rm sine}
\]

to each logical basis vector and constructs the positive return kernels

\[
K_f=A_f^\dagger A_f=P_fH_{01,f}^2P_f,\qquad
K_c=A_c^\dagger A_c=P_cH_{01,c}^2P_c.
\]

The comparison is made in the same symmetry-matched 32 logical coordinates:

\[
D=K_c-K_f.
\]

## Verified numerical result

Workflow run [42704882313](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/42704882313), completed successfully on 2026-11-09.

| Quantity | Fine j=1/2 | Coarse j=1 |
|---|---:|---:|
| Logical dimension | 32 | 32 |
| Rank | 32 | 32 |
| Nullity | 0 | 0 |
| Smallest eigenvalue | 4.306075987001578 | 3.876120509435486 |
| Largest eigenvalue | 13.352781352746614 | 6.378228677546204 |
| Condition number | 3.1009163314938317 | 1.6455186731217297 |
| Operator norm | 13.352781352746614 | 6.378228677546203 |

The spectra are strictly positive in this finite calculation. This certifies positivity and full rank of these two finite return kernels, not positivity/gap of the full physical master constraint.

The dynamical mismatch is substantial:

\[
\|K_c-K_f\|_2=6.986118602340976,
\]

\[
\frac{\|K_c-K_f\|_2}{\|K_f\|_2}
=0.5231957610766957,
\]

\[
\frac{\|K_c-K_f\|_F}{\|K_f\|_F}
=0.45292285811472693,
\qquad
\operatorname{Tr}(K_c-K_f)=-124.22786611330929.
\]

Thus the finite return dynamics is not approximately invariant under this first j=1/2 -> j=1 embedding at the tested regulator. The gate's internal passed flag means its finite-dimensional Hermiticity/PSD/spin-support checks passed; it does not mean the RG mismatch is small or that an RG fixed point was found.

## Regulator and scope caveats

The runner uses doubled-spin cutoffs JMAX2_FINE=3 and JMAX2_COARSE=5 (i.e. walls 3/2 and 5/2). The observed output support obeyed the runner's one-hit bounds: fine maximum spin 1 and coarse maximum spin 3/2. The coarse cutoff is conservative; its wall is 5/2, not the observed maximum 3/2.

This is a one-hit Euclidean return-kernel comparison for the selected two nodes. It is not the full graph-changing master residual

\[
R=M_cW-WM_f,
\]

does not measure off-carrier leakage of the full master, and does not establish a continuum/refinement limit. No beta/kappa/hbar physical coefficient fitting, physical joint Gram, or full physical-kernel computation is certified by this result.

## Next required work

1. Repeat at independently justified regulator walls and demonstrate cutoff convergence of both kernels and of the mismatch.
2. Decompose the mismatch into carrier-projected and off-carrier contributions; verify the exact embedding is isometric in the actual inner product used by the operator.
3. Compute the full graph-changing master on both carriers, then evaluate R=M_c W-W M_f and its leakage decomposition.
4. Only after the full master, physical ordering, and normalization are fixed, compare physical gaps and common kernels.

Until those tests pass, the scientifically correct conclusion is: the first finite return-kernel RG comparison is positive and full-rank at both scales, but exhibits a large (~52.3% operator-relative) dynamical mismatch; no RG fixed point has been demonstrated.
