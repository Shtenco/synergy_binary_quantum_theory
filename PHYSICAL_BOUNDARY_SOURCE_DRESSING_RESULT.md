# Result: enlarged-projector boundary source dressing

Status: **FINITE OPERATOR SOURCE-DRESSING CONTROL — PASS.**

This result records the deterministic output of `scripts/boundary_projector_source_dressing_gate.py`. It is a positive control for the operator order required by the BQG physicalization programme; it is not a BQG physical zero mode or cosmological prediction.

## 1. Enlarged master versus boundary compression

The four-dimensional control master has spectrum

\[
\operatorname{spec}\mathbb M
=\{0,0,2,5\}
\]

up to machine roundoff, so the finite master gap is

\[
\boxed{\Delta_M=2.}
\]

The two-dimensional boundary block has the strictly positive compressed master

\[
\boxed{
B^\dagger\mathbb M B
=\begin{pmatrix}1&0\\0&5/2\end{pmatrix},
}
\]

hence it has no zero vector.

Nevertheless the exact full-space zero projector has boundary overlap

\[
\boxed{
G_0=B^\dagger P_0B=\frac12 I_2.
}
\]

Thus compression before projection would erase the physical overlap.

## 2. Heat-kernel convergence

For

\[
G_\tau=B^\dagger e^{-\tau\mathbb M}B
\]

the exact finite theorem gives

\[
\|G_\tau-G_0\|\le e^{-\tau\Delta_M}.
\]

At the largest frozen projection depth,

\[
\tau=8,
\]

the run obtained

\[
\boxed{
\|G_8-G_0\|
=5.6267587145\times10^{-8},
}
\]

while

\[
\boxed{
e^{-8\Delta_M}
=1.1253517472\times10^{-7}.
}
\]

So the observed boundary convergence lies inside the independently known gap bound.

## 3. Source dressing survives the enlarged projection

The exact physical observables were chosen as Pauli `X/Z` on the enlarged zero sector. After boundary projection and support whitening,

\[
\bar O_X
=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\qquad
\bar O_Z
=\begin{pmatrix}1&0\\0&-1\end{pmatrix}
\]

up to machine precision.

The finite-source check at

\[
(j_X,j_Z)=(0.37,-0.51)
\]

gave exactly

\[
Z_{num}=1.2051545643301038
\]

and

\[
Z_{exact}
=\cosh\sqrt{j_X^2+j_Z^2}
=1.2051545643301038.
\]

Likewise

\[
W_{num}=W_{exact}=0.18660782787044697.
\]

At `tau=8` the whitened zero-source source Hessian is

\[
\begin{pmatrix}
0.9999998874648374&0\\
0&0.9999998874648404
\end{pmatrix},
\]

so

\[
\boxed{
\|W^{(2)}_8-I\|
=1.1253516263\times10^{-7}.
}
\]

## 4. The negative control

If one first compresses the master and only then exponentiates,

\[
e^{-\tau B^\dagger\mathbb M B},
\]

the result tends to zero because the compressed master is positive definite.

At `tau=8` its operator norm is already

\[
3.3546262790\times10^{-4},
\]

whereas the correct full-space projection retains

\[
B^\dagger P_0B=\frac12I_2.
\]

This directly falsifies the shortcut

```text
compress master to q=2 boundary
-> find no boundary zero mode
-> conclude zero physical boundary amplitude.
```

The correct ordering is

\[
\boxed{
\text{enlarged }\mathbb M
\to P_0\text{ / heat filter}
\to B^\dagger P_0B
\to \text{source projection/whitening}
\to Z[J]
\to W[J].
}
\]

## 5. Scientific consequence for BQG

The previously found full-rank logical Euclidean master does not terminate the physicalization programme. It instead proves that the physical state, if it exists with q=2 boundary support, must be a dressed state in the enlarged Peter-Weyl/graph-changing habitat.

The production calculation must therefore ask whether the genuine BQG master has

\[
\boxed{B^\dagger P_{phys}B\ne0}
\]

under a regulator/refinement-controlled sequence.

Only after that can volume/shape/lapse source insertions be called physical and connected interblock Hessians be promoted toward `Gamma_scalar^(2)`.

## 6. Claim boundary

This PASS certifies the operator/source architecture only. It does not establish that the actual BQG enlarged master possesses a zero sector, that the q=2 boundary overlaps it, or that dark matter, dark energy, lensing, a physical graviton frequency or a cosmological vacuum has been derived.
