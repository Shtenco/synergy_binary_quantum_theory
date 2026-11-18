# Exact pushforward of independent local Haar intertwiners under vertex fusion

## Scope

This note studies one explicitly specified finite ensemble:

- two four-valent vertices;
- spin \(j=1/2\) on all six half-links;
- each local tensor is an independent Haar-random unit vector in the
  two-dimensional four-valent SU(2)-invariant intertwiner space;
- the vertices are fused by contracting one pair of half-links with the
  normalized invariant singlet.

This is a mathematical toy measure chosen to test closure under fusion. It is
not identified with the physical BQG history measure.

## Exact linear map

Let \(\mathcal H_4^{\mathrm{inv}}\) be the four-valent intertwiner space and
\(\mathcal H_6^{\mathrm{inv}}\) the six-valent invariant space. Their
dimensions for spin \(1/2\) are

\[
\dim \mathcal H_4^{\mathrm{inv}}=2,
\qquad
\dim \mathcal H_6^{\mathrm{inv}}=5.
\]

The singlet contraction defines a linear map

\[
C:\mathcal H_4^{\mathrm{inv}}\otimes
\mathcal H_4^{\mathrm{inv}}
\longrightarrow \mathcal H_6^{\mathrm{inv}}.
\]

In orthonormal bases, the computed matrix \(M\) satisfies

\[
M^\dagger M=\frac14 I_4.
\]

Thus \(V=2M\) is an isometry. The image has dimension four, strictly smaller
than the five-dimensional six-valent invariant space.

For independent local Haar vectors \(|\psi\rangle,|\phi\rangle\), their
product second moment is \(I_4/4\). Since the fused output is normalized by
\(V\), its exact pushforward covariance is

\[
\rho_{\mathrm{push}}
=V\frac{I_4}{4}V^\dagger
=\frac14 P_{\mathrm{image}}.
\]

A Haar vector drawn from the *full* six-valent invariant space instead has
covariance

\[
\rho_{\mathrm{Haar},6}=\frac15 P_6.
\]

Because \(P_{\mathrm{image}}\) has rank four while \(P_6\) has rank five,
these measures cannot coincide.

## Numerical verification

| Diagnostic | Result |
|---|---:|
| Local input dimension | 4 |
| Full six-valent invariant dimension | 5 |
| Fusion-map rank | 4 |
| Fusion-map singular values | 0.5, 0.5, 0.5, 0.5 |
| \(\|M^\dagger M-I/4\|_{\max}\) | below 2e-16 |
| Pushforward covariance trace | 1 |
| Pushforward covariance nonzero eigenvalues | four values of 0.25 |
| Pushforward covariance purity | 0.25 |
| Full six-valent Haar covariance purity | 0.20 |
| Frobenius distance of covariances | \(\sqrt{0.05}\approx0.2236068\) |
| Fusion image outside invariant subspace | below 2e-15 |

The result is exact up to floating-point linear algebra. The image rank
defect, not a Monte Carlo fit, establishes the non-Haar pushforward.

## Interpretation and limits

This gives a concrete non-closure result for the *independent local Haar
ansatz*: after fusing two vertices, the induced distribution is supported on
a four-dimensional subspace of the five-dimensional six-valent invariant
space. Replacing the fused block by a fresh Haar-random six-valent intertwiner
would therefore change the measure.

It does **not** show that the physical BQG measure fails to close, because
that measure has not yet been specified. The result instead identifies a
specific correction that any candidate RG construction must account for:
the induced block measure retains information about the microscopic
factorization and cannot automatically be replaced by the full invariant-space
Haar measure.

Next extensions should test higher spins, larger blocks, and correlated
microscopic measures, then determine whether a finite family of induced
couplings can parameterize the pushforward measure.

## Higher-spin extension and analytic structure (verified 2026-11-18)

The exact map was extended to \(j=1,3/2,2\), without constructing a dense
six-leg invariant projector. The six-valent singlet multiplicity was computed
by the exact SU(2) tensor-product recursion, while the image columns were
checked directly against all six total generators.

| Spin | Product input dimension \(d^2\) | Six-valent invariant dimension \(D_6\) | Fusion-map rank | Non-Haar covariance distance |
|---|---:|---:|---:|---:|
| j=1/2 | 4 | 5 | 4 | 0.223607 |
| j=1 | 9 | 15 | 9 | 0.210819 |
| j=3/2 | 16 | 34 | 16 | 0.181902 |
| j=2 | 25 | 65 | 25 | 0.156893 |

For every spin, all singular values of the contraction map equal \(1/d\),
where \(d=2j+1\), and its Gram matrix obeys

\[
C^\dagger C=\frac{1}{d^2}I_{d^2}.
\]

This identity also has a short representation-theoretic explanation. For
two four-valent invariant tensors, the partial transition operator on any
single irreducible spin-\(j\) leg commutes with the SU(2) representation.
Schur's lemma therefore makes it proportional to the identity; its trace
fixes the coefficient to the original intertwiner overlap divided by \(d\).
The two ends of the singlet contraction contribute one such factor each,
giving \(C^\dagger C=I/d^2\).

Consequently the normalized map \(V=dC\) is an isometry, and the pushforward
covariance is maximally mixed only on its \(d^2\)-dimensional image. Since
\(D_6>d^2\) in all four tested sectors, it differs from the full
six-valent Haar covariance. In particular,

\[
\operatorname{Tr}(\rho_{\rm push}^2)=\frac{1}{d^2},
\qquad
\operatorname{Tr}(\rho_{\rm Haar,6}^2)=\frac{1}{D_6},
\qquad
\|\rho_{\rm push}-\rho_{\rm Haar,6}\|_F
=\sqrt{\frac{1}{d^2}-\frac{1}{D_6}}.
\]

Verified run:
[44379209043](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44379209043),
commit `50c5de4affb96a2be80c24dae6e716d96126801a`. All 50/50 workflow
steps passed; the dedicated pushforward test module reports 4/4 passing,
and the exported JSON (artifact ID `14463710973`) reports all four spin
controls as `passed: true`.

The conclusion remains restricted to independent local Haar measures and
two-vertex one-edge fusion. It is a concrete non-closure example for that
ansatz, not a claim that the physical BQG measure is Haar or fails to close.
