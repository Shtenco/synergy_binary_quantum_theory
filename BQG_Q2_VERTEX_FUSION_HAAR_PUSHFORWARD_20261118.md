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
