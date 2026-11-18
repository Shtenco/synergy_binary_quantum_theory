# Two-sided endpoint QFI and singlet quotient (2026-11-18)

## Exact local identity

The graph contraction uses the normalized invariant tensor
\[
\epsilon=\frac{1}{\sqrt2}
\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]
For every \(U\in SU(2)\),
\[
U\epsilon U^T=\epsilon,
\]
equivalently the singlet is invariant under the simultaneous fundamental
action on its two tensor indices. In the tensor-index convention used by
the contraction, the infinitesimal action of the same Hermitian generator
on both endpoint legs therefore annihilates the contracted state.

Let \(D_u(e)\) and \(D_v(e)\) denote the projective boundary-state
derivatives for the same Pauli axis applied at the two endpoints of edge
\(e\). The exact relation is
\[
D_v(e)=-D_u(e).
\]
It holds with all other node states and all other graph contractions fixed.

## Full multiparameter metric

For \(E\) edges, let \(F\in\mathbb R^{3E\times3E}\) be the existing
correlated edge QFI, including all cross-edge blocks. Interleave endpoint
coordinates edge by edge as \((x_{e,u},x_{e,v})\). The two-sided QFI is
\[
F_{\rm 2s}=
\begin{pmatrix}
F_{ee'}&-F_{ee'}\\
-F_{ee'}&F_{ee'}
\end{pmatrix}_{e,e'}.
\]
Equivalently, define \(S:\mathbb R^{6E}\to\mathbb R^{3E}\) by
\[
(Sx)_e=x_{e,u}-x_{e,v}.
\]
Then
\[
\boxed{F_{\rm 2s}=S^TFS.}
\]
The common endpoint directions \(x_{e,u}=x_{e,v}\) form a \(3E\)-dimensional
structural kernel. The relative directions \(x_{e,u}=-x_{e,v}\) reproduce
the edge QFI up to the chosen normalization of the relative coordinate.
Any further nullity is inherited from the projective boundary-state
stabilizer already present in \(F\).

This is a coordinate quotient of the two endpoint probes induced by the
singlet invariant. It is not evidence for an independently postulated
edge-frame gauge field.

## Numerical gate

The new executable gate builds the two endpoint injections from the exact
half-link index map, pulls back the full complex projective tangent, and
checks:
- opposite endpoint tangent columns for every internal edge;
- the complete signed-lift identity including cross-edge correlations;
- all \(3E\) common-rotation null vectors;
- PSD up to roundoff;
- equality of two-sided QFI rank and edge-QFI rank;
- rank agreement between the realified tangent and QFI using a
  square-root-consistent SVD threshold.

It runs randomized normalized four-valent spin-\(1/2\) intertwiners on the
three-node path and four-cycle fixtures. JSON output is written to
\`verification_results/Q2_RELATIONAL_TWO_SIDED_ENDPOINT_QFI.json\`.

## Scope

The result is exact for the normalized spin-\(1/2\) singlet contraction and
the current tensor-index convention. It does not yet derive the simultaneous
vertex-gauge action on boundary frames, extend the invariant tensor to
higher representations, or produce BQG history-weighted coarse conductances.
