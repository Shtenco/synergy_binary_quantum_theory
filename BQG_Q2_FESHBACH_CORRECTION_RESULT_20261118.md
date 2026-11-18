# Q2 PL-to-causal finite Feshbach correction — numerical audit (2026-11-18)

## Scope and provenance

This note records an independent numerical reproduction of the frozen generation-0/1/2 combinatorial construction and the matrix-free zero-energy Feshbach calculation now implemented in `scripts/q2_causal_pl_correspondence_gate.py --feshbach`.

The independent calculation rebuilt the PL barycentric flags, causal-edge fibres, and 4-regular PL dual adjacency from the repository definitions, then used SciPy conjugate gradients on the projected operator. It is a finite graph calculation, not a physical QFI result. The newly added repository implementation and tests have been committed; the combined-status endpoint returned no status checks for the latest branch commits at the time of this note, so CI is **pending/unverified** here.

## Definitions

Let (L=4I-A_{PL}) be the unweighted PL tetrahedral-dual graph Laplacian. Let (U) be the isometry whose column for causal edge (e) is constant (1/\sqrt{|F_e|}) on its fibre (F_e=\Phi^{-1}(e)), and zero elsewhere. Set (P=UU^T), (Q=I-P).

The exact off-diagonal block and zero-energy Feshbach correction are
[
B=QLU,\qquad
C=B^T(Q L Q|_{Q\mathcal H})^{-1}B.
]
The calculation solves
[
(QLQ+P)X=B,
]
which is nonsingular on the full space and has (X\in Q\mathcal H) because (B\in Q\mathcal H). Then (C=B^TX). Only the sparse PL adjacency is applied; no dense (N_{PL}\times N_{PL}) matrix is formed.

The Galerkin compression is (A_{PP}=U^TLU), and the zero-energy effective operator is (A_{eff}=A_{PP}-C).

## Numerical results

All operator norms below use the same unweighted PL Laplacian normalization (L=4I-A_{PL}).

| Quantity | Generation 1 | Generation 2 |
|---|---:|---:|
| PL dual vertices | 384 | 9,216 |
| Causal coarse modes | 8 | 64 |
| \(\|B\|_F\) | 2.9814239700 | 9.9889530569 |
| \(\|B\|_2\) | 1.8856180832 | 1.8328101444 |
| \(\|A_{PP}\|_2\) | 4.6666666667 | 4.5162805888 |
| \(\|B\|_2/\|A_{PP}\|_2\) | 0.4040610178 | 0.4058229130 |
| \(\|C\|_2\) | 0.7843137255 | 1.7366734310 |
| \(\|C\|_F\) | 1.1748865479 | 5.7083343485 |
| \(\|C\|_2/\|A_{PP}\|_2\) | 0.1680672269 | 0.3845362122 |
| Lowest Q-sector eigenvalue estimate | 0.1522409350 | 0.0117190632 |
| Condition estimate for \(QLQ+P\) | 52.5483 | 682.6484 |
| Maximum linear-solve residual | \(<10^{-14}\) | \(1.4\times10^{-10}\) |

Generation 1 has exact control values:
[
\|B\|_F^2=\frac{80}{9},\quad
\|C\|_2=\frac{40}{51},\quad
\frac{\|C\|_2}{\|A_{PP}\|_2}=\frac{20}{119}.
]

At both nontrivial levels, (C\) is positive semidefinite up to numerical precision, with one constant zero mode. (A_{eff}\) is also positive semidefinite up to numerical precision and has one zero mode. The largest eigenvalues of (A_{eff}\) are 4.5333333333 and 4.0683984560 at generations 1 and 2 respectively.

Generation zero is a useful null control: all PL tetrahedra belong to one causal fibre, (U) is the normalized constant vector, (LU=0), and therefore (B=C=0).

## Interpretation

1. The finite Galerkin compression is not the exact autonomous coarse dynamics at generations 1 or 2. The actual operator-norm leakage is about 40.4% of the Galerkin operator norm at both levels; it does not show a decreasing trend between these two refinements.
2. The zero-energy Feshbach correction is already 16.8% of the Galerkin operator norm at generation 1 and 38.5% at generation 2. It is not a negligible finite-size correction in this normalization.
3. The Q-sector gap decreases by a factor of about 13 between generations 1 and 2, while the correction ratio grows by a factor of about 2.29. This is consistent with the excluded sector becoming more important, but two refinement steps cannot establish an asymptotic law.
4. The ratios are dimensionless and use a common graph-Laplacian normalization. They are useful finite-level diagnostics, but they do not establish a continuum scaling limit, a physical relational QFI metric, or the correct physical renormalization.

## Claim boundary and next gate

**Established numerically for the two finite levels:** the projected off-diagonal block is nonzero; the zero-energy Feshbach correction is nonzero and positive; its relative operator norm is substantial; the Q-sector solve is numerically well-conditioned enough to converge at these two sizes, though the condition estimate worsens sharply.

**Not established:** convergence or divergence as generation tends to infinity; whether the Q-sector gap closes with a power law; whether a rescaled Feshbach kernel converges to a local or nonlocal continuum operator; whether the correction survives the physical QFI and gauge reductions.

The next required calculation is generation 3 using a streaming PL dual-incidence representation, followed by a refinement study of the normalized Feshbach kernel and its low-eigenvalue spectrum. The current exact correspondence exporter remains capped at generation 2.
