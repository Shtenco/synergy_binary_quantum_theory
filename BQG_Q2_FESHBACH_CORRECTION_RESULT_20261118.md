# Q2 PL-to-causal finite Feshbach correction — numerical audit (2026-11-18)

## Scope and provenance

This note records an independent numerical reproduction of the frozen generation-0/1/2 combinatorial construction and the matrix-free zero-energy Feshbach calculation now implemented in scripts/q2_causal_pl_correspondence_gate.py --feshbach.

The independent calculation rebuilt the PL barycentric flags, causal-edge fibres, and 4-regular PL dual adjacency from the repository definitions, then used SciPy conjugate gradients on the projected operator. It is a finite graph calculation, not a physical QFI result. The newly added repository implementation and tests have been committed; the combined-status endpoint returned no status checks for the latest branch commits at the time of this note, so CI is **pending/unverified** here.

## Definitions

Let \(L=4I-A_{PL}\) be the unweighted PL tetrahedral-dual graph Laplacian. Let \(U\) be the isometry whose column for causal edge \(e\) is constant \(1/\sqrt{|F_e|}\) on its fibre \(F_e=\Phi^{-1}(e)\), and zero elsewhere. Set \(P=UU^T\), \(Q=I-P\).

The exact off-diagonal block and zero-energy Feshbach correction are
\[
B=QLU,\qquad
C=B^T(Q L Q|_{Q\mathcal H})^{-1}B.
\]
The calculation solves
\[
(QLQ+P)X=B,
\]
which is nonsingular on the full space and has \(X\in Q\mathcal H\) because \(B\in Q\mathcal H\). Then \(C=B^TX\). Only the sparse PL adjacency is applied; no dense \(N_{PL}\times N_{PL}\) matrix is formed.

The Galerkin compression is \(A_{PP}=U^TLU\), and the zero-energy effective operator is \(A_{eff}=A_{PP}-C\).

## Numerical results

All operator norms below use the same unweighted PL Laplacian normalization \(L=4I-A_{PL}\).

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
\[
\|B\|_F^2=\frac{80}{9},\quad
\|C\|_2=\frac{40}{51},\quad
\frac{\|C\|_2}{\|A_{PP}\|_2}=\frac{20}{119}.
\]

At both nontrivial levels, \(C\) is positive semidefinite up to numerical precision, with one constant zero mode. \(A_{eff}\) is also positive semidefinite up to numerical precision and has one zero mode. The largest eigenvalues of \(A_{eff}\) are 4.5333333333 and 4.0683984560 at generations 1 and 2 respectively.

Generation zero is a useful null control: all PL tetrahedra belong to one causal fibre, \(U\) is the normalized constant vector, \(LU=0\), and therefore \(B=C=0\).

## Interpretation

1. The finite Galerkin compression is not the exact autonomous coarse dynamics at generations 1 or 2. The actual operator-norm leakage is about 40.4% of the Galerkin operator norm at both levels; it does not show a decreasing trend between these two refinements.
2. The zero-energy Feshbach correction is already 16.8% of the Galerkin operator norm at generation 1 and 38.5% at generation 2. It is not a negligible finite-size correction in this normalization.
3. The Q-sector gap decreases by a factor of about 13 between generations 1 and 2, while the correction ratio grows by a factor of about 2.29. This is consistent with the excluded sector becoming more important, but two refinement steps cannot establish an asymptotic law.
4. The ratios are dimensionless and use a common graph-Laplacian normalization. They are useful finite-level diagnostics, but they do not establish a continuum scaling limit, a physical relational QFI metric, or the correct physical renormalization.

## Energy dependence and induced coarse kinetic metric

The calculation now also evaluates the first derivative of the Feshbach self-energy:
\[
\Sigma(\lambda)=B^T(H_Q-\lambda I)^{-1}B,\qquad
H_Q=QLQ|_{Q\mathcal H},
\]
\[
\Sigma'(0)=B^TH_Q^{-2}B\succeq0.
\]
It is obtained by a second matrix-free solve with the same positive operator. This is not an arbitrary extra term: expanding the exact Schur complement gives
\[
A_{PP}-\lambda I-\Sigma(\lambda)
=
(A_{PP}-\Sigma(0))
-\lambda(I+\Sigma'(0))+O(\lambda^2).
\]
Thus the first energy-dependent approximation is a generalized eigenvalue problem
\[
A_{eff}(0)v=\lambda Z_{eff}v,\qquad
Z_{eff}=I+\Sigma'(0)\succcurlyeq I.
\]
The code reports the spectrum of this positive induced metric and the generalized spectrum. This captures the leading low-energy normalization induced by eliminated PL modes. It does not replace the exact energy-dependent kernel; its validity is restricted to energies small compared with the Q-sector gap. Because that gap is already small at generation 2, the radius of a useful low-energy expansion may shrink with refinement.
For \(|\lambda|<\Delta_Q=\lambda_{min}(H_Q)\), the Neumann expansion gives the explicit conservative remainder bound
\[
\|\Sigma(\lambda)-\Sigma(0)-\lambda\Sigma'(0)\|
\leq
\frac{\lambda^2\|B\|^2}
{\Delta_Q^3(1-|\lambda|/\Delta_Q)}.
\]
The implementation reports the coefficient \(\|B\|^2/\Delta_Q^3\), so the bound can be evaluated at any energy strictly inside the certified finite-level radius. This bound is deliberately conservative and does not assert that the actual remainder saturates it.


## Claim boundary and next gate

**Established numerically for the two finite levels:** the projected off-diagonal block is nonzero; the zero-energy Feshbach correction is nonzero and positive; its relative operator norm is substantial; the Q-sector solve is numerically well-conditioned enough to converge at these two sizes, though the condition estimate worsens sharply.

**Not established:** convergence or divergence as generation tends to infinity; whether the Q-sector gap closes with a power law; whether a rescaled Feshbach kernel converges to a local or nonlocal continuum operator; whether the correction survives the physical QFI and gauge reductions.

The next required calculation is generation 3 using a streaming PL dual-incidence representation, followed by a refinement study of the normalized Feshbach kernel and its low-eigenvalue spectrum. The current exact correspondence exporter remains capped at generation 2.


## Frequency-dependent Feshbach response

The finite implementation also computes the first derivative of the self-energy at zero spectral parameter:
\[
\Sigma(\lambda)=B^T(H-\lambda I_Q)^{-1}B,\qquad
H=QLQ|_{Q\mathcal H},
\]
\[
\Sigma(0)=C,\qquad
\Sigma'(0)=B^TH^{-2}B\succeq0.
\]
It obtains this matrix by a second matrix-free solve: first solve \(HX=B\), then \(HY=X\), and form \(\Sigma'(0)=B^TY\). No dense PL-scale matrix is formed.

For \(|\lambda|<\Delta_Q\), where \(\Delta_Q=\lambda_{\min}(H)\), the resolvent has the convergent expansion
\[
\Sigma(\lambda)=B^TH^{-1}B+\lambda B^TH^{-2}B
+\lambda^2B^TH^{-3}B+\cdots.
\]
The projected eigenvalue equation is
\[
A_{PP}-\lambda I_P-\Sigma(\lambda)=0.
\]
To first order in \(\lambda\), this becomes the generalized eigenvalue problem
\[
A_{\mathrm{eff}}(0)v
=\lambda Z_{\mathrm{eff}}v+O(\lambda^2),\qquad
Z_{\mathrm{eff}}=I_P+\Sigma'(0)\succeq I_P.
\]
Thus the eliminated modes renormalize the coarse spectral normalization, not just the static stiffness. The spectral parameter here is the eigenvalue of the unweighted graph Laplacian; it is not yet a physical frequency or a calibrated spacetime energy.

For \(|\lambda|<\Delta_Q\), where \(\Delta_Q\) is the true Q-sector gap, a rigorous remainder estimate is
\[
\|\Sigma(\lambda)-\Sigma(0)-\lambda\Sigma'(0)\|_2
\leq
\frac{\lambda^2\|B\|_2^2}
{\Delta_Q^3(1-|\lambda|/\Delta_Q)}.
\]
The implementation reports the diagnostic coefficient \(\|B\|_2^2/\widehat\Delta_Q^3\), where \(\widehat\Delta_Q\) is the smallest Ritz value returned by Lanczos. From the recorded finite-level estimates it is approximately
\[
g=1:\ 1.008\times10^3,\qquad
g=2:\ 2.087\times10^6.
\]
This roughly 2,070-fold increase warns that the linearized energy-dependent approximation may be valid only in a progressively narrower spectral window. Important qualification: the Ritz value is an eigenvalue estimate, not a certified lower bound on the true Q-sector gap. Therefore these numerical coefficients are conditioning diagnostics, **not certified upper bounds** on the remainder. The analytic remainder inequality is rigorous when supplied with a proven lower bound \(\Delta_Q\); obtaining such a bound remains open. The coefficient increase is not evidence by itself that the exact self-energy diverges.

The code and tests report the norm, trace, extremal eigenvalues, and linear-solve residual for \(\Sigma'(0)\), as well as the positive induced metric \(Z_{\mathrm{eff}}\) and its generalized spectrum. The code also records Ritz residuals for the smallest and largest eigenvalue estimates. A small Ritz residual checks the returned approximate eigenpair but does not, by itself, certify that no smaller eigenvalue was missed. The true gap lower bound is therefore explicitly marked unavailable, and the numerical remainder coefficient is explicitly marked uncertified. The derivative and generalized-spectrum values are also not yet certified by a completed CI run; they must be read from the generated artifact before being treated as verified numerical results.

## Exact finite variational meaning

The zero-energy Schur complement has an exact minimization interpretation. For every coarse vector \(z\),
\[
z^T A_{\mathrm{eff}}(0)z
=
\min_{q\in Q\mathcal H}(Uz+q)^T L(Uz+q).
\]
The unique minimizing excluded component is
\[
q_*(z)=-(QLQ|_Q)^{-1}QLUz=-H^{-1}Bz.
\]
Completing the square gives
\[
A_{\mathrm{eff}}(0)=A_{PP}-B^TH^{-1}B.
\]
Consequently, at every finite level,
\[
0\preceq A_{\mathrm{eff}}(0)\preceq A_{PP},
\qquad C=B^TH^{-1}B\succeq0.
\]
If the PL dual graph is connected and the constant vector lies in the fibre-constant subspace, then the kernel of \(A_{\mathrm{eff}}(0)\) is exactly the one-dimensional constant coarse mode.

This means the Feshbach operator is not merely a perturbative correction: it is the exact minimum PL Dirichlet energy compatible with fixed coarse data. Galerkin energy is the value obtained by restricting the minimization to \(q=0\), so it is an upper bound. These are finite graph theorems, not numerical conjectures, and they provide direct consistency checks for the implementation.
