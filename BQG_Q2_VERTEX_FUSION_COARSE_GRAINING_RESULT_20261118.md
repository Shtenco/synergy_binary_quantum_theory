# Exact vertex-fusion coarse-graining gate (2026-11-18)

## Definition

This gate performs a genuine finite tensor-network vertex blocking, not merely a restriction of edge-angle coordinates.

Given a partition of the microscopic vertices into disjoint nonempty blocks
\(V=\bigsqcup_\alpha C_\alpha\), each coarse tensor is formed by contracting
all singlet edges whose two endpoints lie in the same block:
\[
T_\alpha =
\operatorname{Contr}_{E_{\mathrm{int}}(C_\alpha)}
\left[
\bigotimes_{v\in C_\alpha}T_v
\otimes
\bigotimes_{e\in E_{\mathrm{int}}(C_\alpha)}S_e
\right].
\]
Every cut half-link and every original external half-link remains an open leg
of the fused tensor. Edges crossing between two blocks remain explicit
inter-block singlet contractions. Boundary legs are ordered by their original
(vertex, half-link) labels, independently of block storage order.

## Exact finite statement

The microscopic and fused networks differ only by the parenthesization of
the same finite tensor contraction. Therefore their unnormalized boundary
amplitude tensors agree exactly (up to floating-point contraction error).
After normalization, the boundary pure states agree up to a global phase.

For a retained inter-block edge, the same infinitesimal SU(2) generator is
inserted on the corresponding open leg of the fused tensor. Thus the
microscopic derivative and coarse derivative are the same contraction
re-associated, and the QFI matrices for retained inter-block edge twists
must agree:
\[
F^{\mathrm{coarse}}_{ab}
=
F^{\mathrm{micro}}_{e_a e_b},
\qquad e_a,e_b\in E_{\mathrm{cut}}.
\]
This is a stronger, testable identity than comparing unrelated spectral
summaries. It does not assert that QFI for twists on edges contracted inside
a block survives as a coarse parameter.

## Coverage

The regression suite covers path, star, and cycle fixtures, two spins
\(j=1/2,1\), singleton blocking, and nontrivial vertex partitions. It checks
boundary-state equality up to global phase, retained-edge QFI equality,
partition validation, and duplicate half-link rejection.

## Interpretation boundary

This is an exact finite-dimensional tensor-network identity. It is not yet a
renormalization-group fixed point, a scale-independent spectral law, a
spin-summed BQG amplitude, a boundary partial trace, or a continuum-limit
result. Those require additional physical choices: a coarse observable map,
a measure/weight on microscopic configurations, treatment of internal
degrees of freedom, and a specified refinement trajectory.

## Verified CI result

GitHub Actions run [44359457472](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44359457472)
completed successfully on commit `3b0b20ec9c37773252d474050c5031ac95668343`.
All 41 workflow steps passed, including the vertex-fusion regression suite
and JSON export. The uploaded `q2-pl-graph-export` artifact was downloaded
and the `Q2_VERTEX_FUSION_COARSE_GRAINING.json` payload was independently
parsed as valid JSON.

For the exported (j=1/2) deterministic cases:

| Graph | Boundary amplitude error | Maximum QFI error |
|---|---:|---:|
| path_3 | (1.497\times10^{-16}) | (3.331\times10^{-16}) |
| path_4 | (1.382\times10^{-16}) | (1.110\times10^{-16}) |
| path_5 | (1.486\times10^{-16}) | (1.821\times10^{-17}) |
| star_4 | (1.299\times10^{-16}) | (1.110\times10^{-16}) |
| cycle_4 | (1.646\times10^{-16}) | (2.220\times10^{-16}) |

These are floating-point residuals for the tested finite cases, not a
continuum or renormalization-group statement. The regression suite additionally
checks (j=1), singleton blocking, invalid partitions, and duplicate half-link
rejection.
