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
partition validation, and duplicate half-link rejection. It also verifies
that two successive blocking steps on path_5 agree with one direct blocking
step, including retained-edge QFI and original edge identity.

## Interpretation boundary

This is an exact finite-dimensional tensor-network identity. It is not yet a
renormalization-group fixed point, a scale-independent spectral law, a
spin-summed BQG amplitude, a boundary partial trace, or a continuum-limit
result. Those require additional physical choices: a coarse observable map,
a measure/weight on microscopic configurations, treatment of internal
degrees of freedom, and a specified refinement trajectory.

## Verified CI result

GitHub Actions run [44360015056](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44360015056)
completed successfully on commit 9a4cc5eef2fc7d238c45a79abc4d3c2a5adf674a.
All 41 workflow steps passed, including the direct and staged vertex-fusion
regression tests and JSON export. The uploaded q2-pl-graph-export artifact
was downloaded and Q2_VERTEX_FUSION_COARSE_GRAINING.json was independently
parsed as valid JSON.

For the exported \(j=1/2\) deterministic cases:

| Graph | Boundary amplitude error | Maximum QFI error |
|---|---:|---:|
| path_3 | \(2.439\times10^{-16}\) | \(5.551\times10^{-16}\) |
| path_4 | \(6.023\times10^{-17}\) | \(2.220\times10^{-16}\) |
| path_5 | \(2.641\times10^{-16}\) | \(2.220\times10^{-16}\) |
| star_4 | \(1.591\times10^{-16}\) | \(2.220\times10^{-16}\) |
| cycle_4 | \(2.202\times10^{-16}\) | \(4.441\times10^{-16}\) |

The two-stage path_5 check gives a boundary-state error of
\(1.629\times10^{-16}\), QFI error of \(1.348\times10^{-32}\), and identical
retained edge set [2] for staged and direct fusion.

These are floating-point residuals for the tested finite cases, not a
continuum or renormalization-group statement. The regression suite additionally
checks \(j=1\), singleton blocking, invalid partitions, and duplicate half-link
rejection.

## Exhaustive partition audit (2026-11-18)

The path_5 partition check has now been promoted from a test-only loop into
the reproducible JSON control exported by the workflow. It enumerates all
(B_5=52) set partitions of the five microscopic vertices and repeats each
case with both the block order and the vertex order inside each block reversed.
For every partition it checks normalized boundary-state equality up to global
phase, equality of the retained-edge QFI against the microscopic calculation,
and invariance of original retained-edge identities under reordering.

Verified workflow:
[44363347302](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44363347302),
commit `27eb01e571f592ca0b7752b1e0371c097b8a3960). All 44 job steps completed
successfully, including the exhaustive vertex-fusion regression test. Its
`q2-pl-graph-export` artifact (ID `14455127392`) was downloaded; the JSON
was parsed independently and reports:

- all five deterministic path/star/cycle controls passed;
- all 52 path_5 partitions passed, including reversed block order;
- maximum normalized boundary-amplitude residual across the exhaustive
  partition control: (3.473\times10^{-16});
- maximum retained-edge QFI residual: (8.882\times10^{-16});
- staged versus direct path_5 fusion: amplitude residual
  (1.629\times10^{-16}), QFI residual (1.348\times10^{-32}), and identical
  retained edge identity ([2]).

The separate CI regression matrix continues to cover both (j=1/2) and
(j=1), including staged-versus-direct fusion. The exported exhaustive
partition control uses (j=1/2); it should not be interpreted as exhaustive
over spin labels or graph topologies.

These checks strengthen the finite algebraic identity and software
implementation. They still do not establish a physical RG transformation:
the map currently contracts internal singlets exactly and retains cut legs,
without tracing or summing internal degrees of freedom, specifying a
coarse-graining measure, or deriving scale-dependent effective couplings.

