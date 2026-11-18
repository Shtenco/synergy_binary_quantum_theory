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

The path_5 partition check is part of the reproducible JSON controls. It
enumerates all 52 set partitions (the fifth Bell number, B_5) of the five
microscopic vertices and repeats each case with both the block order and the
vertex order inside each block reversed. For every partition it checks
normalized boundary-state equality up to global phase, equality of
retained-edge QFI against the microscopic calculation, and invariance of
original retained-edge identities under reordering.

The regression test runs this exhaustive check at both \(j=1/2\) and
\(j=1\). In addition, every fused block tensor is checked against the
residual of the diagonal SU(2) Gauss generators acting on all of its open
legs. The workflow exports a separate control JSON for each spin.
Verified workflow:
[44369135373](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44369135373),
commit `796791744445be2f0155cb92cd4e0c3089e40454`. All 45 job steps
completed successfully. Its `q2-pl-graph-export` artifact (ID
`14457977384`) was downloaded and both JSON files were independently parsed:

| Spin | Partitions | Max amplitude residual | Max QFI residual | Max fused-block Gauss residual | Result |
|---|---:|---:|---:|---:|---|
| \(j=1/2\) | 52 | \(3.472\times10^{-16}\) | \(8.882\times10^{-16}\) | \(8.598\times10^{-16}\) | pass |
| \(j=1\) | 52 | \(3.168\times10^{-16}\) | \(2.665\times10^{-15}\) | \(2.153\times10^{-15}\) | pass |

For both spins, the staged-versus-direct path_5 fusion control also passed:
the retained edge identity is `[2]` in both constructions. For \(j=1/2\),
the amplitude residual is \(1.629\times10^{-16}\) and QFI residual is
\(1.348\times10^{-32}\); for \(j=1\), the corresponding residuals are
\(1.787\times10^{-16}\) and \(4.441\times10^{-16}\).

The separate regression matrix also tests path, star, and cycle fixtures,
singleton blocking, invalid partitions, duplicate half-link rejection, and
ordinary nontrivial partitions at both spins. These results establish
numerical agreement for the tested finite-dimensional cases, not
exhaustiveness over all spin labels, graph topologies, or intertwiner
configurations.

Thus the tested fused tensors remain SU(2)-invariant under the diagonal action
on their retained open legs, within floating-point residuals. This confirms
closure of the local Gauss constraint under this exact finite contraction.

These checks strengthen the finite algebraic identity and software
implementation. They still do not establish a physical RG transformation:
the map currently contracts internal singlets exactly and retains cut legs,
without tracing or summing internal degrees of freedom, specifying a
coarse-graining measure, or deriving scale-dependent effective couplings.

## Two-spin exhaustive partition audit (verified 2026-11-18)

The exhaustive path_5 control was extended to execute at both \(j=1/2\)
and \(j=1\), not merely the smaller spin. The workflow now exports two
separate JSON records and includes both in the artifact:

- \(j=1/2\): all 52 Bell partitions pass; maximum boundary-amplitude
  residual \(3.4724549033\times10^{-16}\), maximum retained-edge QFI
  residual \(8.8817841970\times10^{-16}\).
- \(j=1\): all 52 Bell partitions pass; maximum boundary-amplitude
  residual \(3.1681352919\times10^{-16}\), maximum retained-edge QFI
  residual \(2.6645352591\times10^{-15}\).

The spin-one staged/direct path_5 check also passes: boundary-state residual
\(1.7873771438\times10^{-16}\), QFI residual
\(4.4408920985\times10^{-16}\), with the same retained edge identity [2].

The result was verified from the downloaded and independently parsed
workflow artifact produced by
[GitHub Actions run 44368258342](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44368258342),
commit `b8be1972464190512885baeaf2306ba67a712bc8). All 44 workflow steps
completed successfully; both JSON controls report `passed: true`.

This closes the finite combinatorial partition audit for path_5 at these two
spin labels. It does not close the broader physical coarse-graining problem:
other graph families, larger spins, boundary reductions, internal-state
sums, and a specified scale-dependent effective dynamics remain separate
questions.

## Exhaustive audit across all graph fixtures (verified 2026-11-18)

The exhaustive partition routine was generalized from path_5 to every graph
fixture currently in the gate: path_3, path_4, path_5, star_4, and cycle_4.
For each graph it enumerates the complete Bell set of vertex partitions:
5 + 15 + 52 + 15 + 15 = 102 partitions per spin. Every partition is
re-evaluated after reversing block order and the vertex order inside each
block. The audit also records the diagonal SU(2) Gauss-law residual of each
fused block.

Both spin sectors were run independently:

| Spin | Partitions | Max boundary-state residual | Max retained-edge QFI residual | Max fused-block Gauss residual |
|---|---:|---:|---:|---:|
| j=1/2 | 102 | 2.970e-16 | 6.661e-16 | 8.980e-16 |
| j=1 | 102 | 4.335e-16 | 2.665e-15 | 2.535e-15 |

All 204 canonical partitions passed, as did their reversed-order controls.
Thus the artifact covers 408 partition/order evaluations in total. The
staged-versus-direct path_5 fusion check also passes in both spins and retains
the same original edge identity [2].

Verified run:
[44373496107](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44373496107),
commit `71d568ad2f7d978c2e4a252406e691768cc2266c`. The job completed all 45
steps successfully. The downloaded artifact (ID `14459792831`) contains
both spin-specific JSON exports; both were parsed independently and report
top-level `passed: true`.

This closes exhaustive finite partition coverage for the five current graph
fixtures at these two spin values. It does not prove invariance under
arbitrary graph refinement, arbitrary spin, boundary partial trace, or a
physical renormalization-group flow.
