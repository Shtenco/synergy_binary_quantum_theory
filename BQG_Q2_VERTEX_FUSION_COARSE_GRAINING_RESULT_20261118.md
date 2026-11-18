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

## Higher-spin extension on small-boundary graphs (verified 2026-11-18)

To extend beyond \(j=1\) without constructing prohibitively large dense
boundary vectors, the exhaustive audit now also covers path_3 and cycle_4
at \(j=3/2\) and \(j=2\). These fixtures have at most eight external
half-links in the tested contractions. All set partitions were enumerated,
and each was repeated with reversed block and within-block order:

| Spin | Graph | Partitions | Max boundary-state residual | Max QFI residual | Max Gauss residual |
|---|---|---:|---:|---:|---:|
| j=3/2 | path_3 | 5 | 2.585e-16 | 2.665e-15 | 3.052e-15 |
| j=3/2 | cycle_4 | 15 | 1.794e-16 | 1.776e-15 | 3.645e-15 |
| j=2 | path_3 | 5 | 2.169e-16 | 3.553e-15 | 5.390e-15 |
| j=2 | cycle_4 | 15 | 2.763e-16 | 5.329e-15 | 6.424e-15 |

All 40 canonical partitions and all 40 reversed-order controls passed.
The exported JSON reports 80 partition/order evaluations and top-level
`passed: true`.

Verified run:
[44373952171](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44373952171),
commit `566b0ee860f06b0ef688d4ba76e887af95c06e19`. All 46 workflow steps
passed. The downloaded artifact (ID `14459628478`) was parsed independently;
the higher-spin JSON values above come directly from that artifact.

Combined coverage is now: all 102 partitions of the five graph fixtures at
each of \(j=1/2,1\), plus all partitions of path_3 and cycle_4 at each of
\(j=3/2,2\). This remains a finite set of exact contraction identities,
not an RG fixed-point or continuum result.

## Mixed boundary states and QFI data processing (verified 2026-11-18)

A further finite control now traces out six of the eight external boundary
half-links and retains a two-leg subsystem. It compares the reduced density
matrix, its parameter derivatives, and the symmetric-logarithmic-derivative
(SLD) mixed-state QFI before and after exact vertex fusion. It also checks
the quantum data-processing inequality in matrix form:
\[
F_{\mathrm{pure}}-F_{\mathrm{reduced}}\succeq 0.
\]

The tested controls are path_3 with blocks \(((0,1),(2,))\), retaining
boundary legs \(((0,0),(0,1))\), and cycle_4 with blocks
\(((0,1),(2,3))\), retaining \(((0,1),(0,3))\). Both are evaluated at
\(j=1/2\) and \(j=1\).

| Spin | Graph | Density-matrix error | Derivative error | Mixed-QFI error | Min eigenvalue of pure − mixed QFI |
|---|---|---:|---:|---:|---:|
| j=1/2 | path_3 | 5.551e-17 | 2.194e-17 | 2.261e-33 | 1.000 |
| j=1/2 | cycle_4 | 5.551e-17 | 3.737e-17 | 4.167e-19 | 0.262 |
| j=1 | path_3 | 2.776e-17 | 1.182e-17 | 1.454e-33 | 2.667 |
| j=1 | cycle_4 | 8.557e-17 | 3.273e-17 | 1.110e-16 | 0.827 |

All four controls pass. In path_3, the selected local reduced subsystem is
insensitive to the retained cut-edge twists to numerical precision, so its
mixed QFI is effectively zero; this is consistent with, rather than a
violation of, QFI monotonicity. The cycle_4 subsystem retains nonzero
information, while the pure-minus-mixed QFI matrix remains positive.

The workflow exported `Q2_VERTEX_FUSION_REDUCED_BOUNDARY_QFI.json` in
[run 44374564138](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44374564138),
commit `a12063d35a6cb167ae653f48c1f0c5114891db1d`. The downloaded artifact
(ID `14460102381`) was independently parsed; its reduced-boundary control
reports top-level `passed: true`. The test step and export step both passed.

This adds a verified mixed-state boundary reduction and data-processing
check. It still does not supply a physical measure over microscopic states,
an effective-coupling flow, or an RG fixed point.

## Exhaustive mixed-boundary audit and pure-state limit (verified 2026-11-18)

The reduced-boundary SLD-QFI check was extended to every set partition of
path_3 and cycle_4, for both \(j=1/2\) and \(j=1\), including reversed
block and within-block ordering. This covers 5 + 15 = 20 partitions per
spin, or 40 canonical partitions and 40 reordered controls.

| Spin | Graph | Partitions | Max reduced-density error | Max derivative error | Max mixed-QFI error | Min eigenvalue of pure − mixed QFI |
|---|---|---:|---:|---:|---:|---:|
| j=1/2 | path_3 | 5 | 1.665e-16 | 2.776e-17 | 2.567e-32 | 0.000 |
| j=1/2 | cycle_4 | 15 | 2.220e-16 | 5.554e-17 | 9.080e-19 | 0.000 |
| j=1 | path_3 | 5 | 2.776e-17 | 2.861e-17 | 1.151e-31 | 0.000 |
| j=1 | cycle_4 | 15 | 1.112e-16 | 4.163e-17 | 1.249e-16 | 0.000 |

All controls pass. The minimum eigenvalue is zero to displayed precision,
consistent with the positive-semidefinite data-processing difference; no
strictly positive information loss is claimed for every partition.

A separate pure-limit regression retains all eight path_3 boundary legs at
\(j=1/2\), so no partial trace is performed. The SLD mixed-state QFI agrees
with the direct pure-state QFI within the test tolerance. This independently
checks the mixed-state QFI implementation at its rank-one boundary.

Verified run:
[44375574126](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44375574126),
commit `6b770909a2417bbbb0ba730062aba38ca163730d`. The dedicated vertex-fusion
test step reports **34 passed**; all 47 workflow steps passed. Artifact
`q2-pl-graph-export` (ID `14461372746`) was downloaded and parsed; the
reduced-boundary JSON reports top-level `passed: true`.

This closes the finite reduced-boundary SLD-QFI consistency check for the
tested partitions and spins. It still does not derive a scale-dependent
effective action, a microscopic statistical measure, or a continuum RG flow.

## Full CI verification of the mixed-state extension (verified 2026-11-18)

The latest workflow run includes the corrected pure-limit cross-check
(`SLD QFI = pure-state QFI` when no boundary leg is traced) and exhaustive
mixed-boundary controls over every set partition of path_3 and cycle_4,
including reversed block order and reversed order inside blocks.

The exported mixed-state control contains 80 canonical partition cases:
5 + 15 partitions for each of two graph fixtures and two spin values
\(j=1/2,1\). Reversed-order controls are also evaluated, for 160
partition/order evaluations total. All report `passed: true`.

| Spin | Graph | Partitions | Max reduced-state error | Max derivative error | Max mixed-QFI error |
|---|---|---:|---:|---:|---:|
| j=1/2 | path_3 | 5 | 1.665e-16 | 2.776e-17 | 2.567e-32 |
| j=1/2 | cycle_4 | 15 | 2.220e-16 | 5.554e-17 | 9.080e-19 |
| j=1 | path_3 | 5 | 2.776e-17 | 2.861e-17 | 1.151e-31 |
| j=1 | cycle_4 | 15 | 1.112e-16 | 4.163e-17 | 1.249e-16 |

For each tested reduced subsystem, the minimum eigenvalue of the
pure-minus-mixed QFI matrix was nonnegative at the reported precision
(`0.0` in the JSON summaries), consistent with QFI monotonicity under
partial trace. The independent no-trace test also confirms that the SLD
mixed-state formula reduces to the pure-state QFI formula when the complete
boundary is retained.

Verified CI run:
[44375574126](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44375574126),
commit `6b770909a2417bbbb0ba730062aba38ca163730d`. The workflow completed
47/47 steps successfully, including the dedicated fusion test step and all
four JSON exports. The downloaded artifact (ID `14461372746`) was parsed
independently; all four vertex-fusion JSON files report top-level
`passed: true`.

This validates exact finite tensor-network identities and finite-dimensional
QFI data processing for the enumerated graph fixtures, spins, states, and
subsystems. It does not establish a scale-dependent effective action, a
nontrivial RG transformation on couplings, a fixed point, universality, or
a continuum limit. Those require a specified microscopic ensemble/measure,
coarse variables and induced effective measure, and an explicit scale-flow
observable.

## Exhaustive small-boundary subsystem sweep and batched implementation (verified 2026-11-18)

The reduced-boundary analysis was extended from two selected subsystems to
every one-leg and two-leg subsystem of the eight external half-links in
path_3 and cycle_4. For each subsystem, every vertex partition is tested,
including reversal of block order and the order inside each block. Both
\(j=1/2\) and \(j=1\) are included.

| Quantity | Verified value |
|---|---:|
| Graph fixtures | path_3, cycle_4 |
| Spin sectors | \(j=1/2,1\) |
| External half-links per graph | 8 |
| Retained subsystems per graph and spin | 36 |
| Total retained subsystems | 144 |
| Partition/order evaluations | 2,880 |
| Maximum reduced-density error | 3.886e-16 |
| Maximum reduced-derivative error | 1.403e-16 |
| Maximum mixed-QFI error | 6.106e-16 |
| Minimum eigenvalue of pure-minus-reduced QFI | 0.0 |

All sweep controls report `passed: true`. The implementation batches many
retained subsystems against one fused/microscopic contraction per partition,
rather than rebuilding the same tensor-network state for every subsystem.
A dedicated regression test compares the batched results against the
single-subsystem reference implementation. The complete vertex-fusion test
module reports **36 passed** in 46.02 seconds.

Final verified run:
[44377801973](https://github.com/Shtenco/synergy_binary_quantum_theory/actions/runs/44377801973),
commit `2d2f40b90c58d294dd4ec0f7f7dfa6a0180e3bc9`. All 48/48 workflow
steps passed. The downloaded artifact (ID `14462222520`) was parsed
independently; all five vertex-fusion JSON exports report top-level
`passed: true`.

The nonnegative QFI-gap check is a finite-dimensional data-processing
control, not evidence of a physical RG flow. In particular, exact vertex
fusion is still an exact tensor contraction/reorganization; it does not by
itself integrate out a probabilistically specified set of microscopic
degrees of freedom or define running couplings. The next mathematical
milestone is to specify an ensemble/measure and a coarse-variable map, derive
the induced effective boundary measure, and then test whether any
scale-dependent couplings close under repeated blocking.
