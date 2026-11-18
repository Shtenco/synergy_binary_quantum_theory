# BQG QFI × Q2 emergent-graph audit (2026-11-18)

**Branch:** `research/bqg-qfi-q2-emergent-graph-audit-20261118`  
**Base:** `research/bqg-intrinsic-clock-lorentzian-bridge-20261116` at `a0662f45f2fcfec6acae27a4a5091439cc133d5a`  
**Status:** source audit and falsification specification. No new physical result or PASS is claimed.

## Executive finding

The repository already contains an exact QFI-weighted spectral-dimension calibration in
`BQG_QFI_SPECTRAL_DIMENSION_RESULT.md` and
`scripts/bqg_qfi_spectral_dimension_gate.py`. Repeating its 1D/2D/3D hypercubic
calibration would add no evidence for emergent dimension.

The key limitation is mathematical, not computational:

For a homogeneous graph with edge conductance (g>0),
(L_g=gL_1). Therefore
[
P_g(\tau)=P_1(g\tau),\qquad
d_s^{(g)}(\tau)=d_s^{(1)}(g\tau).
]
Changing a common QFI conductance rescales diffusion time; it does not change
the graph's infrared spectral dimension. On a supplied (D)-dimensional
hypercubic topology the exact result tends to (D), but this does not show
that binary dynamics selected that topology.

**Consequently, QFI calibration + the q=2 causal-volume exponent are not yet
a unified derivation of 3D space.** They are two separate results until the
same dynamically generated graph/state family is used in both calculations.

## Source audit

### Existing results that remain valid within their stated scope

1. `BQG_LOCAL_INFORMATION_RESISTANCE_METRIC.md`: local relative-twist QFI
   defines a positive edge conductance in the tested reduced-state family;
   resistance lengths compose additively on a graph.
2. `BQG_QFI_SPECTRAL_DIMENSION_RESULT.md`: exact homogeneous hypercubic
   heat-kernel formula, infrared calibration (d_s\to D), phase collapse
   under (u=g\tau), and a branching-tree discriminator.
3. `Q2_DIMENSION3_FIXED_POINT_CLOSURE.md`: the frozen q=2 route rewrite has
   exact causal-volume growth (N_g=(4\,8^g+10)/7) and
   (d_g\nearrow3); its canonical PL completion has separately tested
   three-manifold topology.
4. `Q2_RELATIONAL_METRIC_SOURCE_GENERATING_FUNCTIONAL.md`: exact finite
   C8 positive-control source/projector construction, explicitly not the
   physical graviton kernel.

### Missing bridge identified by this audit

The inspected q=2 dimension gate establishes a scalar generation-count
recurrence and the PL completion note establishes finite triangulation/link
properties. Neither is, by itself, a serialized family of the same
q=2-generated adjacency matrices, with edge-to-local-state assignments,
that can be passed unchanged to the QFI Laplacian.

The QFI gate instead supplies regular hypercubic graphs by choosing
(D=1,2,3) in advance. This is a calibration control, not a dimension-blind
test of the q=2 dynamics.

The next implementation must therefore export/construct the actual q=2
route/PL adjacency at each refinement level and preserve vertex/edge
identity through state gluing. It must not replace that graph with a cubic
lattice or infer adjacency from the vertex-count recurrence alone.

## Required common-data pipeline

For each refinement/generation (n), freeze one graph and one quantum state
family:

[
\mathcal G_n=(V_n,E_n),\quad
\rho^{(n)}_{ij},\quad
K^{(n)}_{ij}.
]

For each edge compute the operational QFI
[
F_Q(\rho,K)=
2\sum_{a,b:\lambda_a+\lambda_b>0}
\frac{(\lambda_a-\lambda_b)^2}{\lambda_a+\lambda_b}
|\langle a|K|b\rangle|^2,
\qquad g_{ij}=F_Q/4.
]
Require Hermiticity, unit trace, positive semidefiniteness of every reduced
state, Hermiticity of (K), and (g_{ij}\ge0). A zero-conductance edge
must be reported as a disconnected/degenerate edge; it must not be silently
replaced by an epsilon.

Build the weighted combinatorial Laplacian
[
(L_g)_{ii}=\sum_{j:(ij)\in E_n}g_{ij},\qquad
(L_g)_{ij}=-g_{ij}.
]
Check symmetry, row sums, nonnegative spectrum within numerical tolerance,
and the number of zero modes against graph connected components.

Then compute the normalized heat trace
[
P_n(\tau)=|V_n|^{-1}\operatorname{Tr}e^{-\tau L_g},
\qquad
d_{s,n}(\tau)=-2\,d\log P_n/d\log\tau.
]

For finite connected graphs (d_s(\tau)\to0) at sufficiently large
\(\tau), because the zero mode dominates. Therefore a finite-graph
long-time value is not an infrared dimension. A claimed plateau must lie
between the microscopic cutoff and finite-size/spectral-gap crossover, and
must persist over a preregistered logarithmic window as (n) increases.

## Falsification and acceptance gates

### Gate A — input identity

- The adjacency used for q=2 volume growth, PL checks and QFI diffusion must
  be byte-identical or linked by a documented canonical graph map.
- Record graph hashes, node/edge counts, connected components, degree
  distribution, refinement map and state hashes.
- No coordinates or target dimension may be supplied to the graph generator.

### Gate B — local quantum validity

- Compute reduced density matrices from the actual glued BQG state, not a
  separately chosen XX-chain state.
- Derive the twist generator from the physical constraint/relational
  observable algebra. Until then, label the present (K=(Z_i-Z_j)/2) as a
  reduced-model probe, not a fundamental BQG observable.
- Test gauge covariance of (ho_{ij}), (K_{ij}), and QFI under the
  relevant SU(2) transformations.

### Gate C — spectral-dimension estimator

- Retain exact eigenspectra where tractable; otherwise document Lanczos
  residuals and stochastic trace errors.
- Compare analytic finite-size controls against the numerical estimator.
- Vary derivative method, logarithmic window, graph size and QFI precision.
- Report (P(\tau)), (d_s(\tau)), local slopes, uncertainties and
  spectral gap, not only one averaged number.
- A plateau must be stable across at least three successive refinements and
  not be a boundary/finite-size artifact.

### Gate D — topology discrimination

Run the identical pipeline on:
1. q=2-generated graph family (primary);
2. 1D/2D/3D periodic graphs (estimator controls only);
3. non-amenable branching graph (negative topology control);
4. degree-preserving rewires of the q=2 graph (topology sensitivity);
5. randomized local gluing variants satisfying the same microscopic rules.

The primary result is whether the q=2 family, without a target-dimension
input, develops a stable finite (d_s\simeq3) window. If it does not, the
hypothesis that this route yields a 3D diffusion geometry is falsified in
that regime, even if its causal-volume exponent remains 3.

### Gate E — independent consistency

Compare, without identifying them by definition:
- causal-volume/Hausdorff exponent;
- QFI-Laplacian spectral dimension;
- PL local-link dimension and topology;
- independent physical Lorentzian causal response, once available.

A mismatch is a result to explain, not a value to average away.

## Numerical and interpretive cautions

- A scalar homogeneous (g) cannot select (D); it only rescales (\tau).
- On an irregular graph, edge-dependent (g_{ij}) can change diffusion
  behavior, but a finite plateau near 3 is not automatically universal.
- The resistance metric uses (\ell_{ij}\propto1/g_{ij}), while the
  conductance Laplacian uses (g_{ij}). This is a deliberate network
  convention; its continuum relation to the metric tensor must be derived
  and checked, not assumed.
- The spectral dimension is not automatically the Hausdorff dimension.
- The heat-kernel parameter (\tau) is diffusion/probe scale, not physical
  time. It must never be identified with the intrinsic clock or Lorentzian
  time without a separate derivation.
- The existing C8 Page–Wootters construction remains a finite positive
  control; its clock step is not yet derived from the full BQG constraints.

## Adjacent physicalization blocker from the active research branch

The latest inspected W5 direct all-M SU(2) workflow
`43876911486` (branch `research/bqg-intrinsic-clock-lorentzian-bridge-20261116`)
is marked failed. Its unit-test/compile stages passed, but one production
cell (node 2, input K=0, output K=0) failed during the expensive actual
sine-(H_E) all-M commutator; most other matrix jobs were cancelled by the
superseding-run cleanup. The aggregate workflow failure must not be described
as a general SU(2) violation, and it is not a completed covariance verdict.
The actual production/shadow matrix comparison remains incomplete.

The intrinsic-clock bridge document on that branch correctly keeps
constraint-support order, relational-history order and Lorentzian propagation
separate. The QFI heat parameter cannot fill the missing physical-clock step.

## Decision

**Do not promote BQG to “space-time derived from binary information” yet.**
The most valuable next code task is a canonical q=2 graph exporter plus a
QFI-on-that-exact-graph runner, with deterministic artifacts and the gates
above. The result must be allowed to fail.

No existing PASS is changed by this audit. No physicalization gate is closed
by this document.
