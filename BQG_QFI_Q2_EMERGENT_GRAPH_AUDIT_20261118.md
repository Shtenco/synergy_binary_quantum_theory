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
The graph-export stage has since been implemented as separate PL primal/dual,
causal-route and local-Walsh-QFI controls (see the follow-up sections below).
The next highest-value task is to derive a gauge-covariant two-block reduced
state and relational generator on a shared PL face, and then derive the map
from the scalar causal DAG to the PL dual graph. Only after those inputs exist
should the QFI Laplacian be promoted to a physical spectral-dimension test.

No existing PASS is changed by this audit. No physicalization gate is closed
by this document.


## Follow-up implementation: separate PL and causal graph exporters (2026-11-18)

Two distinct deterministic graph families are now exported rather than
conflated:

1. `scripts/q2_pl_graph_export.py`: the 16-cell boundary and two global
   barycentric refinements. CI run `44269971966` passed all 3 tests and
   uploaded the 70,352-byte adjacency artifact.
2. `scripts/q2_causal_graph_export.py`: the directed scalar q=2 route rewrite
   through generation 5. CI run `44270674765` passed all 3 tests and uploaded
   the 635,543-byte causal adjacency/ancestry artifact.

The causal exporter realizes each directed edge (u\to v) as four
labelled alternatives (u\to m_a\to v), with (a=00,01,10,11). It
reproduces (N_5=18,726), (E_5^{active}=32,768), causal path length 32 and
(4^{31}) source-to-sink paths for this explicit independent-branch
realization.

**New structural warning:** in this literal graph realization the degree of
each retained old vertex is multiplied by four at every refinement; the
boundary degrees are (4^g). The direct limit is therefore not uniformly
locally finite. The exporter records this rather than hiding it. A physical
Laplacian/continuum limit now requires a derived vertex measure, a normalized
refinement rule, or a revised microscopic gluing. Any such change must be
retested against the exact causal-volume scaling.

The two exported graph families have different seed sizes and incidence:
the causal rewrite starts from one directed edge, while the PL completion
starts from the 8-vertex/16-tetrahedron 16-cell boundary. Sharing the four
route labels does not define a global incidence-preserving map. Thus Gate A
remains open; QFI on the physical causal/PL graph remains blocked by the
missing canonical gluing map and physical \\(\\rho_e,K_e\\) assignments.

Detailed results:
- `BQG_Q2_PL_GRAPH_EXPORT_RESULT_20261118.md`
- `BQG_Q2_CAUSAL_GRAPH_EXPORT_AUDIT_20261118.md`

The exact causal-volume exponent remains a theorem about the frozen scalar
recurrence. It has not been upgraded to a QFI spectral-dimension result.


## PL tetrahedral dual graph and local QFI control (2026-11-18)

The PL exporter now also exports the tetrahedral dual adjacency. At levels
0, 1 and 2 it has respectively ((V,E)=(16,32),(384,768),(9216,18432)).
All three dual graphs are connected and exactly 4-regular. Each dual edge is
tagged by the omitted face-rank/route label (00,01,10,11); the two
tetrahedra sharing a face agree on that colour. The updated PL CI run
`44271050772` passed all 5 tests and uploaded the 432,503-byte artifact.

A separate local QFI control now uses the actual four Walsh-derived
tetrahedral face-qubit states from `micro_walsh_qgeom_gate.py`. For a pure
face qubit the SU(2)-frame QFI is (F^{(a)}=I-n_an_a^T), and the four-face
sum is exactly ((8/3)I_3).

The same calculation exposes an important gauge null: after the existing
Gauss-singlet projection (weight (2/9)), the global SU(2) generators
annihilate the state and its global gauge-orbit QFI vanishes to numerical
precision. The reduced product-state relative probe
(K_{ij}=(Z_i-Z_j)/2) gives (F_Q=4/3), hence (F_Q/4=1/3) for all six
face pairs. This (1/3) is **not** a physical edge weight: it is a uniform
reduced-model product-state control, and a uniform scalar only rescales
diffusion time.

The local-QFI CI run `44271370461` passed all 3 tests and uploaded its JSON
artifact. Detailed result: `BQG_Q2_WALSH_LOCAL_QFI_CONTROL_RESULT_20261118.md`.

The next physical QFI step is therefore sharply specified: derive a
gauge-covariant two-block reduced state and a relational generator on each
shared tetrahedral face. A global gauge generator is ruled out as the source
of nonzero conductance; the product-state (1/3) control must not be promoted
to a physical weight. The 4-regular PL dual graph is now available as the
candidate adjacency, but the scalar causal rewrite DAG remains a separate
object until a dynamical gluing map is derived.

## Follow-up: finite causal–PL quotient, relational link QFI, and generation-3 preflight (2026-11-18)

The earlier decision above is historical and has been refined by subsequent
work on the active research branch.

### Finite causal–PL correspondence and Galerkin operator

A recursive many-to-one map from PL tetrahedra to active q=2 causal edges is
now implemented and checked through generations 0–2. It preserves PL-dual
adjacency as endpoint incidence of the image causal edges and is natural under
the tested barycentric/causal refinement. It is not an isomorphism, covering,
or exact Markov lumping. The pushforward counting measure
[
\mu_g(e)=|\Phi_g^{-1}(e)|/24^g
]
has total mass 16 and is conserved parent-to-children at both tested
refinement steps. The induced Galerkin quotient has 8 nodes/16 weighted pairs
at g=1 and 64 nodes/180 weighted pairs at g=2. Its finite generalized spectrum
and heat trace are diagnostic only; no spectral plateau or continuum limit
has been established.

### Two-block relational QFI

A concrete two-node spin-network gluing channel is now specified. Two arbitrary
four-valent spin-1/2 SU(2) intertwiners are glued across one shared face with
the normalized epsilon singlet, retaining the six uncontracted boundary
half-links. The resulting boundary density operator is
rho_e=|Psi_AB><Psi_AB|. The inserted relative twist on B's shared half-link
pulls through the B-node Gauss law to the explicit boundary generator
K_e^a=-sum_{r=1}^3 J^a_{B,r}, acting on B's three retained half-links.
The gluing probability is 1/4 and the pure-state relational QFI is
F_e=I_3, so the normalized metric F_e/4=I_3/4. The script checks the
half-link-to-boundary-generator identity directly.
The result is independent of the two intertwiner states in this fixed
spin-1/2 channel. This closes a *specified local relational-QFI control*,
not the full physical graph metric: boundary gauge reduction, tracing,
network correlations, spin sums, and the BQG history-amplitude measure remain
open. In particular, the microscopic (1/4) must not yet be multiplied into
the Galerkin quotient and called a derived coarse BQG conductance.

### Generation-3 capacity

An allocation-free exact f-vector preflight predicts
[
(f'_0,f'_1,f'_2,f'_3)=(40256,261440,442368,221184),
]
with 442,368 PL-dual edges and 512 causal active edges at g=3. The
correspondence gate now switches to NumPy/LAPACK `eigvalsh` at 128+ quotient
nodes, retaining the stdlib Jacobi method as a fallback. The spectral
diagonalization bottleneck is therefore reduced, but a compact/streaming PL
and dual-incidence representation is still required before attempting a
full g=3 correspondence export; its memory and serialization cost has not
yet been benchmarked.

### Execution status

The new QFI and capacity gates, tests, and workflow steps are committed.
The current GitHub combined-status lookup returns no status checks for the
latest branch commit, so the new repository tests and Actions run remain
**unverified** here. Independent local arithmetic and a numerical
cross-check of the two-block contraction reproduce the stated f-vector,
gluing probability (1/4), and (F^{(e)}=I_3); this is not a substitute for
the repository CI result.

Detailed follow-up files:
- `BQG_Q2_CAUSAL_PL_CORRESPONDENCE_GATE_RESULT_20261118.md`
- `BQG_Q2_RELATIONAL_EDGE_QFI_RESULT_20261118.md`
- `scripts/q2_pl_generation3_preflight.py`
- `scripts/q2_relational_edge_qfi_gate.py`


## Follow-up: two-link network QFI and non-additivity (2026-11-18)

A three-node path A-B-C was added, with two normalized spin-1/2 singlet
contractions A3-B0 and B3-C0. The eight retained half-links form a pure
boundary state. Two independent relative twists on B0 and B3 give a full
6-by-6 multiparameter QFI matrix.

For the tested SU(2)-invariant path family, the two diagonal 3-by-3 link
blocks are each I_3, but the cross block is c I_3. It is not zero in general
and depends on the middle intertwiner. The exact finite-dimensional identity
is

```text
c = (4/3) <iota_B | J_B0 . J_B3 | iota_B>
```

The projected pair-Casimir spectrum in the four-valent invariant space is
{-3/4,+1/4}, so c lies in [-1,1/3]. The full QFI eigenvalues are 1+c and
1-c, each triply degenerate. This proves that summing independent two-node
edge QFIs does not in general reproduce the network multiparameter QFI.
A coarse metric must retain the cross terms or justify a controlled
coarse-graining that removes them.

The new gate checks both signed Gauss-law pullthrough identities, the exact
Casimir formula, isotropy, the QFI spectrum's positivity, and the fixed
two-link gluing probability 1/16. An independent local 32-sample numerical
reproduction found Casimir-formula error below 8e-16 and pullthrough residual
below 3e-16. This is a local numerical cross-check, not a repository test or
CI result.

Files:
- `scripts/q2_relational_network_qfi_gate.py`
- `tests/test_q2_relational_network_qfi_gate.py`
- `BQG_Q2_RELATIONAL_NETWORK_QFI_RESULT_20261118.md`

The workflow now includes the network-QFI tests and JSON artifact. The
repository tests and GitHub Actions remain unverified until an actual run
and status are observed. The result is still a finite spin-1/2 path control,
not a full BQG network measure or continuum conductance.


## Follow-up: four-valent star QFI and exact gauge quotient (2026-11-18)

The two-link chain result extends to all four incident half-links of one
four-valent invariant node. For link labels i,j and SU(2) axes a,b, the
local pure-state QFI factorizes exactly as

```text
F_(i,a),(j,b) = C_ij delta_ab
C_ij = (4/3) <iota | J_i . J_j | iota>
```

The diagonal is one. The Gauss constraint implies C times the all-ones
vector is zero, so the 12-parameter star QFI has three exact null directions:
the common SU(2) rotation. Thus four incident links are not four independent
positive scalar conductances. Their relative twists have a coupled 4-by-4
Casimir matrix.

Fixing the fourth link parameter to zero gives the 3-by-3 principal relative
metric C[0:3,0:3], positive definite for generic invariant intertwiners.
An independent 32-state local numerical reproduction found the largest
Gauss-null residual about 1.1e-15, smallest gauge-fixed eigenvalue about
5.73e-3, and pair-metric eigenvalues nonnegative to roundoff.

This is still a local invariant-node metric. A full BQG graph calculation
must combine these node blocks consistently across shared links, perform
boundary gauge reduction and tracing, and define the history-amplitude
measure before any graph conductance or spectral dimension is claimed.

Files:
- `scripts/q2_relational_star_qfi_gate.py`
- `tests/test_q2_relational_star_qfi_gate.py`
- `BQG_Q2_RELATIONAL_STAR_QFI_RESULT_20261118.md`


## Follow-up: explicit finite graph QFI contractions (2026-11-18)

A graph-level tensor-network gate now contracts normalized four-valent
spin-1/2 intertwiners across explicitly identified half-links using normalized
SU(2) epsilon singlets. Each internal edge has an independent three-axis
twist inserted at its first endpoint. The gate computes the complete pure
boundary-state multiparameter QFI, retaining all cross-edge blocks.

Two fixtures are included: a three-node path and a four-cycle. This is the
first direct construction beyond a single star or a two-edge path. The
cyclic fixture is intentionally not forced to inherit the isolated-edge
identity block: additional contractions can change the reduced boundary
state and its diagonal QFI. The output reports each diagonal block and every
cross-edge contribution separately. It accepts only symmetry and
positive-semidefiniteness of the complete QFI, plus nonzero cross terms for
the tested fixtures.

The graph gate is added to the branch workflow and exports
`Q2_RELATIONAL_GRAPH_QFI.json`. The workflow status is not yet verified.
This is still a finite spin-1/2 graph control, not the general BQG network
metric: the next required layer is a global boundary-gauge quotient and a
canonical pushforward to the causal-PL quotient graph.


### Finite graph contraction check: path versus cycle

The explicit graph tensor-network gate has two fixtures: a three-node path
and a four-cycle. Independent local reproduction (eight random invariant
node-state samples per fixture, seed 20261118) confirms:

- Path: each diagonal edge QFI block remains I_3, while cross-edge terms are
  nonzero (Frobenius norm 0.218–2.336).
- Four-cycle: cross-edge terms are larger (Frobenius norm 1.025–3.236),
  and diagonal edge blocks themselves deviate from I_3 by as much as 0.233.
- Both full QFI matrices remain positive in these samples; the smallest
  eigenvalue observed was 0.0464 for the path and 0.1231 for the cycle.
- The four-cycle contraction probability is state-dependent (about
  0.00365–0.00434), unlike the path's fixed 1/16.

This is direct evidence that network topology and the full boundary
contraction affect both off-diagonal and diagonal QFI blocks. Isolated
two-node QFI is not a universal edge weight. These are independent local
numerical checks, not repository CI results.


### Full graph half-link parameter space and Gauss null modes

The graph gate was extended from one twist triplet per internal edge to all
four half-link twist triplets at every node, including boundary legs. This
produces a 12V-by-12V QFI and makes each node's three SU(2) Gauss directions
explicit. The code checks that all 3V nodewise common-rotation vectors are
null. It then fixes the fourth half-link at each node as a reference and
reports the resulting 9V-dimensional principal metric, including any
additional null modes and its numerical rank.

For the tested path and cycle samples the full QFI is positive semidefinite
up to roundoff and the nodewise gauge-null residuals are at machine
precision. The gauge-fixed metric is also PSD but not full rank: the
independent sample showed ranks 21/27 for the path and 24/36 for the cycle.
This means nodewise Gauss reduction is necessary but not sufficient to
guarantee a nondegenerate physical metric. No epsilon regularization is
introduced.

The graph-level gate now tests both edge-parameter QFI and full half-link
QFI. Its artifact remains a finite spin-1/2 control; the causal-PL
pushforward and BQG history measure are still outstanding.


## Follow-up: global vertex-gauge quotient and two-sided endpoint parameters (2026-11-18)

The finite graph QFI work now has a separate global quotient gate,
`scripts/q2_relational_gauge_quotient_gate.py`. It uses the complete
12V-dimensional half-link parameter space and constructs the 3V-dimensional
vertex Gauss subspace explicitly. The orthogonal quotient has dimension 9V
and metric

```text
F_quot = Q^T F_half Q,
```

where Q spans the orthogonal complement of the vertex gauge directions.
The gate checks the matrix identity F_half G = 0, quotient PSD up to
roundoff, and the quotient nullity that remains after removing the known
Gauss modes.

Two reference sections are compared (fixing half-link 0 versus half-link 3
at every vertex). Their raw coordinate matrices are not compared by
ordinary eigenvalues: the test uses generalized spectra with the induced
quotient Gram matrix. This is the coordinate-invariant comparison. The
result therefore tests whether the quotient spectrum is independent of
which half-link is chosen as a reference.

Internal edges now also have an explicit two-sided endpoint parameter map:
three independent axes at each endpoint, or 6E parameters. Its QFI is
defined as the pullback of the full half-link QFI. It is intentionally not
called the global quotient: boundary half-links must remain in the base
space because the vertex Gauss action rotates all four legs.

The new tests and JSON export are wired into
`.github/workflows/q2-pl-graph-export.yml`. The report
`BQG_Q2_RELATIONAL_GAUGE_QUOTIENT_RESULT_20261118.md` records the
definitions and claim limits. Additional quotient null modes remain an
explicit result, not something removed by regularization. This does not yet
derive a physical edge-frame redundancy, a holonomy dynamics, a BQG history
measure, or the causal-PL QFI pushforward.


## Follow-up: projective stabilizer diagnosis of residual quotient nulls (2026-11-18)

The global quotient gate now builds the tangent map of the normalized
contracted boundary-state ray for every half-link/axis parameter. Each
derivative is projected orthogonally to the full complex state vector, so
both radial normalization and global phase are removed. For the realified
tangent matrix R(D), the implementation checks

```text
F_half = 4 R(D)^T R(D).
```

The quotient tangent SVD is compared with the quotient-QFI spectrum. Thus
the nullity left after the known 3V Gauss directions is identified, within
the tested finite model, as the dimension of the infinitesimal stabilizer
of the projective boundary state. The code records its support over nodes
and half-links. This is a state- and contraction-dependent stabilizer
diagnostic, not a declaration that all residual nulls are fundamental
gauge symmetries.

A regression test checks the full matrix identity above. The focused report
`BQG_Q2_RELATIONAL_GAUGE_QUOTIENT_RESULT_20261118.md` describes the
construction and its interpretation boundary. Numerical rank and
stabilizer-support outputs are intentionally not reported as validated
until the updated workflow has executed.


## Follow-up: SU(2) holonomy versus endpoint-coordinate quotient (2026-11-18)

A separate finite kinematic gate now distinguishes three maps that were
previously easy to conflate:

1. the six endpoint coordinates on an oriented link;
2. the three-dimensional tangent space of its single (SU(2)) holonomy;
3. the graph-wide vertex-gauge orbit in the product of link-holonomy spaces.

With (g_e\mapsto h_u g_e h_v^{-1}) and
(T_a=-i\sigma_a/2), the exact left-trivialized differential is
\[
\delta g_e g_e^{-1}=\xi_u-\operatorname{Ad}_{g_e}\xi_v.
\]
The six-to-three endpoint map is
\[
H_e=[I_3\; -\operatorname{Ad}_{g_e}],
\]
which has rank three and a three-dimensional kernel
(x_{e,u}=\operatorname{Ad}_{g_e}x_{e,v}). This kernel means that two
endpoint-coordinate descriptions can produce the same variation of one
group element. No independent edge-frame field has been introduced, so it
must not be relabelled as a new physical edge-frame gauge symmetry.

Composing this map with the vertex-incidence embedding gives the true
infinitesimal vertex-gauge Jacobian
\[
(A\xi)_e=\xi_u-\operatorname{Ad}_{g_e}\xi_v.
\]
The physical holonomy tangent quotient has dimension
(3|E|-\operatorname{rank}A), which depends on holonomies as well as graph
topology. On a connected tree all link holonomies are gauge removable. On a
generic single cycle, the remaining local modulus is the conjugacy class of
the ordered loop holonomy (one real parameter); the vertex-gauge stabilizer
has dimension one.

New files:
- `scripts/q2_relational_holonomy_quotient_gate.py`
- `tests/test_q2_relational_holonomy_quotient_gate.py`
- `BQG_Q2_RELATIONAL_HOLONOMY_QUOTIENT_RESULT_20261118.md`

The gate checks the adjoint representation, endpoint-map rank/kernel,
vertex-gauge composition, finite cycle-holonomy conjugation, and a
finite-difference check of the infinitesimal law. Its tests and JSON export
are now included in the graph workflow. CI execution is not yet verified.

This remains a kinematic SU(2) control. It does not yet identify the
projective boundary-state stabilizer with the holonomy stabilizer, nor derive
the physical QFI conductance or a continuum limit. The next step is to
compare the differential of the actual contraction map with the holonomy
quotient map, rather than identify their null spaces by assumption.


## Follow-up: holonomy-to-boundary tangent pullback (2026-11-18)

A new explicit linear map now connects the existing internal-edge probe
convention to the full half-link projective tangent. For each oriented
edge (e=(u,l_u,v,l_v)), the current graph-QFI convention inserts its
three generators at the first endpoint ((u,l_u)). The injection
(P:\mathbb R^{3E}\to\mathbb R^{12V}) selects exactly those coordinates.
For the projective boundary-state differential (D_\partial),
\[
D_e=D_\partial P,\qquad
F_e=P^T F_{\rm half}P
=4\mathcal R(D_e)^T\mathcal R(D_e).
\]
The new gate compares this pullback with the direct edge-QFI contraction
on the path and cycle fixtures and checks tangent/QFI rank agreement and
positive semidefiniteness. It reports the kinematic holonomy gauge-orbit
rank separately; it does not subtract it from the boundary tangent.

Files:
- `scripts/q2_relational_holonomy_boundary_tangent_gate.py`
- `tests/test_q2_relational_holonomy_boundary_tangent_gate.py`
- `BQG_Q2_RELATIONAL_HOLONOMY_BOUNDARY_TANGENT_RESULT_20261118.md`

The CI workflow now runs the focused tests and exports
`verification_results/Q2_RELATIONAL_HOLONOMY_BOUNDARY_TANGENT.json`.

**Scope remains deliberately narrow:** the pullback is exact for the
existing first-endpoint probe section. Endpoint-section independence and
the simultaneous vertex action on internal holonomies plus boundary frames
are not yet derived. Consequently the holonomy quotient and projective
boundary-state stabilizer remain distinct mathematical objects; no physical
coarse conductance or continuum result is claimed.


## Follow-up: two-sided endpoint QFI and singlet quotient (2026-11-18)

For the normalized invariant spin-1/2 singlet used on each internal edge,
the simultaneous fundamental action on its two tensor indices leaves the
gluing tensor invariant. In the current contraction convention this yields
the exact projective tangent identity
\[
D_v(e)=-D_u(e)
\]
for equal Pauli-axis parameters applied at the two ends of the same
contracted edge. This identity is independent of the other node intertwiners
and remains valid inside a multi-edge contraction.

Let (F\in\mathbb R^{3E\times3E}) denote the full correlated edge QFI,
including off-diagonal edge blocks. With endpoint coordinates interleaved
edge by edge, define (S:\mathbb R^{6E}\to\mathbb R^{3E}) by
((Sx)_e=x_{e,u}-x_{e,v}). The full two-sided metric is
\[
F_{\rm 2s}=S^TFS
=\begin{pmatrix}F&-F\\-F&F\end{pmatrix}_{e,e'}.
\]
It has at least (3E) structural common-endpoint null directions. Its
rank equals the rank of the relative edge QFI; additional nullity is exactly
the already diagnosed projective boundary-state stabilizer pulled back to
the relative coordinates.

New files:
- `scripts/q2_relational_two_sided_endpoint_qfi_gate.py`
- `tests/test_q2_relational_two_sided_endpoint_qfi_gate.py`
- `BQG_Q2_RELATIONAL_TWO_SIDED_ENDPOINT_QFI_RESULT_20261118.md`

The new CI step runs the focused tests and exports
`verification_results/Q2_RELATIONAL_TWO_SIDED_ENDPOINT_QFI.json`.
This is an exact result for the current normalized spin-1/2 singlet
contraction and tensor-index convention. It is not yet a derivation of a
general higher-spin edge-frame gauge theory or BQG coarse conductance.


## Follow-up: local singlet invariance for arbitrary SU(2) irreps (2026-11-18)

The two-sided endpoint cancellation extends beyond spin one-half at the
level of the invariant two-leg tensor. For (d=2j+1), in the magnetic
basis (m=j,j-1,\ldots,-j), define
\[
(C_j)_{mn}=\frac{(-1)^{j-m}}{\sqrt{2j+1}}\delta_{n,-m}.
\]
The representation generators obey
\[
J_a C_j+C_jJ_a^T=0,
\]
and hence (U_j C_j U_j^T=C_j). The tensor is antisymmetric for
half-integer spin and symmetric for integer spin.

New files:
- `scripts/q2_su2_irrep_singlet_invariance_gate.py`
- `tests/test_q2_su2_irrep_singlet_invariance_gate.py`
- `BQG_Q2_SU2_IRREP_SINGLET_INVARIANCE_RESULT_20261118.md`

The finite-dimensional gate constructs ladder-operator matrices and checks
Lie-algebra and finite-group invariance for (j=1/2,1,3/2,2,5/2,3);
tests additionally include (j=4). It is integrated into the workflow and
exports `verification_results/Q2_SU2_IRREP_SINGLET_INVARIANCE.json`.

This establishes representation-independent local singlet cancellation,
not higher-spin node intertwiners or network QFI. The latter remain the
next required step before a spin-summed BQG edge metric can be claimed.


## CI closure: two-sided endpoint quotient and arbitrary SU(2) irreps

The complete workflow run `44336897934` passed, including all tests and
exports for the graph, gauge, holonomy, projective-boundary, two-sided
endpoint, arbitrary-irrep, Feshbach, and PL generation-2 controls. The
previously contradictory Feshbach test assertion was removed; the positive
certified gap-bound assertions remain.

For six random intertwiner samples per graph, the two-sided endpoint QFI
had rank 6/nullity 6 on the two-edge path and rank 12/nullity 12 on the
four-edge cycle. Thus the only null directions in these edge-probe metrics
were the (3E) common endpoint rotations. The endpoint tangent-sum residual
was exactly zero in exported floating-point results; signed-lift QFI
reconstruction errors were below (1.51\times10^{-15}).

The arbitrary-irrep local singlet tests passed for the six spins
(j=1/2,1,3/2,2,5/2,3), with maximum finite-group invariance residual
(1.59\times10^{-15}); the test suite also covers (j=4).

The network-QFI JSON exporter was made NumPy-scalar safe, and a regression
test now checks its serialized result. All workflow export stages completed.


## Follow-up: higher-spin four-valent network QFI (2026-11-18)

The spin-one-half network calculation has now been generalized to individual
finite SU(2) irreducible representations (j=1/2,1,3/2,2). For each (j),
the four-valent intertwiner space is constructed as the joint kernel of the
three total node generators; its expected dimension is (2j+1). The
normalized invariant two-leg tensor is
\[
(C_j)_{mn}=\frac{(-1)^{j-m}}{\sqrt{2j+1}}\delta_{n,-m}.
\]

For a retained leg of a normalized invariant four-valent node, Schur's
lemma gives (ho_{m leg}=I_{2j+1}/(2j+1)), and therefore the isolated
single-link QFI is
\[
F^{(e)}_{ab}=\frac{4j(j+1)}{3}\delta_{ab}.
\]
For two links meeting at a middle node, the full path QFI has an isotropic
cross block
\[
F^{(12)}_{ab}=\frac43
\langle\iota_B|\mathbf J_{B0}\cdot\mathbf J_{B3}|\iota_B\rangle
\delta_{ab},
\]
with pair-Casimir spectrum
\[
\frac12[k(k+1)-2j(j+1)],\qquad k=0,1,\ldots,2j.
\]
The new gate checks this Casimir identity sample-by-sample, rather than
only checking that cross terms fall within spectral bounds.

New files:
- `scripts/q2_relational_higher_spin_network_qfi_gate.py`
- `tests/test_q2_relational_higher_spin_network_qfi_gate.py`
- `BQG_Q2_RELATIONAL_HIGHER_SPIN_NETWORK_QFI_RESULT_20261118.md`

The workflow runs the higher-spin tests and exports
`verification_results/Q2_RELATIONAL_HIGHER_SPIN_NETWORK_QFI.json`.
This is a finite pure-state path result for individual spins, not a spin
sum or a BQG history-weighted continuum conductance.
