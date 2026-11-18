# Canonical Q2 causal-route graph export — audit result (2026-11-18)

## Implementation

Added `scripts/q2_causal_graph_export.py`, a directed graph realization of
the frozen scalar route rewrite in `scripts/q2_dimension3_fixed_point_gate.py`.

For every active directed edge (u\to v), the rule creates four route vertices
labelled (00,01,10,11) and replaces the parent edge by eight child edges:

[
u\to m_a\to v,\qquad a\in\{00,01,10,11\}.
]

Only current-generation active edges are present in the level-g graph; all
previously created vertices remain as endpoints of the refined edges. The
exporter stores stable parent-edge IDs, route labels, child-edge IDs, graph
hashes, connectivity, DAG status, source-to-sink path counts and degree data.

This construction follows the exact scalar recurrence. It deliberately does
not add the C4 equatorial links between the four route states: those are
spatial/carrier incidence data in the separately selected PL completion and
cannot be silently conflated with directed causal edges.

## Verified finite levels

GitHub Actions run `44270543046` completed successfully:
- all 3 tests passed;
- export through generation 5 succeeded;
- artifact `q2-causal-graph-export` uploaded (635,213 bytes).

| Generation | Vertices | Active directed edges | Causal path length | Source-to-sink paths |
|---:|---:|---:|---:|---:|
| 0 | 2 | 1 | 1 | 1 |
| 1 | 6 | 8 | 2 | 4 |
| 2 | 38 | 64 | 4 | 64 |
| 3 | 294 | 512 | 8 | 16,384 |
| 4 | 2,342 | 4,096 | 16 | 1,073,741,824 |
| 5 | 18,726 | 32,768 | 32 | 4,611,686,018,427,387,904 |

The exact vertex count is
[
N_g=2+4\frac{8^g-1}{7}=\frac{4\cdot8^g+10}{7}.
]
The path count follows a different recurrence. Since every old path has
(2^g) active edges at generation g, refining every edge into four choices
gives
[
P_g=4^{1+2+\cdots+2^{g-1}}=4^{2^g-1}.
]
This is why the number of paths is not merely (4^g).

## Newly exposed structural issue: unbounded degree

The graph export exposes a consequence not visible in the scalar vertex
recurrence. Every old vertex incident to an active edge receives four child
edges in place of that edge. Its degree therefore multiplies by four at each
subsequent refinement. In particular, the source and sink degrees are
[
\deg_g(s)=\deg_g(t)=4^g.
]

The finite-level maximum degree is therefore 1, 4, 16, 64, 256, 1024 at
generations 0 through 5. In the literal direct limit, each retained vertex
has unbounded degree. Thus the raw rewrite is **not uniformly locally finite**.

This is a real continuum-limit blocker, not a coding defect. Before using its
ordinary graph Laplacian as a physical spatial operator, the theory must derive
one of the following from the microscopic dynamics:
1. a nontrivial vertex measure and self-adjoint weighted Laplacian with a
   controlled continuum limit;
2. a normalized edge/vertex refinement rule that preserves finite local
   measure and is compatible with the causal-volume theorem; or
3. a different, explicitly specified gluing/rewrite in which the route
   multiplicity is not realized as four new neighbors at every retained
   endpoint.

Any such modification changes the operator and may change the spectral
dimension. It must be tested rather than assumed harmless.

## Relation to the selected PL completion

The PL completion has a different finite seed: the boundary of the 4D
cross-polytope, with 8 vertices, 24 edges and 16 tetrahedra. Its tetrahedral
dual graph is (Q_4), a 16-vertex 4-regular graph. The scalar causal rewrite
starts from one directed edge and produces a hierarchical DAG with two
boundary vertices and rapidly increasing endpoint degree.

The two structures share the four route labels locally, but their global
incidence data are not identified by the current microscopic rule. A bijection
or incidence-preserving map cannot be inferred from the common label set or
from the matching asymptotic exponent. A separate dynamical gluing map is
still required.

## Decision

- **Closed:** deterministic finite-level directed causal graph realization of
  the frozen q=2 route-counting rule, with ancestry and graph identity.
- **Closed:** exact finite-level checks of vertex count, edge count, path
  multiplicity, connectivity, acyclicity and degree growth.
- **New blocker:** literal graph limit has unbounded vertex degree.
- **Open:** derived measure/normalization, causal-to-PL gluing, physical
  (ρ_e,K_e) assignments, QFI conductances and spectral-dimension plateau.

The exact causal-volume exponent remains (d_H=3) for the frozen scalar
recurrence. The new graph export does not establish a three-dimensional
spectral dimension, and the unbounded-degree result makes that next test more
subtle.
