# Q2 causal-route ↔ PL bridge: structural audit (2026-11-18)

**Scope:** exact comparison of the existing scalar causal-route exporter and the canonical 16-cell PL primal/dual exporters. This note is a structural audit, not a proposed identification.

## 1. Source objects (keep distinct)

The causal exporter `scripts/q2_causal_graph_export.py` implements the frozen rule
`u -> v` becomes four labelled alternatives `u -> m_a -> v`, for
`a ∈ {00,01,10,11}`. It retains old vertices and only the newest active
directed edges. At generation g it reports

- `V_c(g) = 2 + 4(8^g - 1)/7`;
- `E_c,active(g) = 8^g`;
- source-to-sink path length `2^g`;
- path count `4^(2^g - 1)`.

The PL exporter `scripts/q2_pl_graph_export.py` starts instead from the
boundary of the 4D cross-polytope (16-cell), and performs global barycentric
subdivision. It exports both the primal 1-skeleton and tetrahedral dual graph.

| PL generation | primal (V,E) | dual (V,E) |
|---:|---:|---:|
| 0 | (8, 24) | (16, 32) |
| 1 | (80, 464) | (384, 768) |
| 2 | (1696, 10912) | (9216, 18432) |

The source seed, edge semantics, and refinement operations are different.

## 2. Exact obstruction to a levelwise graph isomorphism

A graph isomorphism preserves vertex count, edge count, degree multiset,
connectivity, and (for directed graphs) orientation data.

Already at generation zero:
- causal route graph: 2 vertices and 1 directed active edge;
- PL primal graph: 8 vertices and 24 undirected edges;
- PL dual graph: 16 vertices and 32 undirected edges.

Thus there is no level-zero graph isomorphism to either PL graph. At later
listed levels the vertex/edge counts also disagree. A common four-label
alphabet does not change this invariant.

This rules out a *canonical levelwise graph isomorphism* under the currently
frozen definitions. It does not rule out a non-bijective quotient, incidence
relation, span/correspondence, or a map after changing the microscopic rule.
Any such construction is new mathematical data and needs its own axioms and
tests.

## 3. Local-finiteness obstruction in the literal causal direct limit

The causal exporter explicitly retains old endpoints while every active edge
incident to an endpoint is replaced by four alternatives. Consequently source
and sink degrees are both `4^g`; the maximum degree is also `4^g`.

Therefore
`sup_g degree_g(source) = ∞`. The direct-limit graph is not uniformly
locally finite. A standard bounded-degree graph Laplacian/heat-kernel limit
cannot be assumed from the finite-generation adjacency alone.

A continuum construction must derive, rather than choose ad hoc:
1. a vertex/edge measure and Hilbert space;
2. a normalization or renormalized transition operator compatible with refinement;
3. a self-adjoint operator/domain and a controlled limit;
4. compatibility of that normalization with the causal-volume law.

Changing conductances by a generation-dependent scalar only rescales diffusion
time on a fixed finite graph. It does not cure unbounded incidence or establish
a dimension-three spectral plateau.

## 4. What a legitimate causal-to-PL bridge must specify

A proposed bridge must give, for every generation:

- a source object (causal vertices, active edges, and ancestry);
- a target object (PL simplices, tetrahedra, and face incidences);
- an explicit relation/map, including its domain and codomain;
- treatment of four route labels and the PL shared-face colours;
- preservation claims stated precisely (order, incidence, orientation,
  boundary, or quotient fibres);
- compatibility diagrams between successive refinement maps;
- invariants checked by code, including fibre sizes and local degree bounds.

A precise *local* correspondence is now established in code: at every exported
PL-dual vertex, its four incident dual edges have route colours exactly
{00, 01, 10, 11}, each once. The PL exporter checks this at generations 0–2.
This identifies the four alternatives of a q=2 rewrite cell with the four
coloured face-directions around one tetrahedral dual vertex. It is a local
star correspondence, not a map of the full causal graph: it does not assign
causal vertices/edges to PL tetrahedra/shared faces, nor prove refinement
naturality between the two constructions.

The scalar causal DAG starts from one directed edge, while the PL carrier
starts from a closed 3-manifold triangulation.

## 5. QFI consequence

No physical edge conductance follows from either topology export. The Walsh
local-QFI control gives an isotropic local frame and a reduced product-state
probe value `F_Q/4 = 1/3`, but after Gauss-singlet projection the global
SU(2) gauge-orbit QFI vanishes. The value 1/3 is not an inter-block physical
edge weight.

Until a shared-face relational state `rho_e` and physical relational
generator `K_e` are derived, a QFI Laplacian on the PL dual graph is a
placeholder, not a BQG prediction.

## 6. Updated decision: finite global correspondence found

**Decision:** do not identify the causal DAG with the PL primal or dual graph.
However, a concrete many-to-one correspondence from PL tetrahedra to active
causal edges has now been constructed and implemented in
`scripts/q2_causal_pl_correspondence_gate.py`.

The map uses:
- the route colour of the parent face opposite the first vertex in a
  barycentric flag to select one of the four causal route labels;
- the parent face opposite the final vertex to select the L/R child;
- causal endpoint incidence when the two parent tetrahedra map to distinct
  causal edges;
- opposite oriented PL boundary sides when both parent tetrahedra map to the
  same causal edge.

Through generations 0–2, the gate checks surjectivity, refinement naturality,
and preservation of PL-dual adjacency as endpoint incidence of causal edges.
The map is not an isomorphism or a graph covering. Its fibres are nonuniform
at generation 2: 64–224 PL tetrahedra per causal edge, with mean 144. The
exact histogram and construction are documented in
`BQG_Q2_CAUSAL_PL_CORRESPONDENCE_GATE_RESULT_20261118.md`.

This closes the *finite combinatorial correspondence* gate through generation
2, conditional on the configured CI tests. It does not close the infinite
refinement limit, the measure/operator problem, or physical QFI. The next gate
must derive a refinement-compatible measure from the nonuniform fibres and
test whether a corresponding weighted Laplacian has a controlled limit.
Do not insert inverse-fibre weights as if they were already physically
derived.

## Reproduction sources

- `scripts/q2_causal_graph_export.py`
- `tests/test_q2_causal_graph_export.py`
- `scripts/q2_pl_graph_export.py`
- `tests/test_q2_pl_graph_export.py`
- `BQG_QFI_Q2_EMERGENT_GRAPH_AUDIT_20261118.md`
- `BQG_Q2_WALSH_LOCAL_QFI_CONTROL_RESULT_20261118.md`

This audit makes no claim that spacetime or general relativity has been
derived. It narrows the exact missing theorem and prevents accidental
identification of distinct graph families.
