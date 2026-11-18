# Q2 finite causal-to-PL correspondence gate — result (2026-11-18)

## Result

A concrete, recursive, many-to-one map has been constructed from the tetrahedral vertices of the PL dual graph to active directed edges of the q=2 causal rewrite graph. The finite gate checks that every PL-dual adjacency maps to endpoint incidence of the corresponding causal edges, and that each refinement step maps a child tetrahedron to one of the eight causal children of its mapped parent edge.

This is stronger than a shared route alphabet. It is a finite graph correspondence between the PL tetrahedral dual graph and the line graph of the underlying causal active-edge graph (with source/target endpoint roles retained in the side rule).

It is **not** a graph isomorphism, graph covering, or a physical continuum map.

## Construction

At generation zero, all 16 PL tetrahedra map to the single causal active edge e0.

For each parent tetrahedron T and each of its four vertices v, the face T\{v} has a canonical route colour c(v) supplied by the existing PL shared-face gluing. In a barycentric child flag, let (v1,v2,v3,v4) be the ordered sequence of vertices added along the flag.

- The causal route label is the two-bit label associated with c(v1).
- The L/R side is attached to the parent face T\{v4}.
- If the two parent tetrahedra adjacent across that face map to distinct causal edges, those edges must share exactly one endpoint. The side is forced to select that shared endpoint: L selects the source of the parent edge, R its target.
- If both parent tetrahedra map to the same causal edge, the opposite oriented boundary signs of the two PL tetrahedra select opposite sides.

Thus the rule is recursive and uses only the existing PL face colours, barycentric flags, PL orientation, and causal edge incidence. It introduces no fitted lengths, angles, weights, or quantum states.

## Exact finite counts

| Generation | PL tetrahedra | Causal active edges | Mean fibre size | Fibre sizes |
|---:|---:|---:|---:|---|
| 0 | 16 | 1 | 16 | 16: 1 edge |
| 1 | 384 | 8 | 48 | 48: 8 edges |
| 2 | 9,216 | 64 | 144 | 64: 8; 112: 8; 128: 8; 144: 16; 160: 8; 176: 8; 224: 8 |

At every level the map is surjective and every fibre is nonempty. At generation 2 the fibre multiplicities are nonuniform: minimum 64, maximum 224, mean 144. This nonuniformity is a result, not something to smooth away.

## Incidence and refinement checks

The gate verifies:

1. every PL-dual vertex has one face of each of the four route colours;
2. every PL-dual edge joins two tetrahedra whose assigned causal edges share an endpoint;
3. distinct parent causal edges across a PL shared face are incident at exactly one endpoint;
4. equal parent causal edges use opposite oriented boundary sides;
5. each child tetrahedron maps to one of the eight causal child edges of its own mapped parent edge;
6. the map is surjective at generations 0, 1, and 2.

The PL-dual adjacency counts checked are 32, 768, and 18,432 respectively.

## What this does and does not prove

**Proved by the finite combinatorial construction:** a recursive, surjective, many-to-one adjacency-preserving correspondence exists through generation 2 for the frozen q=2 causal rewrite and the frozen barycentric PL refinement.

**Not proved:**
- extension to all generations or a coherent infinite inverse/direct limit;
- uniform bounded fibres (already false at generation 2);
- a compatible PL/casual vertex or incidence measure;
- a self-adjoint continuum operator or heat-kernel limit;
- a physical inter-block relational state and generator;
- QFI conductances, spectral dimension, semiclassical dynamics, or Einstein equations.

The nonuniform fibre distribution is now the central measure-theoretic issue. Any continuum/QFI construction must derive how these fibres contribute to the measure and operator. Assigning inverse-fibre weights by hand would be an additional model assumption, not a result of this gate.

## Reproduction

- `scripts/q2_causal_pl_correspondence_gate.py`
- `tests/test_q2_causal_pl_correspondence_gate.py`
- CI workflow: `.github/workflows/q2-pl-graph-export.yml`

The workflow is configured to run both the PL exporter tests and this correspondence gate, and to upload the JSON result. The code has been committed; CI execution must be checked independently before describing the tests as passing.
