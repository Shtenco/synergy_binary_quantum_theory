#!/usr/bin/env python3
"""Export the canonical q=2 PL-completion 1-skeleton at successive refinements.

The source complex is the boundary of the 4D cross-polytope (16-cell), exactly
as constructed by bcqg_global_manifold_gate.py. Refinement is global
barycentric subdivision. At level n>0, vertices are the non-empty simplices
of level n-1 and graph edges join comparable simplices (the 1-skeleton of the
order complex). No coordinates or target dimension enter this construction.

IMPORTANT: this exports the PL-completion graph, not the distinct causal
edge-rewrite graph used in the N_g recurrence. It also does not invent local
quantum states or QFI edge weights.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path



def cross_polytope_boundary_4():
    return sorted({
        tuple(sorted(2 * axis + signs[axis] for axis in range(4)))
        for signs in itertools.product((0, 1), repeat=4)
    })


def all_nonempty_faces(tetrahedra):
    faces = set()
    for tet in tetrahedra:
        for size in range(1, 5):
            faces.update(tuple(sorted(face)) for face in itertools.combinations(tet, size))
    return sorted(faces, key=lambda face: (len(face), face))


def barycentric_subdivision(tetrahedra):
    faces = all_nonempty_faces(tetrahedra)
    face_id = {face: index for index, face in enumerate(faces)}
    refined = set()
    for tet in tetrahedra:
        for permutation in itertools.permutations(tet):
            prefix = []
            chain = []
            for vertex in permutation:
                prefix.append(vertex)
                chain.append(face_id[tuple(sorted(prefix))])
            refined.add(tuple(sorted(chain)))
    return sorted(refined), faces


def graph_edges_from_chains(tetrahedra, faces):
    """Return all comparable face pairs, using tetrahedra as local cliques."""
    face_id = {face: index for index, face in enumerate(faces)}
    edges = set()
    for tet in tetrahedra:
        local_faces = [
            tuple(sorted(face))
            for size in range(1, 5)
            for face in itertools.combinations(tet, size)
        ]
        for a, b in itertools.combinations(local_faces, 2):
            if len(a) < len(b) and set(a).issubset(b):
                edges.add((face_id[a], face_id[b]))
            elif len(b) < len(a) and set(b).issubset(a):
                edges.add((face_id[b], face_id[a]))
    return sorted(edges)


def graph_payload(generation, tetrahedra, previous_faces=None, parent_tetrahedra=None):
    if generation == 0:
        vertices = sorted({v for tet in tetrahedra for v in tet})
        edges = sorted({
            tuple(sorted(edge))
            for tet in tetrahedra
            for edge in itertools.combinations(tet, 2)
        })
        node_labels = [{"id": v, "simplex": [v]} for v in vertices]
        refinement_map = None
    else:
        assert previous_faces is not None
        vertices = list(range(len(previous_faces)))
        if parent_tetrahedra is None:
            raise ValueError("parent_tetrahedra required for a refined graph")
        edges = graph_edges_from_chains(parent_tetrahedra, previous_faces)
        node_labels = [
            {"id": i, "simplex": list(face)}
            for i, face in enumerate(previous_faces)
        ]
        refinement_map = [
            {"child_vertex": i, "parent_simplex": list(face)}
            for i, face in enumerate(previous_faces)
        ]

    canonical = json.dumps(
        {"vertices": vertices, "edges": edges},
        sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    degrees = [0] * len(vertices)
    for a, b in edges:
        degrees[a] += 1
        degrees[b] += 1
    return {
        "generation": generation,
        "vertex_count": len(vertices),
        "edge_count": len(edges),
        "component_count": component_count(len(vertices), edges),
        "degree_min": min(degrees) if degrees else 0,
        "degree_max": max(degrees) if degrees else 0,
        "graph_sha256": hashlib.sha256(canonical).hexdigest(),
        "vertices": node_labels,
        "edges": [list(edge) for edge in edges],
        "refinement_map": refinement_map,
        "quantum_state_assignment": None,
        "qfi_edge_conductances": None,
    }


def component_count(n_vertices, edges):
    adjacency = [[] for _ in range(n_vertices)]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    seen = set()
    count = 0
    for root in range(n_vertices):
        if root in seen:
            continue
        count += 1
        stack = [root]
        seen.add(root)
        while stack:
            u = stack.pop()
            for v in adjacency[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    return count


def run(max_generation=2):
    if not 0 <= max_generation <= 2:
        raise ValueError("max_generation must be in [0, 2] for the audited exact PL levels")
    tets = cross_polytope_boundary_4()
    rows = []
    # Generation zero is the ordinary 1-skeleton of the 16-cell boundary.
    rows.append(graph_payload(0, tets))
    for generation in range(1, max_generation + 1):
        parent_tets = tets
        tets, parent_faces = barycentric_subdivision(parent_tets)
        rows.append(graph_payload(generation, tets, parent_faces, parent_tets))
    expected = [(8, 24), (80, 464), (1696, 10912)]
    checks = {
        "counts_match_independent_pl_gate": all(
            (row["vertex_count"], row["edge_count"]) == expected[i]
            for i, row in enumerate(rows)
        ),
        "all_graphs_connected": all(row["component_count"] == 1 for row in rows),
        "stable_hashes_present": all(len(row["graph_sha256"]) == 64 for row in rows),
        "no_fabricated_quantum_data": all(
            row["quantum_state_assignment"] is None
            and row["qfi_edge_conductances"] is None for row in rows
        ),
    }
    return {
        "status": "q2 canonical PL-completion graph export",
        "passed": all(checks.values()),
        "source_complex": "boundary of the 4D cross-polytope (16-cell)",
        "refinement": "global barycentric subdivision",
        "graph_semantics": (
            "1-skeleton of the PL complex at each level; this is not the separate "
            "causal edge-rewrite graph underlying the N_g recurrence"
        ),
        "generations": rows,
        "checks": checks,
        "qfi_readiness": "topology exported; physical reduced states and edge generators are not supplied by the current PL gate",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-generation", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.max_generation)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
