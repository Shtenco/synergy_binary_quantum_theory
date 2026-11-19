#!/usr/bin/env python3
"""Compact direct generation-3 Q2 dual graph and causal-fibre gate.

This avoids the generation-3 Python incidence dictionary and per-edge JSON
objects. NumPy arrays hold child tetrahedra, fibre IDs, dual endpoints and
route colours. The existing generation-0..2 correspondence routines remain
the source of parent orientation, route colours and causal side assignments.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_pl_graph_export import (
    all_nonempty_faces,
    barycentric_subdivision,
    cross_polytope_boundary_4,
    dual_graph_payload,
)
from q2_causal_graph_export import ROUTE_LABELS
from q2_causal_pl_correspondence_gate import (
    _boundary_side,
    _causal_children,
    _parent_face_sides,
    _vertex_face_colours,
)


def _advance_one_generation(tetrahedra, parent_map, orientations,
                            causal_edges, next_vertex, generation,
                            chain_by_tet=None, previous_faces=None):
    dual = dual_graph_payload(
        generation, tetrahedra, chain_by_tet, previous_faces
    )
    vertex_colours = _vertex_face_colours(dual)
    face_sides, _ = _parent_face_sides(
        tetrahedra, parent_map, causal_edges, orientations, dual
    )
    children, faces, child_chains = barycentric_subdivision(tetrahedra)
    parent_index = {tuple(t): i for i, t in enumerate(tetrahedra)}
    child_causal, parent_by_child, next_vertex = _causal_children(
        causal_edges, generation + 1, next_vertex
    )
    child_map = {}
    child_orientations = {}
    for child in children:
        chain = child_chains[tuple(child)]
        parent = tuple(faces[chain[-1]])
        pidx = parent_index[parent]
        permutation = []
        previous = set()
        for face_id in chain:
            current = set(faces[face_id])
            added = current - previous
            if len(added) != 1:
                raise ValueError("invalid maximal barycentric flag")
            permutation.append(next(iter(added)))
            previous = current
        colour = vertex_colours[pidx][permutation[0]]
        parent_face = tuple(sorted(v for v in parent if v != permutation[-1]))
        side = face_sides[(parent, parent_face)]
        causal_id = f"g{generation + 1}:{parent_map[parent]}:{ROUTE_LABELS[colour]}:{side}"
        if causal_id not in child_causal:
            raise ValueError("child flag maps outside the causal child set")
        child_map[tuple(child)] = causal_id
        positions = [tuple(sorted(parent)).index(v) for v in permutation]
        inversions = sum(
            positions[i] > positions[j]
            for i in range(4) for j in range(i + 1, 4)
        )
        child_orientations[tuple(child)] = (
            orientations[parent] * (-1 if inversions % 2 else 1)
        )

    parent_fibres = Counter(parent_map[tuple(t)] for t in tetrahedra)
    child_fibres = Counter(child_map.values())
    child_mass_by_parent = Counter()
    for child_id, count in child_fibres.items():
        child_mass_by_parent[parent_by_child[child_id]] += count
    mass_ok = all(
        child_mass_by_parent[parent_id] == 24 * count
        for parent_id, count in parent_fibres.items()
    )
    if not mass_ok:
        raise ValueError("pushforward parent mass is not conserved")

    return (children, child_map, child_orientations, child_causal,
            next_vertex, faces, child_chains)


def _child_key(parent, permutation, face_id):
    prefix = []
    chain = []
    for vertex in permutation:
        prefix.append(vertex)
        chain.append(face_id[tuple(sorted(prefix))])
    return tuple(sorted(chain))


def _direct_dual(parent_tetrahedra, parent_dual, face_id):
    """Return compact (child_id_a, child_id_b, route) arrays."""
    permutations = [tuple(p) for p in itertools.permutations(range(4))]
    rank = {p: i for i, p in enumerate(permutations)}
    parent_index = {tuple(t): i for i, t in enumerate(parent_tetrahedra)}
    n_children = 24 * len(parent_tetrahedra)
    child_tets = np.empty((n_children, 4), dtype=np.uint32)

    for parent_id, tet in enumerate(parent_tetrahedra):
        for perm_rank, pidx in enumerate(permutations):
            p = tuple(tet[j] for j in pidx)
            child_id = 24 * parent_id + perm_rank
            child_tets[child_id] = _child_key(tet, p, face_id)

    endpoints = np.empty((3 * n_children // 2 + 6 * len(parent_dual["edges"]), 2),
                         dtype=np.uint32)
    routes = np.empty((len(endpoints),), dtype=np.uint8)
    cursor = 0

    # Three adjacent-transposition neighbours inside each parent simplex.
    for parent_id, tet in enumerate(parent_tetrahedra):
        for perm_rank, pidx in enumerate(permutations):
            for route in range(3):
                swapped = list(pidx)
                swapped[route], swapped[route + 1] = swapped[route + 1], swapped[route]
                other_rank = rank[tuple(swapped)]
                a = 24 * parent_id + perm_rank
                b = 24 * parent_id + other_rank
                if a < b:
                    endpoints[cursor] = (a, b)
                    routes[cursor] = route
                    cursor += 1

    # Six cross-parent neighbours for each paired parent triangular face.
    for edge in parent_dual["edges"]:
        ia, ib = edge["source"], edge["target"]
        ta, tb = tuple(parent_tetrahedra[ia]), tuple(parent_tetrahedra[ib])
        shared = tuple(sorted(edge["shared_face"]))
        oa = next(v for v in ta if v not in shared)
        ob = next(v for v in tb if v not in shared)
        for face_perm in itertools.permutations(shared):
            pa = tuple(ta.index(v) for v in (*face_perm, oa))
            pb = tuple(tb.index(v) for v in (*face_perm, ob))
            a = 24 * ia + rank[pa]
            b = 24 * ib + rank[pb]
            endpoints[cursor] = (min(a, b), max(a, b))
            routes[cursor] = 3
            cursor += 1

    if cursor != len(endpoints):
        raise AssertionError(f"edge allocation mismatch: {cursor} != {len(endpoints)}")
    if np.any(endpoints[:, 0] >= endpoints[:, 1]):
        raise AssertionError("dual edge endpoints are not canonical")
    return child_tets, endpoints, routes


def _compare_with_incidence(parent_tetrahedra, parent_dual):
    """Independent small-level comparison using the existing incidence exporter."""
    faces = all_nonempty_faces(parent_tetrahedra)
    face_id = {face: i for i, face in enumerate(faces)}
    permutations = list(itertools.permutations(range(4)))
    direct = set()
    for parent_id, tet in enumerate(parent_tetrahedra):
        for perm_rank, pidx in enumerate(permutations):
            p = tuple(tet[j] for j in pidx)
            key = _child_key(tet, p, face_id)
            for route in range(3):
                swapped = list(p)
                swapped[route], swapped[route + 1] = swapped[route + 1], swapped[route]
                other_key = _child_key(tet, tuple(swapped), face_id)
                if key < other_key:
                    direct.add((key, other_key, route))
    for edge in parent_dual["edges"]:
        ta = tuple(parent_tetrahedra[edge["source"]])
        tb = tuple(parent_tetrahedra[edge["target"]])
        shared = tuple(sorted(edge["shared_face"]))
        oa = next(v for v in ta if v not in shared)
        ob = next(v for v in tb if v not in shared)
        for p in itertools.permutations(shared):
            ka = _child_key(ta, (*p, oa), face_id)
            kb = _child_key(tb, (*p, ob), face_id)
            direct.add((min(ka, kb), max(ka, kb), 3))

    children, child_faces, chain_by_tet = barycentric_subdivision(parent_tetrahedra)
    incidence = dual_graph_payload(
        1, children, chain_by_tet, child_faces
    )
    expected = set()
    for edge in incidence["edges"]:
        ka = tuple(incidence["vertices"][edge["source"]]["tetrahedron"])
        kb = tuple(incidence["vertices"][edge["target"]]["tetrahedron"])
        expected.add((min(ka, kb), max(ka, kb), edge["route_colour"]))
    return {
        "direct_edge_signatures": len(direct),
        "incidence_edge_signatures": len(expected),
        "exact_edge_and_route_match": direct == expected,
        "missing_from_direct": len(expected - direct),
        "extra_in_direct": len(direct - expected),
    }


def run():
    tetrahedra = cross_polytope_boundary_4()
    parent_map = {tuple(t): "e0" for t in tetrahedra}
    orientations = {
        tuple(t): (-1 if sum(v % 2 for v in t) % 2 else 1)
        for t in tetrahedra
    }
    causal_edges = {"e0": (0, 1)}
    next_vertex = 2

    # Compare the direct construction against incidence matching at generations
    # 1 and 2 before using it for the larger generation-3 graph.
    comparison = []
    chain_by_tet = None
    previous_faces = None
    for generation in range(2):
        parent_dual = dual_graph_payload(
            generation, tetrahedra, chain_by_tet, previous_faces
        )
        comparison.append({
            "generation": generation + 1,
            **_compare_with_incidence(tetrahedra, parent_dual),
        })
        (tetrahedra, parent_map, orientations, causal_edges,
         next_vertex, previous_faces, chain_by_tet) = _advance_one_generation(
            tetrahedra, parent_map, orientations, causal_edges,
            next_vertex, generation, chain_by_tet, previous_faces
        )

    generation2_dual = dual_graph_payload(
        2, tetrahedra, chain_by_tet, previous_faces
    )
    if not generation2_dual["all_faces_two_sided_and_colour_consistent"]:
        raise ValueError("generation-2 parent dual face pairing is invalid")

    faces = all_nonempty_faces(tetrahedra)
    face_id = {face: i for i, face in enumerate(faces)}
    child_tets, endpoints, routes = _direct_dual(
        tetrahedra, generation2_dual, face_id
    )
    n_children = len(child_tets)
    degree = np.bincount(endpoints.ravel(), minlength=n_children)
    if not np.all(degree == 4):
        raise ValueError("generation-3 dual graph is not 4-regular")
    unique_edge_count = len(np.unique(endpoints, axis=0))
    if unique_edge_count != len(endpoints):
        raise ValueError("generation-3 direct construction emitted duplicate edges")
    route_degree = np.zeros((n_children, 4), dtype=np.uint8)
    for route in range(4):
        selected = endpoints[routes == route]
        np.add.at(route_degree[:, route], selected[:, 0], 1)
        np.add.at(route_degree[:, route], selected[:, 1], 1)
    if not np.all(route_degree == 1):
        raise ValueError("each child must have exactly one neighbour in each route")

    # Build the generation-3 causal fibres using the already-audited parent
    # route-colour and oriented face-side rules.
    colours = _vertex_face_colours(generation2_dual)
    face_sides, _ = _parent_face_sides(
        tetrahedra, parent_map, causal_edges, orientations, generation2_dual
    )
    child_causal, parent_by_child, _ = _causal_children(
        causal_edges, 3, next_vertex
    )
    causal_ids = sorted(child_causal)
    causal_index = {key: i for i, key in enumerate(causal_ids)}
    fibre_labels = np.empty(n_children, dtype=np.uint16)
    for parent_id, tet in enumerate(tetrahedra):
        for perm_rank, pidx in enumerate(itertools.permutations(range(4))):
            p = tuple(tet[j] for j in pidx)
            colour = colours[parent_id][p[0]]
            face = tuple(sorted(v for v in tet if v != p[-1]))
            side = face_sides[(tuple(tet), face)]
            causal_id = f"g3:{parent_map[tuple(tet)]}:{ROUTE_LABELS[colour]}:{side}"
            if causal_id not in causal_index:
                raise ValueError("generation-3 child references an unknown causal edge")
            fibre_labels[24 * parent_id + perm_rank] = causal_index[causal_id]

    fibre_sizes = np.bincount(fibre_labels, minlength=len(causal_ids))
    fibre_size_histogram = Counter(int(size) for size in fibre_sizes)
    if len(causal_ids) != 512 or np.any(fibre_sizes == 0):
        raise ValueError("generation-3 causal map is not surjective onto 512 fibres")

    # Each parent tetrahedron contributes 24 flags, and each parent causal
    # fibre splits into exactly eight child causal fibres in total.
    parent_counts = Counter(parent_map[tuple(t)] for t in tetrahedra)
    child_counts = Counter(causal_ids[int(i)] for i in fibre_labels)
    child_mass_by_parent = Counter()
    for child_id, count in child_counts.items():
        child_mass_by_parent[parent_by_child[child_id]] += count
    mass_ok = all(
        child_mass_by_parent[parent_id] == 24 * count
        for parent_id, count in parent_counts.items()
    )

    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components
    except ImportError as exc:
        raise RuntimeError("generation-3 connectivity gate requires scipy") from exc
    graph = csr_matrix(
        (np.ones(2 * len(endpoints), dtype=np.uint8),
         (np.concatenate((endpoints[:, 0], endpoints[:, 1])),
          np.concatenate((endpoints[:, 1], endpoints[:, 0])))),
        shape=(n_children, n_children),
    )
    component_count, _ = connected_components(graph, directed=False)

    checks = {
        "generation1_direct_matches_incidence": comparison[0]["exact_edge_and_route_match"],
        "generation2_direct_matches_incidence": comparison[1]["exact_edge_and_route_match"],
        "parent_face_pairing_valid": generation2_dual["all_faces_two_sided_and_colour_consistent"],
        "exact_child_tetrahedra": n_children == 221184,
        "exact_dual_edge_count": len(endpoints) == 442368,
        "within_parent_edge_count": int(np.sum(routes < 3)) == 331776,
        "across_parent_edge_count": int(np.sum(routes == 3)) == 110592,
        "four_regular": bool(np.all(degree == 4)),
        "unique_dual_edges": unique_edge_count == len(endpoints),
        "one_neighbour_per_route_at_every_child": bool(np.all(route_degree == 1)),
        "connected": component_count == 1,
        "causal_fibre_count_512": len(causal_ids) == 512,
        "all_causal_fibres_nonempty": bool(np.all(fibre_sizes > 0)),
        "pushforward_parent_mass_conserved": mass_ok,
    }
    return {
        "status": "Q2 generation-3 compact direct dual and causal-fibre gate",
        "passed": all(checks.values()),
        "counts": {
            "parent_tetrahedra": len(tetrahedra),
            "child_tetrahedra": n_children,
            "dual_edges": len(endpoints),
            "unique_dual_edges": unique_edge_count,
            "route_edge_counts": [int(np.sum(routes == k)) for k in range(4)],
            "within_parent_edges": int(np.sum(routes < 3)),
            "across_parent_edges": int(np.sum(routes == 3)),
            "causal_fibres": len(causal_ids),
            "fibre_size_min": int(fibre_sizes.min()),
            "fibre_size_max": int(fibre_sizes.max()),
            "fibre_size_mean": float(fibre_sizes.mean()),
            "fibre_size_histogram": {
                str(size): fibre_size_histogram[size]
                for size in sorted(fibre_size_histogram)
            },
            "pushforward_measure_total": float(fibre_sizes.sum() / (24 ** 3)),
            "connected_components": int(component_count),
        },
        "small_level_incidence_comparisons": comparison,
        "checks": checks,
        "storage": {
            "child_tetrahedra": str(child_tets.dtype) + " " + str(child_tets.shape),
            "dual_endpoints": str(endpoints.dtype) + " " + str(endpoints.shape),
            "route_colours": str(routes.dtype) + " " + str(routes.shape),
            "causal_fibre_labels": str(fibre_labels.dtype) + " " + str(fibre_labels.shape),
        },
        "scope_boundary": (
            "Finite combinatorial generation-3 graph and causal-fibre gate only; "
            "this does not identify PL adjacency with physical QFI conductance "
            "or establish a continuum spacetime limit."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
