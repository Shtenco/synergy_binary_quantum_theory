#!/usr/bin/env python3
"""Bounded-chunk generation-3 Q2 dual graph exporter.

Child IDs are stable: child_id = 24 * parent_tetrahedron_id + lexicographic
permutation rank. Large child/edge/fibre arrays are written to NumPy memmaps
in parent chunks. Validation is performed from the persisted arrays, not from
an in-memory copy of the generated graph.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import resource
from pathlib import Path

import numpy as np

from q2_pl_generation3_stream import (
    _advance_one_generation,
    _child_key,
)
from q2_pl_graph_export import (
    all_nonempty_faces,
    cross_polytope_boundary_4,
    dual_graph_payload,
)
from q2_causal_graph_export import ROUTE_LABELS
from q2_causal_pl_correspondence_gate import (
    _causal_children,
    _parent_face_sides,
    _vertex_face_colours,
)


PERMUTATIONS = tuple(itertools.permutations(range(4)))
PERMUTATION_RANK = {p: i for i, p in enumerate(PERMUTATIONS)}


def stable_child_id(parent_id: int, permutation_rank: int) -> int:
    if parent_id < 0 or not 0 <= permutation_rank < 24:
        raise ValueError("invalid stable child identity")
    return 24 * int(parent_id) + int(permutation_rank)


def _sha256_memmap(array, chunk_rows: int) -> str:
    digest = hashlib.sha256()
    for start in range(0, len(array), chunk_rows):
        digest.update(np.ascontiguousarray(array[start:start + chunk_rows]).tobytes())
    return digest.hexdigest()


def _parent_state():
    tetrahedra = cross_polytope_boundary_4()
    parent_map = {tuple(t): "e0" for t in tetrahedra}
    orientations = {
        tuple(t): (-1 if sum(v % 2 for v in t) % 2 else 1)
        for t in tetrahedra
    }
    causal_edges = {"e0": (0, 1)}
    next_vertex = 2
    for generation in range(2):
        (tetrahedra, parent_map, orientations, causal_edges,
         next_vertex, previous_faces, chain_by_tet) = _advance_one_generation(
            tetrahedra, parent_map, orientations, causal_edges,
            next_vertex, generation,
            None if generation == 0 else chain_by_tet,
            None if generation == 0 else previous_faces,
        )
    dual = dual_graph_payload(2, tetrahedra, chain_by_tet, previous_faces)
    if not dual["all_faces_two_sided_and_colour_consistent"]:
        raise ValueError("generation-2 parent dual face pairing is invalid")
    return tetrahedra, parent_map, orientations, causal_edges, next_vertex, dual


def export(output_dir: Path, *, chunk_size: int = 4096) -> dict:
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)

    parents, parent_map, orientations, causal_edges, next_vertex, parent_dual = _parent_state()
    n_parent = len(parents)
    n_child = 24 * n_parent
    n_within = 36 * n_parent
    n_cross = 6 * len(parent_dual["edges"])
    n_edges = n_within + n_cross
    if n_child != 221184 or n_edges != 442368:
        raise ValueError("unexpected generation-3 allocation dimensions")

    faces = all_nonempty_faces(parents)
    face_id = {face: i for i, face in enumerate(faces)}
    colours = _vertex_face_colours(parent_dual)
    face_sides, _ = _parent_face_sides(
        parents, parent_map, causal_edges, orientations, parent_dual
    )
    child_causal, parent_by_child, _ = _causal_children(
        causal_edges, 3, next_vertex
    )
    causal_ids = sorted(child_causal)
    causal_index = {key: i for i, key in enumerate(causal_ids)}
    if len(causal_ids) != 512:
        raise ValueError("generation-3 must contain exactly 512 causal fibres")

    child_tets = np.lib.format.open_memmap(
        output_dir / "g3_child_tetrahedra.npy", mode="w+", dtype=np.uint32,
        shape=(n_child, 4),
    )
    endpoints = np.lib.format.open_memmap(
        output_dir / "g3_dual_endpoints.npy", mode="w+", dtype=np.uint32,
        shape=(n_edges, 2),
    )
    routes = np.lib.format.open_memmap(
        output_dir / "g3_route_colours.npy", mode="w+", dtype=np.uint8,
        shape=(n_edges,),
    )
    fibre_labels = np.lib.format.open_memmap(
        output_dir / "g3_fibre_labels.npy", mode="w+", dtype=np.uint16,
        shape=(n_child,),
    )

    # Each parent owns exactly 24 child rows and 36 internal dual edges.
    for parent_start in range(0, n_parent, chunk_size):
        parent_stop = min(n_parent, parent_start + chunk_size)
        for parent_id in range(parent_start, parent_stop):
            tet = tuple(parents[parent_id])
            for rank, permutation in enumerate(PERMUTATIONS):
                ordered = tuple(tet[j] for j in permutation)
                child_id = stable_child_id(parent_id, rank)
                child_tets[child_id] = _child_key(tet, ordered, face_id)

                colour = colours[parent_id][ordered[0]]
                face = tuple(sorted(v for v in tet if v != ordered[-1]))
                side = face_sides[(tet, face)]
                causal_id = f"g3:{parent_map[tet]}:{ROUTE_LABELS[colour]}:{side}"
                if causal_id not in causal_index:
                    raise ValueError("child flag maps outside the causal fibre set")
                fibre_labels[child_id] = causal_index[causal_id]

            local_edge = 0
            for rank, permutation in enumerate(PERMUTATIONS):
                for route in range(3):
                    swapped = list(permutation)
                    swapped[route], swapped[route + 1] = swapped[route + 1], swapped[route]
                    other_rank = PERMUTATION_RANK[tuple(swapped)]
                    a = stable_child_id(parent_id, rank)
                    b = stable_child_id(parent_id, other_rank)
                    if a < b:
                        edge_id = 36 * parent_id + local_edge
                        endpoints[edge_id] = (a, b)
                        routes[edge_id] = route
                        local_edge += 1
            if local_edge != 36:
                raise AssertionError("parent did not emit exactly 36 within-parent edges")

    # Cross-parent neighbours have a stable suffix after all within-parent edges.
    for parent_edge_id, edge in enumerate(parent_dual["edges"]):
        ia, ib = edge["source"], edge["target"]
        ta, tb = tuple(parents[ia]), tuple(parents[ib])
        shared = tuple(sorted(edge["shared_face"]))
        oa = next(v for v in ta if v not in shared)
        ob = next(v for v in tb if v not in shared)
        for local_id, face_perm in enumerate(itertools.permutations(shared)):
            pa = tuple(ta.index(v) for v in (*face_perm, oa))
            pb = tuple(tb.index(v) for v in (*face_perm, ob))
            a = stable_child_id(ia, PERMUTATION_RANK[pa])
            b = stable_child_id(ib, PERMUTATION_RANK[pb])
            endpoints[n_within + 6 * parent_edge_id + local_id] = (min(a, b), max(a, b))
            routes[n_within + 6 * parent_edge_id + local_id] = 3

    for array in (child_tets, endpoints, routes, fibre_labels):
        array.flush()

    # Reconcile from the persisted memmaps in bounded edge chunks.
    degree = np.zeros(n_child, dtype=np.uint8)
    route_degree = np.zeros((n_child, 4), dtype=np.uint8)
    for start in range(0, n_edges, chunk_size * 36):
        stop = min(n_edges, start + chunk_size * 36)
        block = np.asarray(endpoints[start:stop])
        if np.any(block[:, 0] >= block[:, 1]):
            raise ValueError("non-canonical edge endpoint order")
        if np.any(block >= n_child):
            raise ValueError("edge endpoint outside generation-3 graph")
        for route in range(4):
            selected = block[np.asarray(routes[start:stop]) == route]
            if len(selected):
                np.add.at(route_degree[:, route], selected[:, 0], 1)
                np.add.at(route_degree[:, route], selected[:, 1], 1)
        np.add.at(degree, block[:, 0], 1)
        np.add.at(degree, block[:, 1], 1)

    canonical = np.sort(np.asarray(endpoints), axis=1)
    unique_edges = len(np.unique(canonical, axis=0))
    fibre_sizes = np.bincount(np.asarray(fibre_labels), minlength=len(causal_ids))
    mass_by_parent = {}
    for fibre_id, count in enumerate(fibre_sizes):
        parent_id = parent_by_child[causal_ids[fibre_id]]
        mass_by_parent[parent_id] = mass_by_parent.get(parent_id, 0) + int(count)
    expected_parent_mass = {
        parent_id: 24 * count
        for parent_id, count in __import__("collections").Counter(
            parent_map[tuple(t)] for t in parents
        ).items()
    }
    mass_ok = mass_by_parent == expected_parent_mass

    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components
    except ImportError as exc:
        raise RuntimeError("scipy is required for the connectivity reconciliation") from exc
    graph = csr_matrix(
        (np.ones(2 * n_edges, dtype=np.uint8),
         (np.concatenate((endpoints[:, 0], endpoints[:, 1])),
          np.concatenate((endpoints[:, 1], endpoints[:, 0])))),
        shape=(n_child, n_child),
    )
    component_count, _ = connected_components(graph, directed=False)

    checks = {
        "exact_child_count": n_child == 221184,
        "exact_edge_count": n_edges == 442368,
        "unique_edges": unique_edges == n_edges,
        "four_regular": bool(np.all(degree == 4)),
        "one_neighbour_per_route": bool(np.all(route_degree == 1)),
        "connected": int(component_count) == 1,
        "exact_fibre_count": len(causal_ids) == 512,
        "all_fibres_nonempty": bool(np.all(fibre_sizes > 0)),
        "pushforward_parent_mass_conserved": mass_ok,
    }
    if not all(checks.values()):
        raise ValueError(f"generation-3 chunked reconciliation failed: {checks}")

    files = {
        "child_tetrahedra": child_tets,
        "dual_endpoints": endpoints,
        "route_colours": routes,
        "fibre_labels": fibre_labels,
    }
    manifest = {
        "schema": "q2-g3-chunked-export-v1",
        "status": "PASS_G3_CHUNKED_EXPORT",
        "passed": True,
        "generation": 3,
        "parent_tetrahedra": n_parent,
        "child_tetrahedra": n_child,
        "dual_edges": n_edges,
        "within_parent_edges": n_within,
        "across_parent_edges": n_cross,
        "causal_fibres": len(causal_ids),
        "chunk_size_parent_tetrahedra": chunk_size,
        "stable_child_id_rule": "24 * parent_id + lexicographic_permutation_rank",
        "checks": checks,
        "measure_total": float(fibre_sizes.sum() / (24 ** 3)),
        "fibre_size_min": int(fibre_sizes.min()),
        "fibre_size_max": int(fibre_sizes.max()),
        "connected_components": int(component_count),
        "array_files": {
            name: {
                "file": filename,
                "dtype": str(array.dtype),
                "shape": list(array.shape),
                "sha256": _sha256_memmap(array, chunk_size * 36),
            }
            for name, filename, array in (
                ("child_tetrahedra", "g3_child_tetrahedra.npy", child_tets),
                ("dual_endpoints", "g3_dual_endpoints.npy", endpoints),
                ("route_colours", "g3_route_colours.npy", routes),
                ("fibre_labels", "g3_fibre_labels.npy", fibre_labels),
            )
        },
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "storage_mode": "numpy open_memmap; parent-chunk writes; validation reads bounded edge slices",
        "scope_boundary": "Finite combinatorial G3 export only; not physical QFI conductance or continuum closure.",
    }
    (output_dir / "G3_CHUNKED_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=4096)
    args = parser.parse_args()
    result = export(args.output_dir, chunk_size=args.chunk_size)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
