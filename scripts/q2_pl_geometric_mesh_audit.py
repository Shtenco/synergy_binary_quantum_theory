#!/usr/bin/env python3
"""Build a reference-PL geometric audit and weighted dual finite-volume operator.

The embedding is a declared Euclidean reference model (the boundary of the
4D cross-polytope in R^4), not an assertion about physical SBQT spacetime.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

try:
    from .q2_pl_graph_export import barycentric_subdivision, cross_polytope_boundary_4
    from .q2_pl_generation3_chunked_export import export as export_g3_chunked
    from .q2_dirichlet_form_convergence import GraphDirichletForm
except ImportError:
    from q2_pl_graph_export import barycentric_subdivision, cross_polytope_boundary_4
    from q2_pl_generation3_chunked_export import export as export_g3_chunked
    from q2_dirichlet_form_convergence import GraphDirichletForm


def initial_reference_embedding():
    """Embed cross-polytope vertices as the eight points +/- e_i in R^4."""
    coordinates = np.zeros((8, 4), dtype=float)
    for axis in range(4):
        coordinates[2 * axis, axis] = -1.0
        coordinates[2 * axis + 1, axis] = 1.0
    tetrahedra = cross_polytope_boundary_4()
    return coordinates, tetrahedra


def mesh_metrics(coordinates, tetrahedra):
    tets = np.asarray(tetrahedra, dtype=np.int64)
    if tets.ndim != 2 or tets.shape[1] != 4:
        raise ValueError("tetrahedra must have shape (n, 4)")
    points = np.asarray(coordinates, dtype=float)[tets]
    edge_vectors = points[:, 1:, :] - points[:, :1, :]
    gram = np.einsum("nki,nli->nkl", edge_vectors, edge_vectors)
    determinants = np.linalg.det(gram)
    if np.any(determinants <= 0) or not np.all(np.isfinite(determinants)):
        raise ValueError("degenerate or non-finite PL tetrahedron")
    volumes = np.sqrt(determinants) / 6.0
    pair_lengths = [
        np.linalg.norm(points[:, i, :] - points[:, j, :], axis=1)
        for i in range(4) for j in range(i + 1, 4)
    ]
    diameters = np.maximum.reduce(pair_lengths)
    return {
        "tetrahedra": int(len(tets)),
        "vertices": int(len(coordinates)),
        "pl_volume_total": float(np.sum(volumes)),
        "cell_volume_min": float(np.min(volumes)),
        "cell_volume_max": float(np.max(volumes)),
        "cell_volume_uniformity_relative": float(
            (np.max(volumes) - np.min(volumes)) / np.mean(volumes)
        ),
        "mesh_diameter_max": float(np.max(diameters)),
        "mesh_diameter_min": float(np.min(diameters)),
    }, volumes


def refine_embedding(coordinates, tetrahedra):
    children, faces, _chain_by_tet = barycentric_subdivision(tetrahedra)
    old = np.asarray(coordinates, dtype=float)
    new_coordinates = np.empty((len(faces), old.shape[1]), dtype=float)
    for face_id, face in enumerate(faces):
        new_coordinates[face_id] = np.mean(old[np.asarray(face, dtype=np.int64)], axis=0)
    return new_coordinates, children


def geometric_dual_weights(coordinates, tetrahedra, edge_endpoints, cell_volumes,
                            *, chunk_rows=50000):
    """Conductance = shared-face area / barycentric dual path length.

    The centroid-to-face-centroid distance inside a tetrahedron is one quarter
    of its altitude, so d_dual = (altitude_left + altitude_right)/4.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    tetrahedra = np.asarray(tetrahedra, dtype=np.int64)
    endpoints = np.asarray(edge_endpoints, dtype=np.int64)
    volumes = np.asarray(cell_volumes, dtype=float)
    weights = np.empty(len(endpoints), dtype=float)
    face_areas = np.empty(len(endpoints), dtype=float)
    dual_distances = np.empty(len(endpoints), dtype=float)

    for start in range(0, len(endpoints), chunk_rows):
        stop = min(len(endpoints), start + chunk_rows)
        block = endpoints[start:stop]
        left = tetrahedra[block[:, 0]]
        right = tetrahedra[block[:, 1]]
        matches = left[:, :, None] == right[:, None, :]
        if not np.all(matches.sum(axis=(1, 2)) == 3):
            raise ValueError("each dual edge must correspond to exactly one shared face")
        shared_mask = np.any(matches, axis=2)
        shared = np.sort(np.where(shared_mask, left, -1), axis=1)[:, -3:]
        triangles = coordinates[shared]
        v1 = triangles[:, 1] - triangles[:, 0]
        v2 = triangles[:, 2] - triangles[:, 0]
        g11 = np.einsum("ni,ni->n", v1, v1)
        g22 = np.einsum("ni,ni->n", v2, v2)
        g12 = np.einsum("ni,ni->n", v1, v2)
        area = 0.5 * np.sqrt(np.maximum(g11 * g22 - g12 * g12, 0.0))
        if np.any(area <= 0) or not np.all(np.isfinite(area)):
            raise ValueError("degenerate shared triangular face")
        altitude_left = 3.0 * volumes[block[:, 0]] / area
        altitude_right = 3.0 * volumes[block[:, 1]] / area
        distance = 0.25 * (altitude_left + altitude_right)
        conductance = area / distance
        if np.any(distance <= 0) or np.any(conductance <= 0):
            raise ValueError("non-positive dual distance or conductance")
        face_areas[start:stop] = area
        dual_distances[start:stop] = distance
        weights[start:stop] = conductance

    return weights, face_areas, dual_distances


def _sha256(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def build_geometry(max_generation=3):
    coordinates, tetrahedra = initial_reference_embedding()
    levels = []
    for generation in range(max_generation + 1):
        metrics, volumes = mesh_metrics(coordinates, tetrahedra)
        metrics.update({
            "generation": generation,
            "combinatorial_count_measure": float(len(tetrahedra) / (24 ** generation)),
            "expected_pl_volume": 16.0 / 3.0,
            "pl_volume_error": abs(metrics["pl_volume_total"] - 16.0 / 3.0),
        })
        if metrics["pl_volume_error"] > 1e-10:
            raise ValueError(f"PL volume is not conserved at generation {generation}")
        levels.append(metrics)
        if generation < max_generation:
            coordinates, tetrahedra = refine_embedding(coordinates, tetrahedra)
    for previous, current in zip(levels, levels[1:]):
        if current["mesh_diameter_max"] >= previous["mesh_diameter_max"]:
            raise ValueError("maximum PL mesh diameter did not decrease under refinement")
    return coordinates, np.asarray(tetrahedra, dtype=np.uint32), volumes, levels


def run(output_dir: Path, *, chunk_size=4096, chunk_rows=50000):
    output_dir.mkdir(parents=True, exist_ok=True)
    coordinates, tetrahedra, volumes, levels = build_geometry(3)
    chunk_dir = output_dir / "g3_chunked"
    chunk_manifest = export_g3_chunked(chunk_dir, chunk_size=chunk_size)
    if not chunk_manifest.get("passed"):
        raise ValueError("G3 combinatorial exporter did not pass")
    persisted_tets = np.load(chunk_dir / "g3_child_tetrahedra.npy", mmap_mode="r")
    endpoints = np.load(chunk_dir / "g3_dual_endpoints.npy", mmap_mode="r")
    if not np.array_equal(tetrahedra, persisted_tets):
        raise ValueError("reference PL tetrahedron ordering differs from G3 stable export")
    if tetrahedra.shape != (221184, 4) or endpoints.shape != (442368, 2):
        raise ValueError("unexpected G3 geometric audit dimensions")

    weights, face_areas, dual_distances = geometric_dual_weights(
        coordinates, tetrahedra, endpoints, volumes, chunk_rows=chunk_rows
    )
    form = GraphDirichletForm.build(
        len(tetrahedra), endpoints, weights, volumes, scale=1.0
    )
    degree = np.bincount(
        np.concatenate((endpoints[:, 0], endpoints[:, 1])),
        minlength=len(tetrahedra),
    )
    if not np.all(degree == 4):
        raise ValueError("geometric dual graph is not four-regular")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("geometric dual conductances must be finite and positive")

    arrays = {
        "g3_pl_vertex_coordinates.npy": coordinates,
        "g3_pl_tetrahedra.npy": tetrahedra,
        "g3_pl_cell_volumes.npy": volumes,
        "g3_geometric_dual_endpoints.npy": np.asarray(endpoints),
        "g3_geometric_dual_conductances.npy": weights,
        "g3_geometric_dual_face_areas.npy": face_areas,
        "g3_geometric_dual_distances.npy": dual_distances,
    }
    for filename, array in arrays.items():
        np.save(output_dir / filename, array, allow_pickle=False)

    manifest = {
        "schema": "q2-pl-geometric-mesh-audit-v1",
        "status": "PASS_G3_PL_GEOMETRIC_MESH_AUDIT",
        "passed": True,
        "reference_embedding": "cross-polytope boundary vertices +/-e_i in R^4; all refined vertices are barycentres of faces",
        "levels": levels,
        "g3": {
            "vertices_as_dual_cells": int(len(tetrahedra)),
            "dual_edges": int(len(endpoints)),
            "pl_vertex_count": int(len(coordinates)),
            "pl_volume_total": float(np.sum(volumes)),
            "combinatorial_count_measure": float(len(tetrahedra) / (24 ** 3)),
            "combinatorial_to_pl_volume_ratio": float((len(tetrahedra) / (24 ** 3)) / np.sum(volumes)),
            "degree_min": int(degree.min()),
            "degree_max": int(degree.max()),
            "conductance_min": float(weights.min()),
            "conductance_max": float(weights.max()),
            "shared_face_area_min": float(face_areas.min()),
            "shared_face_area_max": float(face_areas.max()),
            "dual_distance_min": float(dual_distances.min()),
            "dual_distance_max": float(dual_distances.max()),
            "finite_volume_form": form.summary(),
        },
        "arrays": {
            filename: {"shape": list(array.shape), "dtype": str(array.dtype),
                       "sha256": _sha256(array)}
            for filename, array in arrays.items()
        },
        "scope_boundary": (
            "Reference piecewise-flat geometry and a barycentric-dual finite-volume "
            "candidate. This is not asserted to be the physical SBQT metric and "
            "does not establish convergence to a continuum operator."
        ),
    }
    (output_dir / "Q2_PL_GEOMETRIC_MESH_AUDIT.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--chunk-rows", type=int, default=50000)
    args = parser.parse_args()
    result = run(args.output_dir, chunk_size=args.chunk_size,
                 chunk_rows=args.chunk_rows)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
