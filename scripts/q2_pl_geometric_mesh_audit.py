#!/usr/bin/env python3
"""Build a reference-PL geometric audit and weighted dual finite-volume operator.

The embedding is a declared Euclidean reference model (the boundary of the
4D cross-polytope in R^4), not an assertion about physical SBQT spacetime.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

try:
    from .q2_pl_graph_export import barycentric_subdivision, cross_polytope_boundary_4, dual_graph_payload
    from .q2_pl_generation3_stream import _advance_one_generation
    from .q2_pl_generation3_chunked_export import export as export_g3_chunked
    from .q2_dirichlet_form_convergence import GraphDirichletForm
except ImportError:
    from q2_pl_graph_export import barycentric_subdivision, cross_polytope_boundary_4, dual_graph_payload
    from q2_pl_generation3_stream import _advance_one_generation
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
    edge_vectors_all = [
        points[:, i, :] - points[:, j, :]
        for i in range(4) for j in range(i + 1, 4)
    ]
    edge_lengths = [np.linalg.norm(vector, axis=1) for vector in edge_vectors_all]
    edge_squared_sum = np.sum(np.stack([length * length for length in edge_lengths]), axis=0)
    diameters = np.maximum.reduce(edge_lengths)
    face_areas = []
    for omitted in range(4):
        triangle = np.delete(points, omitted, axis=1)
        first = triangle[:, 1, :] - triangle[:, 0, :]
        second = triangle[:, 2, :] - triangle[:, 0, :]
        g11 = np.einsum("ni,ni->n", first, first)
        g22 = np.einsum("ni,ni->n", second, second)
        g12 = np.einsum("ni,ni->n", first, second)
        face_areas.append(0.5 * np.sqrt(np.maximum(g11 * g22 - g12 * g12, 0.0)))
    face_areas = np.stack(face_areas, axis=1)
    surface_areas = np.sum(face_areas, axis=1)
    inradii = 3.0 * volumes / surface_areas
    shape_quality = 12.0 * np.power(3.0 * volumes, 2.0 / 3.0) / edge_squared_sum
    if (np.any(inradii <= 0) or np.any(shape_quality <= 0)
            or not np.all(np.isfinite(shape_quality))):
        raise ValueError("invalid tetrahedral shape-quality diagnostic")
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
        "shape_quality_min": float(np.min(shape_quality)),
        "shape_quality_p01": float(np.quantile(shape_quality, 0.01)),
        "shape_quality_median": float(np.median(shape_quality)),
        "shape_quality_max": float(np.max(shape_quality)),
        "inradius_min": float(np.min(inradii)),
        "max_diameter_to_inradius": float(np.max(diameters / inradii)),
    }, volumes


def refine_embedding(coordinates, tetrahedra):
    children, faces, _chain_by_tet = barycentric_subdivision(tetrahedra)
    old = np.asarray(coordinates, dtype=float)
    new_coordinates = np.empty((len(faces), old.shape[1]), dtype=float)
    for face_id, face in enumerate(faces):
        new_coordinates[face_id] = np.mean(old[np.asarray(face, dtype=np.int64)], axis=0)
    return new_coordinates, children, faces


def stable_child_order(tetrahedra, faces):
    """Match the exporter's stable ID: 24*parent_id + permutation rank."""
    face_id = {tuple(face): index for index, face in enumerate(faces)}
    ordered = np.empty((24 * len(tetrahedra), 4), dtype=np.uint32)
    for parent_id, tetrahedron in enumerate(tetrahedra):
        tet = tuple(int(v) for v in tetrahedron)
        for rank, permutation in enumerate(itertools.permutations(range(4))):
            prefix = []
            chain = []
            for position in permutation:
                prefix.append(tet[position])
                chain.append(face_id[tuple(sorted(prefix))])
            ordered[24 * parent_id + rank] = np.sort(chain)
    return ordered


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


def coarse_finite_volume_operator(edge_endpoints, conductances, fibre_labels,
                                   cell_volumes):
    """Mass-normalized Galerkin compression to fibre-constant functions."""
    import scipy.sparse as sp

    edges = np.asarray(edge_endpoints, dtype=np.int64)
    weights = np.asarray(conductances, dtype=float)
    labels = np.asarray(fibre_labels, dtype=np.int64)
    volumes = np.asarray(cell_volumes, dtype=float)
    if labels.ndim != 1 or volumes.shape != labels.shape:
        raise ValueError("fibre labels and cell volumes must be matching vectors")
    if edges.ndim != 2 or edges.shape[1] != 2 or weights.shape != (len(edges),):
        raise ValueError("dual edges and conductances have incompatible shapes")
    if np.any(labels < 0) or np.any(volumes <= 0) or np.any(weights <= 0):
        raise ValueError("labels, volumes, and conductances must be positive/valid")
    n_coarse = int(labels.max()) + 1
    mass = np.bincount(labels, weights=volumes, minlength=n_coarse)
    if np.any(mass <= 0):
        raise ValueError("every coarse fibre must have positive geometric mass")
    left = labels[edges[:, 0]]
    right = labels[edges[:, 1]]
    crossing = left != right
    a, b, w = left[crossing], right[crossing], weights[crossing]
    rows = np.concatenate((a, b, a, b))
    cols = np.concatenate((a, b, b, a))
    data = np.concatenate((w, w, -w, -w))
    stiffness = sp.coo_matrix((data, (rows, cols)),
                              shape=(n_coarse, n_coarse)).tocsr()
    dense = stiffness.toarray()
    inverse_sqrt_mass = 1.0 / np.sqrt(mass)
    operator = (inverse_sqrt_mass[:, None] * dense) * inverse_sqrt_mass[None, :]
    operator = 0.5 * (operator + operator.T)
    return operator, mass


def coarse_spectral_report(generation, operator, mass):
    eigenvalues = np.linalg.eigvalsh(operator)
    tolerance = max(1e-10, float(np.max(np.abs(eigenvalues), initial=0.0)) * 1e-12)
    positive = eigenvalues[eigenvalues > tolerance]
    heat_times = (1e-5, 1e-4, 1e-3, 1e-2)
    return {
        "generation": int(generation),
        "coarse_dimension": int(len(eigenvalues)),
        "coarse_mass_total": float(np.sum(mass)),
        "coarse_mass_min": float(np.min(mass)),
        "coarse_mass_max": float(np.max(mass)),
        "zero_modes_at_relative_tolerance": int(np.count_nonzero(np.abs(eigenvalues) <= tolerance)),
        "zero_mode_tolerance": tolerance,
        "smallest_positive_eigenvalue": float(positive[0]) if len(positive) else None,
        "largest_eigenvalue": float(eigenvalues[-1]),
        "trace": float(np.sum(eigenvalues)),
        "symmetry_error": float(np.max(np.abs(operator - operator.T), initial=0.0)),
        "minimum_eigenvalue": float(eigenvalues[0]),
        "first_12_eigenvalues": eigenvalues[:12].tolist(),
        "heat_trace": {
            str(t): float(np.sum(np.exp(-t * eigenvalues)))
            for t in heat_times
        },
        "scope": "finite-volume geometric operator compressed to causal-fibre-constant functions; not a continuum convergence result",
    }


def coarse_space_approximation_report(coordinates, tetrahedra, fibre_labels,
                                       cell_volumes):
    """Test whether causal fibres approximate simple PL coordinate functions."""
    points = np.asarray(coordinates, dtype=float)
    tets = np.asarray(tetrahedra, dtype=np.int64)
    labels = np.asarray(fibre_labels, dtype=np.int64)
    volumes = np.asarray(cell_volumes, dtype=float)
    if labels.shape != volumes.shape or len(labels) != len(tets):
        raise ValueError("coarse-space approximation inputs have incompatible shapes")
    n_fibres = int(labels.max()) + 1
    cell_centres = points[tets].mean(axis=1)
    fibre_mass = np.bincount(labels, weights=volumes, minlength=n_fibres)
    fibre_centres = np.column_stack([
        np.bincount(labels, weights=volumes * cell_centres[:, axis],
                    minlength=n_fibres) / fibre_mass
        for axis in range(points.shape[1])
    ])
    coordinate_errors = []
    for axis in range(points.shape[1]):
        values = cell_centres[:, axis]
        projected = fibre_centres[labels, axis]
        denominator = float(np.sum(volumes * values * values))
        numerator = float(np.sum(volumes * (values - projected) ** 2))
        coordinate_errors.append(
            float(np.sqrt(numerator / denominator)) if denominator > 0 else 0.0
        )
    displacement = cell_centres - fibre_centres[labels]
    squared_radius = np.sum(displacement * displacement, axis=1)
    rms_radius_by_fibre = np.sqrt(
        np.bincount(labels, weights=volumes * squared_radius,
                    minlength=n_fibres) / fibre_mass
    )
    max_radius_by_fibre = np.zeros(n_fibres, dtype=float)
    np.maximum.at(max_radius_by_fibre, labels, np.sqrt(squared_radius))
    return {
        "test_functions": "four ambient PL coordinate functions x_0,...,x_3 sampled at tetrahedron barycentres",
        "relative_l2_projection_errors": coordinate_errors,
        "rms_fibre_radius_min": float(rms_radius_by_fibre.min()),
        "rms_fibre_radius_median": float(np.median(rms_radius_by_fibre)),
        "rms_fibre_radius_max": float(rms_radius_by_fibre.max()),
        "max_fibre_radius": float(max_radius_by_fibre.max()),
        "coordinate_modes_resolved": bool(max(coordinate_errors) < 0.5),
        "scope": "finite approximation-space diagnostic; coordinate tests do not by themselves prove density or disprove all possible limits",
    }


def build_coarse_spectral_sequence():
    """Build G0--G2 geometric dual operators and their causal-fibre compressions."""
    coordinates, tetrahedra = initial_reference_embedding()
    parent_map = {tuple(t): "e0" for t in tetrahedra}
    orientations = {
        tuple(t): (-1 if sum(v % 2 for v in t) % 2 else 1)
        for t in tetrahedra
    }
    causal_edges = {"e0": (0, 1)}
    next_vertex = 2
    chain_by_tet = None
    previous_faces = None
    reports = []

    for generation in range(3):
        _metrics, volumes = mesh_metrics(coordinates, tetrahedra)
        dual = dual_graph_payload(generation, tetrahedra, chain_by_tet, previous_faces)
        endpoints = np.asarray(
            [(edge["source"], edge["target"]) for edge in dual["edges"]],
            dtype=np.uint32,
        )
        weights, _areas, _distances = geometric_dual_weights(
            coordinates, tetrahedra, endpoints, volumes, chunk_rows=50000
        )
        fibre_ids = sorted(set(parent_map.values()))
        fibre_index = {key: index for index, key in enumerate(fibre_ids)}
        labels = np.asarray(
            [fibre_index[parent_map[tuple(tet)]] for tet in tetrahedra],
            dtype=np.uint16,
        )
        if len(fibre_ids) != (1, 8, 64)[generation]:
            raise ValueError(f"unexpected causal-fibre count at G{generation}")
        operator, mass = coarse_finite_volume_operator(
            endpoints, weights, labels, volumes
        )
        report = coarse_spectral_report(generation, operator, mass)
        report.update({
            "fine_vertices": int(len(tetrahedra)),
            "fine_edges": int(len(endpoints)),
            "fine_graph_components": int(dual["component_count"]),
            "coarse_space_approximation": coarse_space_approximation_report(
                coordinates, tetrahedra, labels, volumes
            ),
        })
        reports.append(report)

        if generation < 2:
            (children, child_map, child_orientations, child_causal,
             next_vertex, previous_faces, chain_by_tet) = _advance_one_generation(
                tetrahedra, parent_map, orientations, causal_edges,
                next_vertex, generation, chain_by_tet, previous_faces,
            )
            coordinates, generated_children, _faces = refine_embedding(
                coordinates, tetrahedra
            )
            if not np.array_equal(np.asarray(children), np.asarray(generated_children)):
                raise ValueError("refinement helper disagrees with causal exporter ordering")
            tetrahedra = generated_children
            parent_map = child_map
            orientations = child_orientations
            causal_edges = child_causal

    return reports


def _sha256(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def assess_shape_regularity(levels):
    quality = [float(row["shape_quality_min"]) for row in levels]
    aspect = [float(row["max_diameter_to_inradius"]) for row in levels]
    quality_degrades = any(b < a - 1e-12 for a, b in zip(quality, quality[1:]))
    aspect_grows = any(b > a + 1e-12 for a, b in zip(aspect, aspect[1:]))
    return {
        "minimum_shape_quality_by_generation": quality,
        "maximum_diameter_to_inradius_by_generation": aspect,
        "shape_quality_degradation_detected": bool(quality_degrades),
        "aspect_ratio_growth_detected": bool(aspect_grows),
        "uniform_shape_regularity_established": False,
        "warning": (
            "Barycentric refinement shows deteriorating tetrahedral shape quality; "
            "standard shape-regular finite-element estimates cannot be assumed "
            "without a separate proof or a better refinement family."
        ) if quality_degrades and aspect_grows else (
            "Finite generations do not establish uniform shape regularity."
        ),
    }


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
            coordinates, children, faces = refine_embedding(coordinates, tetrahedra)
            tetrahedra = (
                stable_child_order(tetrahedra, faces)
                if generation + 1 == max_generation
                else children
            )
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
    fibre_labels = np.load(chunk_dir / "g3_fibre_labels.npy", mmap_mode="r")
    g3_coarse_operator, g3_coarse_mass = coarse_finite_volume_operator(
        endpoints, weights, fibre_labels, volumes
    )
    g3_coarse_report = coarse_spectral_report(
        3, g3_coarse_operator, g3_coarse_mass
    )
    coarse_spectral_sequence = build_coarse_spectral_sequence()
    coarse_spectral_sequence.append(g3_coarse_report)

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
        "g3_geometric_coarse_operator.npy": g3_coarse_operator,
        "g3_geometric_coarse_eigenvalues.npy": np.linalg.eigvalsh(g3_coarse_operator),
    }
    for filename, array in arrays.items():
        np.save(output_dir / filename, array, allow_pickle=False)

    manifest = {
        "schema": "q2-pl-geometric-mesh-audit-v1",
        "status": "PASS_G3_PL_GEOMETRIC_MESH_AUDIT",
        "passed": True,
        "reference_embedding": "cross-polytope boundary vertices +/-e_i in R^4; all refined vertices are barycentres of faces",
        "levels": levels,
        "shape_regularity_assessment": assess_shape_regularity(levels),
        "coarse_geometric_spectral_sequence": coarse_spectral_sequence,
        "g3": {
            "coarse_spectral": g3_coarse_report,
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
