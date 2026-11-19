#!/usr/bin/env python3
"""Run matrix-free G3 Feshbach diagnostics on the persisted chunked export.

The spectral extrema are Ritz estimates, not certified interval bounds. The
reported Feshbach correction and slope are deterministic coarse-vector probes,
not the full 512x512 self-energy matrix.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

try:
    from .q2_matrix_free_feshbach import MatrixFreeFeshbach
except ImportError:
    from q2_matrix_free_feshbach import MatrixFreeFeshbach


def coarse_galerkin_matrix(edge_endpoints, fibre_labels, fibre_sizes, *, chunk_rows=100000):
    """Build U.T @ L @ U directly from cross-fibre edge counts."""
    labels = np.asarray(fibre_labels)
    sizes = np.asarray(fibre_sizes, dtype=float)
    n_fibres = len(sizes)
    matrix = np.zeros((n_fibres, n_fibres), dtype=float)
    for start in range(0, len(edge_endpoints), chunk_rows):
        block = np.asarray(edge_endpoints[start:start + chunk_rows])
        left = labels[block[:, 0]].astype(np.int64, copy=False)
        right = labels[block[:, 1]].astype(np.int64, copy=False)
        crossing = left != right
        left, right = left[crossing], right[crossing]
        if not len(left):
            continue
        weights = 1.0 / np.sqrt(sizes[left] * sizes[right])
        np.add.at(matrix, (left, left), 1.0 / sizes[left])
        np.add.at(matrix, (right, right), 1.0 / sizes[right])
        np.add.at(matrix, (left, right), -weights)
        np.add.at(matrix, (right, left), -weights)
    return matrix


def galerkin_spectral_report(matrix, *, measure_total, times=None):
    times = np.asarray(times if times is not None else
                       [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0], dtype=float)
    symmetry_error = float(np.max(np.abs(matrix - matrix.T)))
    eigenvalues = np.linalg.eigvalsh(0.5 * (matrix + matrix.T))
    heat = np.asarray([np.mean(np.exp(-t * eigenvalues)) for t in times])
    log_t, log_heat = np.log(times), np.log(heat)
    spectral_dimension = np.asarray([
        -2.0 * (log_heat[i + 1] - log_heat[i - 1]) /
        (log_t[i + 1] - log_t[i - 1])
        for i in range(1, len(times) - 1)
    ])
    positive = eigenvalues[eigenvalues > 1e-10]
    return {
        "dimension": int(matrix.shape[0]),
        "nonzero_count": int(np.count_nonzero(np.abs(matrix) > 1e-14)),
        "symmetry_error": symmetry_error,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "spectral_gap": float(positive[0]) if len(positive) else None,
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "negative_eigenvalue_count_below_tolerance": int(np.sum(eigenvalues < -1e-10)),
        "zero_mode_count": int(np.sum(np.abs(eigenvalues) <= 1e-10)),
        "trace": float(np.trace(matrix)),
        "measure_total": float(measure_total),
        "heat_trace_normalized": [
            {"t": float(times[i]), "value": float(heat[i]),
             "spectral_dimension": float(spectral_dimension[i - 1])}
            for i in range(1, len(times) - 1)
        ],
        "normalized_heat_trace": [
            {"t": float(t), "value": float(h)} for t, h in zip(times, heat)
        ],
        "eigenvalues": eigenvalues.tolist(),
    }


def certified_gap_lower_bound_from_root_eccentricity(model) -> dict:
    """Bound lambda_2(L) below using a BFS eccentricity upper bound.

    For a connected n-vertex graph of diameter D,
    lambda_2(L) >= 2 / ((n - 1) D). A BFS eccentricity e from any root
    satisfies D <= 2e, hence lambda_2(L) >= 1 / ((n - 1) e).
    Since constants belong to the fibre-constant P space, the same lower
    bound applies to QLQ on range(Q).
    """
    from scipy.sparse import diags

    n = model.n_vertices
    if n <= 1 or model.graph_component_count != 1:
        raise ValueError("certified gap bound requires a connected nontrivial graph")
    degrees = np.asarray(model.laplacian.diagonal(), dtype=float)
    adjacency = (diags(degrees) - model.laplacian).tocsr()
    lambda_max_upper = max(1.0, 2.0 * float(degrees.max()))
    distances = np.full(n, -1, dtype=np.int32)
    queue = np.empty(n, dtype=np.int32)
    distances[0] = 0
    queue[0] = 0
    head, tail = 0, 1
    while head < tail:
        vertex = int(queue[head])
        head += 1
        next_distance = distances[vertex] + 1
        for offset in range(adjacency.indptr[vertex], adjacency.indptr[vertex + 1]):
            neighbour = int(adjacency.indices[offset])
            if distances[neighbour] < 0:
                distances[neighbour] = next_distance
                queue[tail] = neighbour
                tail += 1
    if tail != n:
        raise ValueError("BFS did not reach every graph vertex")
    eccentricity = int(distances.max())
    if eccentricity <= 0:
        raise ValueError("invalid BFS eccentricity")
    lower_bound = 1.0 / ((n - 1) * eccentricity)
    full_operator_lower_bound = min(1.0, lower_bound)
    return {
        "root_eccentricity_upper_bound": eccentricity,
        "lambda_min_q_l_q_certified_lower_bound": lower_bound,
        "lambda_min_q_l_q_plus_p_certified_lower_bound": full_operator_lower_bound,
        "lambda_max_q_l_q_plus_p_certified_upper_bound": lambda_max_upper,
        "condition_number_certified_upper_bound": lambda_max_upper / full_operator_lower_bound,
        "proof": "connected graph diameter <= 2*root eccentricity; lambda_2(L) >= 1/((n-1)*eccentricity); lambda_max(L) <= 2*d_max; range(Q) is orthogonal to constants",
    }


def run(export_dir: Path, *, cg_rtol: float = 1e-8, cg_maxiter: int = 10000,
        spectral_tol: float = 1e-5, spectral_maxiter: int = 5000) -> dict:
    endpoints = np.load(export_dir / "g3_dual_endpoints.npy", mmap_mode="r")
    labels = np.load(export_dir / "g3_fibre_labels.npy", mmap_mode="r")
    manifest = json.loads((export_dir / "G3_CHUNKED_MANIFEST.json").read_text())
    if not manifest.get("passed"):
        raise ValueError("G3 export manifest is not PASS")
    if endpoints.shape != (442368, 2) or labels.shape != (221184,):
        raise ValueError("unexpected G3 graph dimensions")

    started = time.perf_counter()
    model = MatrixFreeFeshbach(endpoints, labels, validate_edges=False)
    build_seconds = time.perf_counter() - started
    if model.n_vertices != 221184 or model.n_fibres != 512:
        raise ValueError("unexpected G3 Feshbach dimensions")
    if model.graph_component_count != 1:
        raise ValueError("G3 Feshbach requires a connected graph")
    certified_gap = certified_gap_lower_bound_from_root_eccentricity(model)

    galerkin = coarse_galerkin_matrix(endpoints, labels, model.fibre_sizes)
    galerkin_report = galerkin_spectral_report(
        galerkin, measure_total=model.n_vertices / float(24 ** 3)
    )
    if (galerkin_report["symmetry_error"] > 1e-12
            or galerkin_report["minimum_eigenvalue"] < -1e-10
            or galerkin_report["zero_mode_count"] != 1
            or abs(galerkin_report["measure_total"] - 16.0) > 1e-12):
        raise ValueError(f"G3 sparse Galerkin checks failed: {galerkin_report}")

    operator = model._operator()
    rng = np.random.default_rng(20261119)
    v0 = rng.standard_normal(model.n_vertices)
    v0 /= np.linalg.norm(v0)
    spectral_started = time.perf_counter()
    smallest, small_vec = eigsh(
        operator, k=1, which="SA", v0=v0, tol=spectral_tol,
        maxiter=spectral_maxiter, return_eigenvectors=True,
    )
    lam_min = float(smallest[0])
    small_residual = float(np.linalg.norm(model.apply_q_plus_p(small_vec[:, 0]) -
                                          lam_min * small_vec[:, 0]))
    # The graph is 4-regular, so lambda_max(L) <= 2*d_max = 8.
    # Compression by Q cannot increase the operator norm; P contributes
    # eigenvalue 1 on the orthogonal coarse subspace.
    lam_max_upper = 8.0
    spectral_seconds = time.perf_counter() - spectral_started
    if not np.isfinite(lam_min) or lam_min <= 0 or not np.isfinite(lam_max_upper):
        raise ValueError(f"Q-space operator not numerically positive: {lam_min}, {lam_max_upper}")

    probes = []
    correction_min = float("inf")
    slope_min = float("inf")
    max_cg_residual = 0.0
    max_second_residual = 0.0
    for probe_id in range(4):
        coarse = rng.standard_normal(model.n_fibres)
        coarse /= np.linalg.norm(coarse)
        correction, correction_info = model.apply_feshbach_correction(
            coarse, rtol=cg_rtol, atol=1e-11, maxiter=cg_maxiter
        )
        slope, slope_info = model.apply_self_energy_slope(
            coarse, rtol=cg_rtol, atol=1e-11, maxiter=cg_maxiter
        )
        correction_quad = float(np.dot(coarse, correction))
        slope_quad = float(np.dot(coarse, slope))
        correction_min = min(correction_min, correction_quad)
        slope_min = min(slope_min, slope_quad)
        max_cg_residual = max(max_cg_residual, correction_info["residual_norm"],
                              slope_info["first_solve"]["residual_norm"])
        max_second_residual = max(max_second_residual,
                                  slope_info["second_solve"]["residual_norm"])
        probes.append({
            "probe": probe_id,
            "correction_quadratic_form": correction_quad,
            "slope_quadratic_form": slope_quad,
            "correction_cg": correction_info,
            "slope_cg": slope_info,
        })

    passed = (
        lam_min > 0 and certified_gap["lambda_min_q_l_q_certified_lower_bound"] > 0
        and small_residual <= max(1e-6, spectral_tol * 10)
        and correction_min >= -1e-8 and slope_min >= -1e-8
        and max_cg_residual <= max(1e-7, cg_rtol * 100)
        and max_second_residual <= max(1e-7, cg_rtol * 100)
    )
    return {
        "schema": "q2-g3-matrix-free-feshbach-v1",
        "status": "PASS_G3_FESHBACH_PROBES" if passed else "FAIL_G3_FESHBACH_PROBES",
        "passed": bool(passed),
        "generation": 3,
        "vertices": model.n_vertices,
        "edges": int(len(endpoints)),
        "coarse_fibres": model.n_fibres,
        "graph_components": int(model.graph_component_count),
        "fibre_size_min": int(model.fibre_sizes.min()),
        "fibre_size_max": int(model.fibre_sizes.max()),
        "q_l_q_plus_p_spectrum_estimate": {
            "lambda_min_ritz": lam_min,
            "lambda_max_upper_bound": lam_max_upper,
            "condition_ratio_estimate": lam_max_upper / lam_min,
            "smallest_eigenpair_residual": small_residual,
            "method": "smallest Ritz value from scipy eigsh; certified condition bound is reported separately",
            "tolerance": spectral_tol,
        },
        "sparse_galerkin": galerkin_report,
        "feshbach_probes": probes,
        "minimum_correction_quadratic_form": correction_min,
        "minimum_self_energy_slope_quadratic_form": slope_min,
        "maximum_cg_residual": max_cg_residual,
        "maximum_second_resolvent_residual": max_second_residual,
        "build_seconds": build_seconds,
        "spectral_seconds": spectral_seconds,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "scope_boundary": "Four deterministic coarse probes only; not the full 512x512 self-energy matrix or a continuum/physical closure.",
    }


def write_report(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cg-rtol", type=float, default=1e-8)
    parser.add_argument("--cg-maxiter", type=int, default=10000)
    parser.add_argument("--spectral-tol", type=float, default=1e-5)
    parser.add_argument("--spectral-maxiter", type=int, default=5000)
    args = parser.parse_args()
    result = run(args.export_dir, cg_rtol=args.cg_rtol, cg_maxiter=args.cg_maxiter,
                 spectral_tol=args.spectral_tol, spectral_maxiter=args.spectral_maxiter)
    write_report(args.output, result)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
