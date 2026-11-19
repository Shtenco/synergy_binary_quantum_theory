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
    if not np.isfinite(lam_min) or lam_min <= 0 or not np.isfinite(lam_max):
        raise ValueError(f"Q-space operator not numerically positive: {lam_min}, {lam_max}")

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
        lam_min > 0 and small_residual <= max(1e-6, spectral_tol * 10)
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
            "method": "scipy eigsh on matrix-free QLQ+P; Ritz estimates, not certified bounds",
            "tolerance": spectral_tol,
        },
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
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\\n",
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
