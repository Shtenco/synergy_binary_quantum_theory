#!/usr/bin/env python3
"""Compute the complete finite G3 Feshbach correction and self-energy matrices.

For every coarse basis vector e_j this evaluates
B.T (Q L Q|_Q)^-1 B e_j and B.T (Q L Q|_Q)^-2 B e_j using projected CG.
The resulting 512x512 matrices are finite geometric graph operators; they are
not physical QFI conductance matrices.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np

try:
    from .q2_matrix_free_feshbach import MatrixFreeFeshbach
except ImportError:
    from q2_matrix_free_feshbach import MatrixFreeFeshbach


def full_feshbach_matrices(model: MatrixFreeFeshbach, *, rtol=1e-8,
                           atol=1e-11, maxiter=10000) -> tuple[np.ndarray, np.ndarray, dict]:
    """Assemble full coarse correction and slope matrices column by column."""
    m = model.n_fibres
    correction = np.empty((m, m), dtype=float)
    slope = np.empty((m, m), dtype=float)
    max_first_residual = 0.0
    max_second_residual = 0.0
    max_first_iterations = 0
    max_second_iterations = 0
    total_iterations = 0
    started = time.perf_counter()

    for column in range(m):
        coarse = np.zeros(m, dtype=float)
        coarse[column] = 1.0
        rhs = model.apply_b(coarse)
        first, first_info = model.solve_q(
            rhs, rtol=rtol, atol=atol, maxiter=maxiter
        )
        second, second_info = model.solve_q(
            first, rtol=rtol, atol=atol, maxiter=maxiter
        )
        correction[:, column] = model.apply_bt(first)
        slope[:, column] = model.apply_bt(second)
        max_first_residual = max(max_first_residual, first_info["residual_norm"])
        max_second_residual = max(max_second_residual, second_info["residual_norm"])
        max_first_iterations = max(max_first_iterations, first_info["iterations"])
        max_second_iterations = max(max_second_iterations, second_info["iterations"])
        total_iterations += first_info["iterations"] + second_info["iterations"]

    elapsed = time.perf_counter() - started
    return correction, slope, {
        "columns_solved": m,
        "max_first_solve_residual": max_first_residual,
        "max_second_solve_residual": max_second_residual,
        "max_first_solve_iterations": max_first_iterations,
        "max_second_solve_iterations": max_second_iterations,
        "total_cg_iterations": total_iterations,
        "solve_seconds": elapsed,
    }


def matrix_diagnostics(matrix: np.ndarray, *, tolerance=1e-7) -> dict:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Feshbach matrix must be square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Feshbach matrix contains non-finite entries")
    symmetry_error = float(np.max(np.abs(matrix - matrix.T)))
    eigenvalues = np.linalg.eigvalsh(0.5 * (matrix + matrix.T))
    return {
        "dimension": int(matrix.shape[0]),
        "symmetry_error": symmetry_error,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "negative_eigenvalue_count_below_tolerance": int(
            np.sum(eigenvalues < -tolerance)
        ),
        "trace": float(np.trace(matrix)),
        "frobenius_norm": float(np.linalg.norm(matrix, ord="fro")),
        "minimum_diagonal": float(np.min(np.diag(matrix))),
        "maximum_diagonal": float(np.max(np.diag(matrix))),
        "eigenvalues": eigenvalues.tolist(),
    }


def run(export_dir: Path, output_dir: Path, *, rtol=1e-8, atol=1e-11,
        maxiter=10000) -> dict:
    endpoints = np.load(export_dir / "g3_dual_endpoints.npy", mmap_mode="r")
    labels = np.load(export_dir / "g3_fibre_labels.npy", mmap_mode="r")
    manifest = json.loads(
        (export_dir / "G3_CHUNKED_MANIFEST.json").read_text(encoding="utf-8")
    )
    if not manifest.get("passed"):
        raise ValueError("G3 chunked export manifest is not PASS")
    if endpoints.shape != (442368, 2) or labels.shape != (221184,):
        raise ValueError("unexpected G3 graph dimensions")

    model = MatrixFreeFeshbach(endpoints, labels, validate_edges=False)
    if model.graph_component_count != 1 or model.n_fibres != 512:
        raise ValueError("full G3 Feshbach requires connected graph and 512 fibres")

    correction, slope, solve_info = full_feshbach_matrices(
        model, rtol=rtol, atol=atol, maxiter=maxiter
    )
    correction_report = matrix_diagnostics(correction)
    slope_report = matrix_diagnostics(slope)
    symmetry_tolerance = max(1e-6, rtol * 1000)
    passed = (
        correction_report["symmetry_error"] <= symmetry_tolerance
        and slope_report["symmetry_error"] <= symmetry_tolerance
        and correction_report["minimum_eigenvalue"] >= -1e-7
        and slope_report["minimum_eigenvalue"] >= -1e-7
        and solve_info["max_first_solve_residual"] <= max(1e-7, rtol * 100)
        and solve_info["max_second_solve_residual"] <= max(1e-7, rtol * 100)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    np.save(output_dir / "G3_FESHBACH_CORRECTION_MATRIX.npy", correction,
            allow_pickle=False)
    np.save(output_dir / "G3_SELF_ENERGY_SLOPE_MATRIX.npy", slope,
            allow_pickle=False)
    report = {
        "schema": "q2-g3-full-feshbach-matrices-v1",
        "status": "PASS_G3_FULL_FESHBACH_MATRICES" if passed else "FAIL_G3_FULL_FESHBACH_MATRICES",
        "passed": bool(passed),
        "generation": 3,
        "vertices": model.n_vertices,
        "edges": int(len(endpoints)),
        "coarse_fibres": model.n_fibres,
        "graph_components": int(model.graph_component_count),
        "fibre_size_min": int(model.fibre_sizes.min()),
        "fibre_size_max": int(model.fibre_sizes.max()),
        "correction_matrix": correction_report,
        "self_energy_slope_matrix": slope_report,
        "solver": {
            **solve_info,
            "rtol": rtol,
            "atol": atol,
            "maxiter": maxiter,
        },
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "matrix_files": {
            "correction": "G3_FESHBACH_CORRECTION_MATRIX.npy",
            "self_energy_slope": "G3_SELF_ENERGY_SLOPE_MATRIX.npy",
        },
        "scope_boundary": (
            "Complete finite 512x512 geometric Feshbach correction and slope "
            "matrices. These are induced by the unweighted PL dual graph, not "
            "physical QFI conductances and not a continuum-limit proof."
        ),
    }
    (output_dir / "G3_FULL_FESHBACH_MANIFEST.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-11)
    parser.add_argument("--maxiter", type=int, default=10000)
    args = parser.parse_args()
    report = run(args.export_dir, args.output_dir, rtol=args.rtol,
                 atol=args.atol, maxiter=args.maxiter)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
