#!/usr/bin/env python3
"""Pull back the projective boundary-state tangent to internal holonomy probes.

An edge holonomy probe is represented in the existing contraction by a twist
on the first endpoint half-link. This gate constructs that injection explicitly
and verifies D_edge = D_half P and F_edge = P^T F_half P. It compares ranks and
nullspaces without identifying the graph-wide vertex-gauge quotient with the
projective boundary-state stabilizer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_relational_graph_qfi_gate import (
    graph_fixtures,
    graph_qfi,
    full_half_link_qfi,
    random_intertwiner,
)
from q2_relational_gauge_quotient_gate import projective_tangent_jacobian
from q2_relational_holonomy_quotient_gate import (
    endpoint_to_holonomy_jacobian,
    fixture_holonomies,
    vertex_gauge_embedding,
    holonomy_gauge_jacobian,
)


def edge_probe_injection(
    node_count: int, edges: list[tuple[int, int, int, int]]
) -> np.ndarray:
    """P: R^(3E) -> R^(12V), twist each edge's first endpoint only."""
    result = np.zeros((12 * node_count, 3 * len(edges)), dtype=float)
    occupied = set()
    for edge_index, (u, leg_u, v, leg_v) in enumerate(edges):
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("invalid edge endpoints")
        if not (0 <= leg_u < 4 and 0 <= leg_v < 4):
            raise ValueError("invalid half-link index")
        if (u, leg_u) in occupied or (v, leg_v) in occupied:
            raise ValueError("half-link reused by multiple edges")
        occupied.add((u, leg_u))
        occupied.add((v, leg_v))
        for axis in range(3):
            result[12 * u + 3 * leg_u + axis, 3 * edge_index + axis] = 1.0
    return result


def realify(matrix: np.ndarray) -> np.ndarray:
    return np.concatenate([matrix.real, matrix.imag], axis=0)


def rank_and_spectrum(matrix: np.ndarray, relative_tolerance: float = 1e-9):
    singular = np.linalg.svd(matrix, compute_uv=False)
    scale = float(singular[0]) if singular.size else 0.0
    tolerance = relative_tolerance * max(1.0, scale)
    return int(np.sum(singular > tolerance)), singular, tolerance


def compare_one(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
) -> dict[str, object]:
    node_count = len(node_states)
    _p, full_qfi, _sym, _boundary = full_half_link_qfi(node_states, edges)
    half_tangent, _probability = projective_tangent_jacobian(node_states, edges)
    injection = edge_probe_injection(node_count, edges)
    pulled_tangent = half_tangent @ injection
    edge_rank, edge_singular, edge_tol = rank_and_spectrum(realify(pulled_tangent))
    edge_qfi_from_pullback = injection.T @ full_qfi @ injection
    edge_qfi_from_pullback = 0.5 * (edge_qfi_from_pullback + edge_qfi_from_pullback.T)
    _probability2, direct_edge_qfi, direct_symmetry, _boundary2 = graph_qfi(node_states, edges)
    qfi_eigenvalues = np.linalg.eigvalsh(edge_qfi_from_pullback)
    qfi_scale = max(1.0, float(np.max(np.abs(qfi_eigenvalues))))
    qfi_rank = int(np.sum(qfi_eigenvalues > 1e-9 * qfi_scale))
    reconstructed = 4.0 * realify(pulled_tangent).T @ realify(pulled_tangent)

    # The holonomy-space gauge Jacobian is deliberately reported separately.
    # It is not silently subtracted from the boundary-state tangent: the
    # boundary frame transforms too, and that action is not yet implemented.
    fixture_name = None
    for name, (n, fixture_edges) in graph_fixtures().items():
        if n == node_count and fixture_edges == edges:
            fixture_name = name
            break
    holonomy_gauge = {}
    if fixture_name is not None:
        holonomies = fixture_holonomies(fixture_name)
        endpoint_map = endpoint_to_holonomy_jacobian(holonomies)
        vertex_map = vertex_gauge_embedding(node_count, [(e[0], e[2]) for e in edges])
        gauge_jac = endpoint_map @ vertex_map
        holonomy_gauge = {
            "vertex_gauge_orbit_rank": int(np.linalg.matrix_rank(gauge_jac, tol=1e-9)),
            "holonomy_moduli_dimension": int(3 * len(edges) - np.linalg.matrix_rank(gauge_jac, tol=1e-9)),
            "not_identified_with_boundary_tangent_kernel": True,
        }

    return {
        "node_count": node_count,
        "edge_count": len(edges),
        "half_link_tangent_shape": list(half_tangent.shape),
        "edge_injection_shape": list(injection.shape),
        "edge_projective_tangent_shape": list(pulled_tangent.shape),
        "edge_projective_tangent_rank": edge_rank,
        "edge_projective_tangent_nullity": 3 * len(edges) - edge_rank,
        "edge_projective_tangent_singular_values": edge_singular.tolist(),
        "edge_projective_tangent_tolerance": edge_tol,
        "edge_qfi_pullback_rank": qfi_rank,
        "edge_qfi_pullback_minimum_eigenvalue": float(qfi_eigenvalues[0]),
        "edge_qfi_reconstruction_error": float(np.linalg.norm(reconstructed - edge_qfi_from_pullback)),
        "direct_graph_qfi_pullback_error": float(np.linalg.norm(direct_edge_qfi - edge_qfi_from_pullback)),
        "direct_graph_qfi_symmetry_error": direct_symmetry,
        "holonomy_gauge_diagnostics": holonomy_gauge,
        "checks": {
            "qfi_equals_four_projective_tangent_gram": np.allclose(
                reconstructed, edge_qfi_from_pullback, atol=3e-10, rtol=3e-10
            ),
            "direct_edge_qfi_matches_half_link_pullback": np.allclose(
                direct_edge_qfi, edge_qfi_from_pullback, atol=3e-10, rtol=3e-10
            ),
            "tangent_rank_matches_qfi_rank": edge_rank == qfi_rank,
            "edge_qfi_psd_up_to_roundoff": qfi_eigenvalues[0] >= -1e-9 * qfi_scale,
        },
    }


def run(seed: int = 20261118, trials: int = 4) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    results = {}
    for name, (node_count, edges) in graph_fixtures().items():
        rows = [
            compare_one([random_intertwiner(rng) for _ in range(node_count)], edges)
            for _ in range(trials)
        ]
        results[name] = {
            "nodes": node_count,
            "edges": len(edges),
            "trials": trials,
            "edge_projective_tangent_rank_min": min(
                row["edge_projective_tangent_rank"] for row in rows
            ),
            "edge_projective_tangent_nullity_max": max(
                row["edge_projective_tangent_nullity"] for row in rows
            ),
            "max_qfi_reconstruction_error": max(
                row["edge_qfi_reconstruction_error"] for row in rows
            ),
            "max_direct_graph_qfi_pullback_error": max(
                row["direct_graph_qfi_pullback_error"] for row in rows
            ),
            "holonomy_gauge_diagnostics": rows[0]["holonomy_gauge_diagnostics"],
            "checks": {
                key: all(row["checks"][key] for row in rows)
                for key in rows[0]["checks"]
            },
            "samples": rows,
        }
        results[name]["passed"] = all(results[name]["checks"].values())
    return {
        "status": "finite holonomy-to-projective-boundary tangent pullback",
        "passed": all(row["passed"] for row in results.values()),
        "seed": seed,
        "trials_per_graph": trials,
        "map": "D_edge = D_half_link P; F_edge = P^T F_half_link P",
        "important_scope": (
            "The internal-edge twist is inserted at the first endpoint, matching "
            "the current graph_qfi convention. This proves the tangent/QFI pullback "
            "for that chosen probe section. It does not yet implement the simultaneous "
            "vertex action on internal holonomies AND boundary frames; therefore the "
            "holonomy gauge quotient and projective boundary tangent are reported "
            "separately and are not identified."
        ),
        "graphs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.seed, args.trials)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
