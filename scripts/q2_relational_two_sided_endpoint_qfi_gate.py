#!/usr/bin/env python3
"""Two-sided endpoint QFI and exact singlet-gluing redundancy.

For every internal edge, the contraction uses the invariant spin-1/2 singlet.
In the tensor-index convention of contract_graph, equal infinitesimal SU(2)
rotations on the two contracted endpoint indices leave the singlet invariant.
Consequently D_v(e) = -D_u(e) for the normalized boundary ray, and the
two-sided QFI has the block form [[F,-F],[-F,F]], including all cross-edge
correlations. This is a coordinate redundancy of the two-end probe, not a
claim of an independent edge-frame gauge field.
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
    random_intertwiner,
)
from q2_relational_gauge_quotient_gate import projective_tangent_jacobian


def two_sided_endpoint_injection(
    node_count: int, edges: list[tuple[int, int, int, int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return P_u, P_v and interleaved P_uv, each endpoint has 3 axes per edge."""
    if node_count < 1:
        raise ValueError("node_count must be positive")
    p_u = np.zeros((12 * node_count, 3 * len(edges)), dtype=float)
    p_v = np.zeros_like(p_u)
    p_uv = np.zeros((12 * node_count, 6 * len(edges)), dtype=float)
    occupied: set[tuple[int, int]] = set()
    for edge_index, edge in enumerate(edges):
        if len(edge) != 4:
            raise ValueError("edge must be (u, leg_u, v, leg_v)")
        u, leg_u, v, leg_v = map(int, edge)
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("invalid edge endpoints")
        if not (0 <= leg_u < 4 and 0 <= leg_v < 4):
            raise ValueError("invalid half-link index")
        if (u, leg_u) in occupied or (v, leg_v) in occupied:
            raise ValueError("a half-link cannot belong to more than one edge")
        occupied.add((u, leg_u))
        occupied.add((v, leg_v))
        for axis in range(3):
            row_u = 12 * u + 3 * leg_u + axis
            row_v = 12 * v + 3 * leg_v + axis
            p_u[row_u, 3 * edge_index + axis] = 1.0
            p_v[row_v, 3 * edge_index + axis] = 1.0
            p_uv[row_u, 6 * edge_index + axis] = 1.0
            p_uv[row_v, 6 * edge_index + 3 + axis] = 1.0
    return p_u, p_v, p_uv


def lifted_two_sided_qfi(edge_qfi: np.ndarray, edge_count: int) -> np.ndarray:
    """Build interleaved [[F,-F],[-F,F]] edge/endpoint blocks."""
    result = np.zeros((6 * edge_count, 6 * edge_count), dtype=float)
    for e in range(edge_count):
        for f in range(edge_count):
            block = edge_qfi[3*e:3*e+3, 3*f:3*f+3]
            result[6*e:6*e+3, 6*f:6*f+3] = block
            result[6*e:6*e+3, 6*f+3:6*f+6] = -block
            result[6*e+3:6*e+6, 6*f:6*f+3] = -block
            result[6*e+3:6*e+6, 6*f+3:6*f+6] = block
    return result


def analyze_sample(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
) -> dict[str, object]:
    node_count, edge_count = len(node_states), len(edges)
    tangent, probability = projective_tangent_jacobian(node_states, edges)
    p_u, p_v, p_uv = two_sided_endpoint_injection(node_count, edges)
    d_u, d_v = tangent @ p_u, tangent @ p_v
    d_uv = tangent @ p_uv
    _prob, edge_qfi, edge_symmetry, _boundary = graph_qfi(node_states, edges)
    two_qfi = 4.0 * np.real(d_uv.conj().T @ d_uv)
    two_qfi = 0.5 * (two_qfi + two_qfi.T)
    expected = lifted_two_sided_qfi(edge_qfi, edge_count)
    eigenvalues = np.linalg.eigvalsh(two_qfi)
    edge_eigenvalues = np.linalg.eigvalsh(edge_qfi)
    scale = max(1.0, float(np.max(np.abs(eigenvalues))))
    tolerance = 1e-9 * scale
    two_rank = int(np.sum(eigenvalues > tolerance))
    edge_scale = max(1.0, float(np.max(np.abs(edge_eigenvalues))))
    edge_rank = int(np.sum(edge_eigenvalues > 1e-9 * edge_scale))
    kernel_vectors = []
    for edge in range(edge_count):
        for axis in range(3):
            direction = np.zeros(6 * edge_count)
            direction[6 * edge + axis] = 1.0
            direction[6 * edge + 3 + axis] = 1.0
            kernel_vectors.append(direction)
    kernel_residuals = [float(np.linalg.norm(two_qfi @ v)) for v in kernel_vectors]
    singular = np.linalg.svd(np.concatenate([d_uv.real, d_uv.imag], axis=0), compute_uv=False)
    tangent_tolerance = float(np.sqrt(tolerance / 4.0))
    tangent_rank = int(np.sum(singular > tangent_tolerance))
    checks = {
        "singlet_endpoint_tangents_are_opposite": bool(np.linalg.norm(d_u + d_v) < 3e-9),
        "two_sided_qfi_matches_lifted_edge_qfi": bool(np.allclose(two_qfi, expected, atol=3e-9, rtol=3e-9)),
        "every_edge_common_endpoint_rotation_is_null": bool(max(kernel_residuals, default=0.0) < 3e-9),
        "two_sided_qfi_is_psd_up_to_roundoff": bool(eigenvalues[0] >= -1e-9 * scale),
        "two_sided_rank_matches_edge_qfi_rank": two_rank == edge_rank,
        "tangent_rank_matches_two_sided_qfi_rank": tangent_rank == two_rank,
        "edge_qfi_symmetric": bool(edge_symmetry < 1e-10),
    }
    return {
        "gluing_probability": probability,
        "edge_qfi_shape": list(edge_qfi.shape),
        "two_sided_qfi_shape": list(two_qfi.shape),
        "edge_qfi_rank": edge_rank,
        "two_sided_qfi_rank": two_rank,
        "two_sided_qfi_nullity": 6 * edge_count - two_rank,
        "structural_common_rotation_nullity": 3 * edge_count,
        "additional_nullity_from_boundary_stabilizers": 6 * edge_count - two_rank - 3 * edge_count,
        "two_sided_minimum_eigenvalue": float(eigenvalues[0]),
        "max_endpoint_tangent_sum_residual": float(np.linalg.norm(d_u + d_v)),
        "max_common_rotation_qfi_residual": max(kernel_residuals, default=0.0),
        "max_lifted_qfi_error": float(np.linalg.norm(two_qfi - expected)),
        "max_tangent_singular_value": float(singular[0]) if singular.size else 0.0,
        "tangent_rank_tolerance": tangent_tolerance,
        "checks": checks,
    }


def run(seed: int = 20261118, trials: int = 6) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    results = {}
    for name, (node_count, edges) in graph_fixtures().items():
        rows = [
            analyze_sample([random_intertwiner(rng) for _ in range(node_count)], edges)
            for _ in range(trials)
        ]
        checks = {
            key: all(row["checks"][key] for row in rows)
            for key in rows[0]["checks"]
        }
        results[name] = {
            "nodes": node_count,
            "edges": len(edges),
            "trials": trials,
            "edge_qfi_rank_min": min(row["edge_qfi_rank"] for row in rows),
            "two_sided_qfi_rank_min": min(row["two_sided_qfi_rank"] for row in rows),
            "two_sided_qfi_nullity_max": max(row["two_sided_qfi_nullity"] for row in rows),
            "max_endpoint_tangent_sum_residual": max(row["max_endpoint_tangent_sum_residual"] for row in rows),
            "max_lifted_qfi_error": max(row["max_lifted_qfi_error"] for row in rows),
            "checks": checks,
            "passed": all(checks.values()),
            "samples": rows,
        }
    return {
        "status": "two-sided SU(2) singlet endpoint QFI quotient control",
        "passed": all(row["passed"] for row in results.values()),
        "seed": seed,
        "trials_per_graph": trials,
        "endpoint_convention": "same Pauli generators on both tensor indices of normalized epsilon singlet",
        "exact_structure": "D_v(e)=-D_u(e); F_two-sided is the signed lift [[F,-F],[-F,F]]",
        "quotient": "identify endpoint pairs along common-rotation null directions; relative triplets reproduce edge QFI",
        "claim_boundary": (
            "This is an exact kinematic consequence of the chosen invariant singlet "
            "contraction and the tested finite spin-1/2 graph states. It does not "
            "establish endpoint-frame gauge dynamics, arbitrary representations, "
            "BQG history amplitudes, or continuum conductances."
        ),
        "graphs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=6)
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
