#!/usr/bin/env python3
"""Joint quotient by vertex Gauss rotations and singlet-edge common rotations.

The full half-link QFI has two structurally distinct exact null families:
(1) one SU(2) Gauss rotation on all four legs at each node; and
(2) one common same-axis rotation of the two half-links contracted by each
normalized invariant singlet. This gate constructs both generator families,
checks their combined rank, and tests whether any projective stabilizers
remain after quotienting by their joint span.
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
    full_half_link_qfi,
    graph_fixtures,
    random_intertwiner,
)
from q2_relational_gauge_quotient_gate import (
    projective_tangent_jacobian,
    vertex_gauge_generators,
)


def edge_singlet_common_generators(
    node_count: int, edges: list[tuple[int, int, int, int]]
) -> np.ndarray:
    """One 12V-vector per internal edge/axis, +1 on both contracted legs."""
    if node_count < 1:
        raise ValueError("node_count must be positive")
    result = np.zeros((12 * node_count, 3 * len(edges)), dtype=float)
    occupied = set()
    for edge_index, edge in enumerate(edges):
        if len(edge) != 4:
            raise ValueError("edge must be (u, leg_u, v, leg_v)")
        u, lu, v, lv = map(int, edge)
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("invalid edge endpoints")
        if not (0 <= lu < 4 and 0 <= lv < 4):
            raise ValueError("invalid half-link indices")
        for endpoint in ((u, lu), (v, lv)):
            if endpoint in occupied:
                raise ValueError("a half-link cannot belong to multiple edges")
            occupied.add(endpoint)
        for axis in range(3):
            result[12 * u + 3 * lu + axis, 3 * edge_index + axis] = 1.0
            result[12 * v + 3 * lv + axis, 3 * edge_index + axis] = 1.0
    return result


def joint_quotient_basis(
    node_count: int, edges: list[tuple[int, int, int, int]]
) -> tuple[np.ndarray, dict[str, int]]:
    gauss = vertex_gauge_generators(node_count)
    edge_common = edge_singlet_common_generators(node_count, edges)
    generators = np.column_stack([gauss, edge_common])
    u, singular, _vh = np.linalg.svd(generators, full_matrices=True)
    scale = float(singular[0]) if singular.size else 1.0
    tolerance = max(generators.shape) * np.finfo(float).eps * max(1.0, scale) * 100.0
    rank = int(np.sum(singular > tolerance))
    quotient = u[:, rank:]
    return quotient, {
        "full_parameter_count": 12 * node_count,
        "vertex_gauss_generator_count": 3 * node_count,
        "edge_singlet_common_generator_count": 3 * len(edges),
        "joint_generator_rank": rank,
        "joint_quotient_dimension": int(quotient.shape[1]),
        "expected_independent_generator_count": 3 * node_count + 3 * len(edges),
        "generator_families_independent": bool(
            rank == 3 * node_count + 3 * len(edges)
        ),
    }


def analyze_sample(states, edges):
    node_count = len(states)
    _probability, qfi, symmetry, _boundary = full_half_link_qfi(states, edges)
    tangent, _ = projective_tangent_jacobian(states, edges)
    gauss = vertex_gauge_generators(node_count)
    edge_common = edge_singlet_common_generators(node_count, edges)
    quotient, dimensions = joint_quotient_basis(node_count, edges)
    gauss_residual = float(np.linalg.norm(qfi @ gauss))
    edge_residual = float(np.linalg.norm(qfi @ edge_common))
    joint_residual = float(np.linalg.norm(qfi @ np.column_stack([gauss, edge_common])))
    quotient_qfi = quotient.T @ qfi @ quotient
    quotient_qfi = 0.5 * (quotient_qfi + quotient_qfi.T)
    eigenvalues = np.linalg.eigvalsh(quotient_qfi)
    scale = max(1.0, float(np.max(np.abs(eigenvalues))))
    qfi_tolerance = 1e-9 * scale
    qfi_rank = int(np.sum(eigenvalues > qfi_tolerance))
    qfi_nullity = int(np.sum(np.abs(eigenvalues) <= qfi_tolerance))
    quotient_tangent = tangent @ quotient
    realified = np.concatenate([quotient_tangent.real, quotient_tangent.imag], axis=0)
    singular = np.linalg.svd(realified, compute_uv=False)
    tangent_tolerance = float(np.sqrt(qfi_tolerance / 4.0))
    tangent_rank = int(np.sum(singular > tangent_tolerance))
    checks = {
        "vertex_gauss_generators_are_qfi_null": bool(gauss_residual < 2e-9),
        "edge_singlet_common_generators_are_qfi_null": bool(edge_residual < 2e-9),
        "joint_null_generators_are_qfi_null": bool(joint_residual < 3e-9),
        "gauss_and_edge_generators_are_independent": dimensions["generator_families_independent"],
        "joint_quotient_qfi_is_psd": bool(eigenvalues[0] > -2e-9 * scale),
        "joint_quotient_has_no_remaining_null_modes": bool(qfi_nullity == 0),
        "joint_tangent_and_qfi_ranks_agree": bool(tangent_rank == qfi_rank),
        "joint_quotient_is_full_rank": bool(qfi_rank == quotient.shape[1]),
        "full_half_link_qfi_is_symmetric": bool(symmetry < 1e-9),
    }
    return {
        **dimensions,
        "gauss_annihilation_residual": gauss_residual,
        "edge_common_annihilation_residual": edge_residual,
        "joint_annihilation_residual": joint_residual,
        "joint_quotient_minimum_eigenvalue": float(eigenvalues[0]),
        "joint_quotient_maximum_eigenvalue": float(eigenvalues[-1]),
        "joint_quotient_qfi_rank": qfi_rank,
        "joint_quotient_qfi_nullity": qfi_nullity,
        "joint_quotient_tangent_rank": tangent_rank,
        "tangent_svd_tolerance": tangent_tolerance,
        "checks": checks,
        "passed": all(checks.values()),
    }


def run(seed=20261118, trials=4):
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    graphs = {}
    for name, (node_count, edges) in graph_fixtures().items():
        samples = [
            analyze_sample([random_intertwiner(rng) for _ in range(node_count)], edges)
            for _ in range(trials)
        ]
        graphs[name] = {
            "nodes": node_count,
            "internal_edges": len(edges),
            "trials": trials,
            "joint_quotient_dimension": samples[0]["joint_quotient_dimension"],
            "joint_generator_rank": samples[0]["joint_generator_rank"],
            "minimum_joint_quotient_eigenvalue": min(
                row["joint_quotient_minimum_eigenvalue"] for row in samples
            ),
            "maximum_joint_quotient_nullity": max(
                row["joint_quotient_qfi_nullity"] for row in samples
            ),
            "maximum_gauss_annihilation_residual": max(
                row["gauss_annihilation_residual"] for row in samples
            ),
            "maximum_edge_common_annihilation_residual": max(
                row["edge_common_annihilation_residual"] for row in samples
            ),
            "maximum_joint_annihilation_residual": max(
                row["joint_annihilation_residual"] for row in samples
            ),
            "all_checks_pass": all(row["passed"] for row in samples),
            "samples": samples,
        }
    return {
        "status": "joint vertex-Gauss and singlet-edge common-rotation quotient",
        "passed": all(row["all_checks_pass"] for row in graphs.values()),
        "seed": seed,
        "graphs": graphs,
        "interpretation": (
            "The full half-link QFI has exact null generators from nodewise "
            "Gauss rotations and from simultaneous same-axis rotations of the "
            "two endpoints of each normalized invariant singlet. The quotient "
            "by their joint span is tested for residual projective stabilizers."
        ),
        "claim_boundary": (
            "This classifies exact nulls of the present finite spin-1/2 "
            "contraction. It does not establish that the singlet-edge nulls "
            "are an independent fundamental gauge symmetry, nor does it "
            "derive a general BQG edge-frame bundle or continuum dynamics."
        ),
    }


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=4)
    args = parser.parse_args()
    result = run(trials=args.trials)
    rendered = json.dumps(result, indent=2, default=_json_default)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
