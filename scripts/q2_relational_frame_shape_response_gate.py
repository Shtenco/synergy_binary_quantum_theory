#!/usr/bin/env python3
"""Finite response diagnostic for frame and logical-intertwiner variations.

Uses the same path_3/cycle_4 contraction as T75. In addition to the exact
frame-tangent construction, it differentiates multilinearly with respect to
two local logical-ray directions per node. The resulting rank check is
floating-point and exploratory; it is not an exact combined minor certificate.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.linalg import svdvals

import q2_relational_graph_qfi_gate as graph_gate
from q2_relational_gauge_quotient_gate import projective_tangent_jacobian
from q2_relational_joint_gauge_edge_quotient_gate import joint_quotient_basis
from spatial_qubit_geometry_gate import logical_basis


WITNESSES = {
    "path_3": [[1, -1 + 2j], [-1 + 2j, 2 - 2j], [-1, 1 - 1j]],
    # This cycle witness also keeps every node away from the zero-volume
    # boundary of the tetrahedral shape disk.
    "cycle_4": [[1 + 1j, 2], [1, 2j], [2, -1 - 1j], [1, 1 + 1j]],
}


def analyze(name: str) -> dict[str, object]:
    node_count, edges = graph_gate.graph_fixtures()[name]
    logical = logical_basis()
    coefficients = [np.asarray(row, dtype=complex) for row in WITNESSES[name]]
    rays = [row / np.linalg.norm(row) for row in coefficients]
    states = [logical @ ray for ray in rays]

    frame_tangent, _probability = projective_tangent_jacobian(states, edges)
    frame_quotient, dimensions = joint_quotient_basis(node_count, edges)
    frame_tangent = frame_tangent @ frame_quotient
    raw, _ = graph_gate.contract_graph(states, edges)
    raw_norm = np.linalg.norm(raw)
    ray = raw / raw_norm

    # For normalized c=(a,b), u=(-conj(b),conj(a)) is a normalized complex
    # orthogonal tangent. The real directions u and i*u span the local CP^1
    # tangent plane. Contracting with u is the exact first derivative because
    # the graph contraction is multilinear in every node tensor.
    shape_columns = []
    for node, coefficients_at_node in enumerate(rays):
        a, b = coefficients_at_node
        tangent = np.array([-np.conjugate(b), np.conjugate(a)])
        for direction in (tangent, 1j * tangent):
            varied = list(states)
            varied[node] = logical @ direction
            derivative, _ = graph_gate.contract_graph(varied, edges)
            derivative = derivative / raw_norm
            derivative -= ray * np.vdot(ray, derivative)
            shape_columns.append(derivative)

    shape_tangent = np.column_stack(shape_columns)
    frame_real = np.vstack((frame_tangent.real, frame_tangent.imag))
    shape_real = np.vstack((shape_tangent.real, shape_tangent.imag))
    combined = np.column_stack((frame_real, shape_real))
    threshold = 1e-8
    frame_singular = svdvals(frame_real)
    shape_singular = svdvals(shape_real)
    combined_singular = svdvals(combined)
    frame_rank = int(np.sum(frame_singular > threshold))
    shape_rank = int(np.sum(shape_singular > threshold))
    combined_rank = int(np.sum(combined_singular > threshold))
    expected = dimensions["joint_quotient_dimension"] + 2 * node_count
    return {
        "fixture": name,
        "frame_quotient_dimension": dimensions["joint_quotient_dimension"],
        "frame_tangent_rank": frame_rank,
        "node_shape_dimension": 2 * node_count,
        "node_shape_tangent_rank": shape_rank,
        "combined_rank": combined_rank,
        "combined_expected_rank": expected,
        "minimum_shape_singular_value": float(shape_singular[-1]),
        "minimum_combined_singular_value": float(combined_singular[-1]),
        "rank_threshold": threshold,
        "derivative_method": "exact multilinear substitution; floating-point contraction and SVD",
        "interpretation": "finite response diagnostic only; no exact combined-minor or refinement-gap certificate",
    }


def main() -> int:
    result = {"fixtures": [analyze(name) for name in WITNESSES]}
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

