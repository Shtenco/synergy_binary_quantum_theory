#!/usr/bin/env python3
"""Finite SU(2) holonomy quotient: endpoint coordinates vs vertex gauge.

For an oriented edge e=(u,v), g_e -> h_u g_e h_v^{-1}.
With anti-Hermitian generators T_a=-i sigma_a/2 and left-trivialized
variation, delta g g^{-1} = x_u - Ad_{g_e} x_v.

Six endpoint coordinates map to three tangent coordinates. Their kernel is
an endpoint-coordinate redundancy (two descriptions of the same delta g),
not an additional physical edge-frame gauge symmetry. Vertex gauge is a
separate graph-wide orbit. The quotient dimension depends on holonomy
stabilizers and graph cycle structure.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)
GENERATORS = tuple(-0.5j * p for p in PAULI)


def su2_from_vector(vector: np.ndarray) -> np.ndarray:
    """Exponential exp(x^a T_a), with T_a=-i sigma_a/2."""
    x = np.asarray(vector, dtype=float)
    if x.shape != (3,):
        raise ValueError("SU(2) exponential coordinate must have shape (3,)")
    angle = float(np.linalg.norm(x))
    if angle < 1e-14:
        return np.eye(2, dtype=complex) + sum(
            x[a] * GENERATORS[a] for a in range(3)
        )
    generator = sum(x[a] * GENERATORS[a] for a in range(3))
    return math.cos(angle / 2.0) * np.eye(2) + (
        2.0 * math.sin(angle / 2.0) / angle
    ) * generator


def adjoint_matrix(g: np.ndarray) -> np.ndarray:
    """Return Ad_g in the stated anti-Hermitian basis."""
    matrix = np.asarray(g, dtype=complex)
    if matrix.shape != (2, 2):
        raise ValueError("holonomy must be a 2x2 matrix")
    if not np.allclose(matrix.conj().T @ matrix, np.eye(2), atol=1e-10):
        raise ValueError("holonomy must be unitary")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1e-10):
        raise ValueError("holonomy must have determinant one")
    result = np.zeros((3, 3), dtype=float)
    for a, ta in enumerate(GENERATORS):
        for b, tb in enumerate(GENERATORS):
            result[a, b] = float(np.real(-2.0 * np.trace(ta @ matrix @ tb @ matrix.conj().T)))
    if not np.allclose(result.T @ result, np.eye(3), atol=2e-10):
        raise AssertionError("SU(2) adjoint representation is not orthogonal")
    return result


def endpoint_to_holonomy_jacobian(holonomies: list[np.ndarray]) -> np.ndarray:
    """Block diagonal H: R^(6E) endpoint parameters -> R^(3E) link tangents."""
    result = np.zeros((3 * len(holonomies), 6 * len(holonomies)), dtype=float)
    for edge, g in enumerate(holonomies):
        result[3*edge:3*edge+3, 6*edge:6*edge+3] = np.eye(3)
        result[3*edge:3*edge+3, 6*edge+3:6*edge+6] = -adjoint_matrix(g)
    return result


def vertex_gauge_embedding(
    node_count: int, edges: list[tuple[int, int]]
) -> np.ndarray:
    """Embed vertex Lie-algebra parameters into 6E endpoint coordinates."""
    if node_count < 1:
        raise ValueError("node_count must be positive")
    result = np.zeros((6 * len(edges), 3 * node_count), dtype=float)
    for edge, (u, v) in enumerate(edges):
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("edge must join two distinct valid vertices")
        result[6*edge:6*edge+3, 3*u:3*u+3] = np.eye(3)
        result[6*edge+3:6*edge+6, 3*v:3*v+3] = np.eye(3)
    return result


def holonomy_gauge_jacobian(
    node_count: int,
    edges: list[tuple[int, int]],
    holonomies: list[np.ndarray],
) -> np.ndarray:
    if len(edges) != len(holonomies):
        raise ValueError("one holonomy is required for every edge")
    return endpoint_to_holonomy_jacobian(holonomies) @ vertex_gauge_embedding(
        node_count, edges
    )


def numerical_rank(matrix: np.ndarray, relative_tolerance: float = 1e-9) -> tuple[int, np.ndarray, float]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    scale = float(singular[0]) if singular.size else 0.0
    tolerance = relative_tolerance * max(1.0, scale)
    return int(np.sum(singular > tolerance)), singular, tolerance


def transform_holonomies(
    holonomies: list[np.ndarray], edges: list[tuple[int, int]],
    vertex_transforms: list[np.ndarray],
) -> list[np.ndarray]:
    return [
        vertex_transforms[u] @ g @ vertex_transforms[v].conj().T
        for g, (u, v) in zip(holonomies, edges)
    ]


def loop_holonomy(holonomies: list[np.ndarray]) -> np.ndarray:
    """Ordered product for a directed closed cycle."""
    result = np.eye(2, dtype=complex)
    for g in holonomies:
        result = result @ g
    return result


def graph_fixtures():
    return {
        "path_3": (3, [(0, 1), (1, 2)]),
        "cycle_4": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    }


def fixture_holonomies(name: str) -> list[np.ndarray]:
    if name == "path_3":
        vectors = [(0.31, -0.17, 0.22), (-0.19, 0.28, 0.13)]
    elif name == "cycle_4":
        # Non-commuting generic holonomies; their product has a U(1) stabilizer.
        vectors = [
            (0.31, -0.17, 0.22), (-0.19, 0.28, 0.13),
            (0.11, 0.23, -0.29), (-0.27, -0.08, 0.34),
        ]
    else:
        raise KeyError(name)
    return [su2_from_vector(np.asarray(vector)) for vector in vectors]


def diagnostics(
    name: str, node_count: int, edges: list[tuple[int, int]],
    holonomies: list[np.ndarray],
) -> dict[str, object]:
    endpoint = endpoint_to_holonomy_jacobian(holonomies)
    embedding = vertex_gauge_embedding(node_count, edges)
    gauge = endpoint @ embedding
    endpoint_rank, endpoint_singular, _ = numerical_rank(endpoint)
    gauge_rank, gauge_singular, _ = numerical_rank(gauge)
    gauge_stabilizer = 3 * node_count - gauge_rank
    physical_tangent_dimension = 3 * len(edges) - gauge_rank
    loop_count = len(edges) - node_count + 1
    return {
        "graph": name,
        "nodes": node_count,
        "edges": len(edges),
        "cycle_rank": loop_count,
        "endpoint_coordinate_count": 6 * len(edges),
        "link_holonomy_tangent_dimension": 3 * len(edges),
        "endpoint_map_rank": endpoint_rank,
        "endpoint_map_kernel_dimension": 6 * len(edges) - endpoint_rank,
        "endpoint_kernel_expected": 3 * len(edges),
        "vertex_gauge_parameter_count": 3 * node_count,
        "vertex_gauge_orbit_rank": gauge_rank,
        "vertex_gauge_stabilizer_dimension": gauge_stabilizer,
        "holonomy_mod_vertex_gauge_tangent_dimension": physical_tangent_dimension,
        "endpoint_map_singular_values": endpoint_singular.tolist(),
        "vertex_gauge_singular_values": gauge_singular.tolist(),
        "checks": {
            "six_to_three_endpoint_reduction": endpoint_rank == 3 * len(edges),
            "endpoint_kernel_is_three_per_edge": 6 * len(edges) - endpoint_rank == 3 * len(edges),
            "gauge_orbit_within_holonomy_space": gauge_rank <= 3 * len(edges),
            "quotient_dimension_nonnegative": physical_tangent_dimension >= 0,
        },
    }


def run() -> dict[str, object]:
    results = {}
    for name, (node_count, edges) in graph_fixtures().items():
        holonomies = fixture_holonomies(name)
        row = diagnostics(name, node_count, edges, holonomies)
        results[name] = row

    # Finite covariance and infinitesimal derivative test on the cycle.
    name = "cycle_4"
    node_count, edges = graph_fixtures()[name]
    holonomies = fixture_holonomies(name)
    vertex_x = [
        np.array([0.13, -0.07, 0.19]),
        np.array([-0.11, 0.17, 0.05]),
        np.array([0.09, 0.04, -0.16]),
        np.array([-0.06, -0.12, 0.08]),
    ]
    transforms = [su2_from_vector(x) for x in vertex_x]
    transformed = transform_holonomies(holonomies, edges, transforms)
    old_loop = loop_holonomy(holonomies)
    new_loop = loop_holonomy(transformed)
    expected_loop = transforms[0] @ old_loop @ transforms[0].conj().T

    eps = 1e-7
    endpoint_direction = vertex_gauge_embedding(node_count, edges) @ np.concatenate(vertex_x)
    predicted = holonomy_gauge_jacobian(node_count, edges, holonomies) @ np.concatenate(vertex_x)
    finite_differences = []
    for edge, (u, v) in enumerate(edges):
        varied = (
            su2_from_vector(eps * vertex_x[u]) @ holonomies[edge]
            @ su2_from_vector(-eps * vertex_x[v])
        )
        delta = (varied - holonomies[edge]) / eps @ holonomies[edge].conj().T
        measured = np.array([
            float(np.real(-2.0 * np.trace(GENERATORS[a] @ delta)))
            for a in range(3)
        ])
        finite_differences.extend(measured.tolist())
    measured = np.asarray(finite_differences)
    finite_difference_error = float(np.linalg.norm(measured - predicted))
    loop_covariance_error = float(np.linalg.norm(new_loop - expected_loop))

    # A generic SU(2) holonomy has a one-dimensional Lie-algebra centralizer.
    cycle_gauge_stabilizer = results["cycle_4"]["vertex_gauge_stabilizer_dimension"]
    tree_gauge_stabilizer = results["path_3"]["vertex_gauge_stabilizer_dimension"]
    results["cycle_4"]["loop_trace"] = float(np.real(np.trace(old_loop)))
    results["cycle_4"]["loop_conjugacy_class_tangent_dimension"] = (
        results["cycle_4"]["holonomy_mod_vertex_gauge_tangent_dimension"]
    )

    checks = {
        "finite_vertex_gauge_covariance": loop_covariance_error < 2e-10,
        "infinitesimal_endpoint_law": finite_difference_error < 2e-7,
        "tree_has_no_local_holonomy_moduli": results["path_3"]["holonomy_mod_vertex_gauge_tangent_dimension"] == 0,
        "generic_single_cycle_has_one_conjugacy_modulus": results["cycle_4"]["holonomy_mod_vertex_gauge_tangent_dimension"] == 1,
        "tree_global_stabilizer_dimension": tree_gauge_stabilizer == 3,
        "generic_cycle_stabilizer_dimension": cycle_gauge_stabilizer == 1,
        "all_endpoint_maps_have_expected_kernel": all(
            row["checks"]["endpoint_kernel_is_three_per_edge"] for row in results.values()
        ),
    }
    return {
        "status": "finite SU(2) holonomy and vertex-gauge quotient control",
        "passed": all(checks.values()),
        "conventions": {
            "group_action": "g_e -> h_u g_e h_v^{-1}",
            "lie_algebra_basis": "T_a = -i sigma_a/2",
            "left_trivialized_variation": "delta g_e g_e^{-1} = xi_u - Ad_{g_e} xi_v",
            "endpoint_coordinates": "g_e(t) = exp(t x_{e,u}) g_e exp(-t x_{e,v})",
            "endpoint_to_holonomy_map": "[I_3, -Ad_{g_e}]",
            "vertex_gauge_embedding": "x_{e,u}=xi_u and x_{e,v}=xi_v for every incidence",
            "endpoint_kernel": "x_{e,u}=Ad_{g_e} x_{e,v}; a coordinate redundancy, not an independent physical gauge symmetry",
            "physical_moduli": "holonomy tangent space modulo the image of the vertex gauge action",
        },
        "finite_difference_error": finite_difference_error,
        "loop_covariance_error": loop_covariance_error,
        "graphs": results,
        "checks": checks,
        "claim_boundary": (
            "These are finite-dimensional SU(2) kinematics on a path and a cycle. "
            "No independent edge-frame field has been introduced, so no additional "
            "edge-frame gauge symmetry is asserted. No BQG dynamics, QFI weighting, "
            "boundary-state stabilizer identification, or continuum limit follows."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
