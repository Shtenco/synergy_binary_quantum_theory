#!/usr/bin/env python3
"""General finite SU(2) spin-network graph QFI from explicit singlet contractions.

Each vertex carries a normalized four-valent spin-1/2 intertwiner. Each graph
edge contracts one half-link at each endpoint with the normalized epsilon
singlet. Remaining half-links form the boundary. One independent relative
twist parameter triplet is inserted on the first endpoint of each graph edge.
The complete pure-boundary-state multiparameter QFI is computed without
assuming edge additivity.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import spatial_qubit_geometry_gate as SQ

EPSILON = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex) / math.sqrt(2.0)
PAULI = SQ.PAULI


def node_generator(pauli: np.ndarray, leg: int) -> np.ndarray:
    return SQ.kron_all([
        0.5 * pauli if index == leg else SQ.I2
        for index in range(4)
    ])


def validate_graph(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
) -> list[tuple[int, int]]:
    if not node_states:
        raise ValueError("at least one node state is required")
    occupied: set[tuple[int, int]] = set()
    for node, state in enumerate(node_states):
        vector = np.asarray(state, dtype=complex)
        if vector.shape != (16,):
            raise ValueError(f"node {node} must be a 16-component four-qubit state")
        if not np.isclose(np.linalg.norm(vector), 1.0, atol=1e-11):
            raise ValueError(f"node {node} state must be normalized")
    normalized_edges = []
    for edge in edges:
        if len(edge) != 4:
            raise ValueError("each edge must be (u, leg_u, v, leg_v)")
        u, lu, v, lv = map(int, edge)
        if u == v or not (0 <= u < len(node_states)) or not (0 <= v < len(node_states)):
            raise ValueError("edge endpoints must be distinct valid node indices")
        if not (0 <= lu < 4 and 0 <= lv < 4):
            raise ValueError("half-link slots must be in [0,3]")
        for endpoint in ((u, lu), (v, lv)):
            if endpoint in occupied:
                raise ValueError("a half-link cannot belong to more than one graph edge")
            occupied.add(endpoint)
        normalized_edges.append((u, lu, v, lv))
    return normalized_edges


def contract_graph(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
    twist: tuple[int, int, np.ndarray] | None = None,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Contract all graph edges; optional twist=(edge_index, axis, generator)."""
    edges = validate_graph(node_states, edges)
    tensors = [np.asarray(state, dtype=complex).reshape((2,) * 4).copy()
               for state in node_states]
    if twist is not None:
        edge_index, _axis, generator = twist
        if not 0 <= edge_index < len(edges):
            raise ValueError("twist edge index is out of range")
        u, leg_u, _v, _leg_v = edges[edge_index]
        tensors[u] = (-1j * (node_generator(generator, leg_u) @ tensors[u].reshape(16))).reshape((2,) * 4)

    labels = [[None] * 4 for _ in node_states]
    next_label = 0
    epsilon_tensors = []
    for u, lu, v, lv in edges:
        label_a, label_b = next_label, next_label + 1
        next_label += 2
        labels[u][lu] = label_a
        labels[v][lv] = label_b
        epsilon_tensors.append((EPSILON, [label_a, label_b]))
    boundary = []
    for node in range(len(node_states)):
        for leg in range(4):
            if labels[node][leg] is None:
                labels[node][leg] = next_label
                boundary.append((node, leg))
                next_label += 1
    operands = []
    for node, tensor in enumerate(tensors):
        operands.extend([tensor, labels[node]])
    for tensor, tensor_labels in epsilon_tensors:
        operands.extend([tensor, tensor_labels])
    operands.append([labels[node][leg] for node, leg in boundary])
    amplitude = np.einsum(*operands, optimize=True).reshape(-1)
    return amplitude, boundary


def graph_qfi(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
) -> tuple[float, np.ndarray, float, list[tuple[int, int]]]:
    edges = validate_graph(node_states, edges)
    raw, boundary = contract_graph(node_states, edges)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("graph singlet contraction has zero support")
    state = raw / math.sqrt(probability)
    derivatives = []
    for edge_index in range(len(edges)):
        for axis, pauli in enumerate(PAULI):
            derivative_raw, derivative_boundary = contract_graph(
                node_states, edges, (edge_index, axis, pauli)
            )
            if derivative_boundary != boundary:
                raise AssertionError("boundary ordering changed under twist")
            radial = float(np.vdot(raw, derivative_raw).real) / probability
            derivatives.append(derivative_raw / math.sqrt(probability) - radial * state)
    overlaps = [np.vdot(state, derivative) for derivative in derivatives]
    qfi = np.zeros((3 * len(edges), 3 * len(edges)), dtype=float)
    for i, j in itertools.product(range(len(derivatives)), repeat=2):
        qfi[i, j] = 4.0 * float(np.real(
            np.vdot(derivatives[i], derivatives[j])
            - np.conjugate(overlaps[i]) * overlaps[j]
        ))
    return probability, qfi, float(np.linalg.norm(qfi - qfi.T)), boundary


def random_intertwiner(rng: np.random.Generator) -> np.ndarray:
    coefficients = rng.normal(size=2) + 1j * rng.normal(size=2)
    coefficients /= np.linalg.norm(coefficients)
    return SQ.logical_basis() @ coefficients


def graph_fixtures():
    # Three-node path: two adjacent edges share a middle intertwiner.
    path_edges = [(0, 3, 1, 0), (1, 3, 2, 0)]
    # Four-cycle: each node uses two half-links internally and retains two.
    cycle_edges = [(0, 0, 1, 1), (1, 2, 2, 1),
                   (2, 2, 3, 1), (3, 2, 0, 2)]
    return {"path_3": (3, path_edges), "cycle_4": (4, cycle_edges)}


def run(seed: int = 20261118, trials: int = 12) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    rows = {}
    for name, (node_count, edges) in graph_fixtures().items():
        probabilities = []
        symmetry = []
        minimum_eigenvalues = []
        cross_norms = []
        diagonal_errors = []
        for _ in range(trials):
            states = [random_intertwiner(rng) for _ in range(node_count)]
            probability, qfi, sym_error, _boundary = graph_qfi(states, edges)
            probabilities.append(probability)
            symmetry.append(sym_error)
            minimum_eigenvalues.append(float(np.linalg.eigvalsh(qfi)[0]))
            block_diagonal = np.zeros_like(qfi)
            for edge in range(len(edges)):
                sl = slice(3 * edge, 3 * edge + 3)
                block_diagonal[sl, sl] = qfi[sl, sl]
                diagonal_errors.append(float(np.linalg.norm(qfi[sl, sl] - np.eye(3))))
            cross_norms.append(float(np.linalg.norm(qfi - block_diagonal)))
        rows[name] = {
            "nodes": node_count,
            "edges": len(edges),
            "qfi_shape": [3 * len(edges), 3 * len(edges)],
            "boundary_qubits": 4 * node_count - 2 * len(edges),
            "gluing_probability_min": min(probabilities),
            "gluing_probability_max": max(probabilities),
            "max_symmetry_error": max(symmetry),
            "minimum_qfi_eigenvalue": min(minimum_eigenvalues),
            "max_diagonal_link_block_error_from_identity": max(diagonal_errors),
            "cross_link_qfi_frobenius_norm_min": min(cross_norms),
            "cross_link_qfi_frobenius_norm_max": max(cross_norms),
            "checks": {
                "qfi_symmetric": max(symmetry) < 1e-10,
                "qfi_positive_semidefinite": min(minimum_eigenvalues) > -1e-10,
                "diagonal_link_blocks_identity": max(diagonal_errors) < 1e-9,
                "cross_link_terms_present": max(cross_norms) > 1e-3,
            },
        }
        rows[name]["passed"] = all(rows[name]["checks"].values())
    return {
        "status": "general finite spin-network graph relational QFI control",
        "passed": all(row["passed"] for row in rows.values()),
        "seed": seed,
        "trials_per_graph": trials,
        "edge_parameterization": "one relative SU(2) twist triplet per graph edge, inserted at the first endpoint",
        "state_construction": "explicit contraction of normalized four-valent SU(2) intertwiners with normalized epsilon singlets",
        "graphs": rows,
        "claim_boundary": (
            "Only the explicit path and four-cycle finite spin-1/2 graph "
            "families are tested. This does not yet define arbitrary-graph "
            "boundary gauge reduction, higher-spin sums, BQG history measure, "
            "causal-PL coarse graining, or continuum spectral dimension."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=12)
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
