#!/usr/bin/env python3
"""Higher-spin finite-graph relational QFI for path and cycle contractions."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_relational_higher_spin_network_qfi_gate import (
    four_valent_intertwiner_basis,
    invariant_singlet_tensor,
    node_generator,
    random_intertwiner,
    spin_generators,
)


def graph_fixtures():
    return {
        "path_3": (3, [(0, 3, 1, 0), (1, 3, 2, 0)]),
        "cycle_4": (4, [(0, 0, 1, 1), (1, 2, 2, 1),
                       (2, 2, 3, 1), (3, 2, 0, 2)]),
    }


def validate_graph(states, edges, dimension):
    if not states:
        raise ValueError("at least one node state is required")
    occupied = set()
    for node, state in enumerate(states):
        if np.asarray(state).shape != (dimension**4,):
            raise ValueError(f"node {node} has an invalid four-valent state shape")
        if not np.isclose(np.linalg.norm(state), 1.0, atol=1e-10):
            raise ValueError(f"node {node} state is not normalized")
    normalized = []
    for edge in edges:
        if len(edge) != 4:
            raise ValueError("edges must be (u, leg_u, v, leg_v)")
        u, lu, v, lv = map(int, edge)
        if u == v or not (0 <= u < len(states) and 0 <= v < len(states)):
            raise ValueError("invalid graph edge endpoints")
        if not (0 <= lu < 4 and 0 <= lv < 4):
            raise ValueError("invalid half-link slots")
        for endpoint in ((u, lu), (v, lv)):
            if endpoint in occupied:
                raise ValueError("a half-link is used by multiple edges")
            occupied.add(endpoint)
        normalized.append((u, lu, v, lv))
    return normalized


def contract_graph(states, edges, singlet, twist=None):
    dimension = singlet.shape[0]
    edges = validate_graph(states, edges, dimension)
    tensors = [np.asarray(state, complex).reshape((dimension,) * 4).copy()
               for state in states]
    if twist is not None:
        edge_index, generator = twist
        if not 0 <= edge_index < len(edges):
            raise ValueError("twist edge index is out of range")
        node, leg, _v, _other_leg = edges[edge_index]
        tensors[node] = (
            -1j * (node_generator(generator, leg, dimension) @ tensors[node].reshape(-1))
        ).reshape((dimension,) * 4)
    labels = [[None] * 4 for _ in states]
    next_label = 0
    operands = []
    for node, tensor in enumerate(tensors):
        operands.extend([tensor, labels[node]])
    for u, lu, v, lv in edges:
        left, right = next_label, next_label + 1
        next_label += 2
        labels[u][lu], labels[v][lv] = left, right
        operands.extend([singlet, [left, right]])
    boundary = []
    for node in range(len(states)):
        for leg in range(4):
            if labels[node][leg] is None:
                labels[node][leg] = next_label
                boundary.append((node, leg))
                next_label += 1
    # Rebuild node operands after all labels have been assigned.
    operands = []
    for node, tensor in enumerate(tensors):
        operands.extend([tensor, labels[node]])
    for edge_index, (u, lu, v, lv) in enumerate(edges):
        # Each edge's pair labels are recoverable from the node label table.
        operands.extend([singlet, [labels[u][lu], labels[v][lv]]])
    operands.append([labels[node][leg] for node, leg in boundary])
    return np.einsum(*operands, optimize=True).reshape(-1), boundary


def graph_qfi(states, edges, singlet, generators):
    dimension = singlet.shape[0]
    raw, boundary = contract_graph(states, edges, singlet)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("graph contraction has zero support")
    state = raw / math.sqrt(probability)
    derivatives = []
    for edge_index in range(len(edges)):
        for generator in generators:
            derivative_raw, derivative_boundary = contract_graph(
                states, edges, singlet, (edge_index, generator)
            )
            if derivative_boundary != boundary:
                raise AssertionError("boundary ordering changed under twist")
            derivative = derivative_raw / math.sqrt(probability)
            derivative -= state * np.vdot(state, derivative)
            derivatives.append(derivative)
    qfi = np.zeros((3 * len(edges), 3 * len(edges)), dtype=float)
    for i in range(len(derivatives)):
        for k in range(len(derivatives)):
            qfi[i, k] = 4.0 * float(np.real(np.vdot(derivatives[i], derivatives[k])))
    return probability, qfi, float(np.linalg.norm(qfi - qfi.T)), boundary


def analyze_graph(two_j, name, node_count, edges, trials=3, seed=20261118):
    rng = np.random.default_rng(seed + two_j + len(edges))
    basis, _ = four_valent_intertwiner_basis(two_j)
    singlet = invariant_singlet_tensor(two_j)
    generators = spin_generators(two_j)
    j = two_j / 2.0
    isolated_qfi = 4.0 * j * (j + 1.0) / 3.0
    probabilities, symmetries, minima, cross_norms, diagonal_errors = [], [], [], [], []
    path_casimir_errors = []
    pair_op = sum(
        (node_generator(g, 0, two_j + 1) @ node_generator(g, 3, two_j + 1)
         for g in generators),
        np.zeros((two_j + 1)**4, dtype=complex),
    )
    for _ in range(trials):
        states = [random_intertwiner(rng, basis) for _ in range(node_count)]
        probability, qfi, symmetry, _boundary = graph_qfi(states, edges, singlet, generators)
        probabilities.append(probability)
        symmetries.append(symmetry)
        minima.append(float(np.linalg.eigvalsh(qfi)[0]))
        block_diagonal = np.zeros_like(qfi)
        for edge_index in range(len(edges)):
            sl = slice(3 * edge_index, 3 * edge_index + 3)
            block_diagonal[sl, sl] = qfi[sl, sl]
            diagonal_errors.append(float(np.linalg.norm(qfi[sl, sl] - isolated_qfi * np.eye(3))))
        cross_norms.append(float(np.linalg.norm(qfi - block_diagonal)))
        if name == "path_3":
            predicted = (4.0 / 3.0) * float(
                np.vdot(states[1], pair_op @ states[1]).real
            )
            path_casimir_errors.append(float(np.linalg.norm(
                qfi[:3, 3:] - predicted * np.eye(3)
            )))
    return {
        "nodes": node_count,
        "edges": len(edges),
        "qfi_shape": [3 * len(edges), 3 * len(edges)],
        "boundary_half_links": 4 * node_count - 2 * len(edges),
        "gluing_probability_min": min(probabilities),
        "gluing_probability_max": max(probabilities),
        "max_qfi_symmetry_error": max(symmetries),
        "minimum_qfi_eigenvalue": min(minima),
        "max_diagonal_block_error_from_isolated_spin_qfi": max(diagonal_errors),
        "cross_edge_qfi_frobenius_norm_min": min(cross_norms),
        "cross_edge_qfi_frobenius_norm_max": max(cross_norms),
        "max_path_middle_casimir_formula_error": max(path_casimir_errors) if path_casimir_errors else None,
        "checks": {
            "qfi_is_symmetric": max(symmetries) < 2e-9,
            "qfi_is_psd": min(minima) > -2e-8,
            "cross_edge_correlations_are_present": max(cross_norms) > 1e-5,
            "path_diagonal_blocks_match_isolated_link": (
                max(diagonal_errors) < 2e-8 if name == "path_3" else True
            ),
            "path_cross_block_matches_middle_casimir": (
                max(path_casimir_errors) < 2e-8 if name == "path_3" else True
            ),
        },
    }


def run(two_j_values=(1, 2), trials=3):
    rows = {}
    for two_j in two_j_values:
        rows[str(two_j / 2.0)] = {}
        for name, (nodes, edges) in graph_fixtures().items():
            rows[str(two_j / 2.0)][name] = analyze_graph(
                two_j, name, nodes, edges, trials=trials
            )
    return {
        "status": "higher-spin finite graph relational QFI",
        "passed": all(
            all(result["checks"].values())
            for spin_rows in rows.values() for result in spin_rows.values()
        ),
        "tested_spins": [int(value) / 2.0 for value in two_j_values],
        "trials_per_graph_and_spin": trials,
        "graphs": rows,
        "interpretation": (
            "Path and cycle graph contractions are evaluated directly. The path "
            "has the exact middle-node Casimir cross block and isolated-link "
            "diagonal scale. Cycle metrics are not assumed additive or equal to "
            "the path; their complete QFI is computed from the contracted state."
        ),
        "claim_boundary": (
            "Finite pure-state higher-spin graph contractions only. No spin sum, "
            "boundary tracing, BQG history amplitudes, or continuum graph weight."
        ),
    }


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=3)
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
