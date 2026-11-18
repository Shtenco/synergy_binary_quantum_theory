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
        "path_4": (4, [(0, 3, 1, 0), (1, 3, 2, 0), (2, 3, 3, 0)]),
        "path_5": (5, [(0, 3, 1, 0), (1, 3, 2, 0), (2, 3, 3, 0), (3, 3, 4, 0)]),
        "star_4": (4, [(0, 0, 1, 0), (0, 1, 2, 0), (0, 2, 3, 0)]),
        "cycle_4": (4, [(0, 0, 1, 1), (1, 2, 2, 1),
                       (2, 2, 3, 1), (3, 2, 0, 2)]),
    }


def _validated_edge(edge, node_count: int) -> tuple[int, int, int, int]:
    if len(edge) != 4:
        raise ValueError("edges must be (u, leg_u, v, leg_v)")
    if any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        for value in edge
    ):
        raise ValueError("edge endpoint and half-link indices must be integers")
    u, lu, v, lv = map(int, edge)
    if u == v or not (0 <= u < node_count and 0 <= v < node_count):
        raise ValueError("invalid graph edge endpoints")
    if not (0 <= lu < 4 and 0 <= lv < 4):
        raise ValueError("invalid half-link slots")
    return u, lu, v, lv


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
        u, lu, v, lv = _validated_edge(edge, len(states))
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


def collective_edge_blocks(name, edge_count):
    """A fixed, explicit partition used to restrict edge twists to blockwise modes."""
    partitions = {
        "path_3": ((0, 1),),
        "path_4": ((0, 1), (2,)),
        "path_5": ((0, 1), (2, 3)),
        "star_4": ((0, 1), (2,)),
        "cycle_4": ((0, 1), (2, 3)),
    }
    if name not in partitions:
        raise ValueError(f"no collective edge partition for graph {name!r}")
    blocks = partitions[name]
    flat = [edge for block in blocks for edge in block]
    if (
        any(
            isinstance(edge, (bool, np.bool_))
            or not isinstance(edge, (int, np.integer))
            for edge in flat
        )
        or any(len(block) == 0 for block in blocks)
        or sorted(map(int, flat)) != list(range(edge_count))
    ):
        raise ValueError("collective edge blocks must partition every edge exactly once")
    return tuple(tuple(map(int, block)) for block in blocks)


def collective_edge_parameter_map(edge_count, blocks, orthonormalize_blocks=False):
    """Embed block angles; optionally give each collective mode unit Euclidean norm."""
    flat = [edge for block in blocks for edge in block]
    if (
        any(
            isinstance(edge, (bool, np.bool_))
            or not isinstance(edge, (int, np.integer))
            for edge in flat
        )
        or any(len(block) == 0 for block in blocks)
        or sorted(map(int, flat)) != list(range(edge_count))
    ):
        raise ValueError("collective edge blocks must partition every edge exactly once")
    mapping = np.zeros((3 * edge_count, 3 * len(blocks)), dtype=float)
    for block_index, block in enumerate(blocks):
        for edge in block:
            coefficient = 1.0 / math.sqrt(len(block)) if orthonormalize_blocks else 1.0
            mapping[3 * int(edge):3 * int(edge) + 3,
                    3 * block_index:3 * block_index + 3] = coefficient * np.eye(3)
    return mapping


def collective_edge_pullback_qfi(micro_qfi, blocks, orthonormalize_blocks=False):
    """Pull back QFI using shared-angle or unit-norm collective block coordinates."""
    if np.asarray(micro_qfi).ndim != 2 or micro_qfi.shape[0] != micro_qfi.shape[1]:
        raise ValueError("micro QFI must be a square matrix")
    if micro_qfi.shape[0] % 3:
        raise ValueError("micro QFI dimension must be divisible by three")
    edge_count = micro_qfi.shape[0] // 3
    mapping = collective_edge_parameter_map(
        edge_count, blocks, orthonormalize_blocks=orthonormalize_blocks
    )
    return mapping.T @ np.asarray(micro_qfi, dtype=float) @ mapping, mapping


def analyze_graph(two_j, name, node_count, edges, trials=3, seed=20261118):
    rng = np.random.default_rng(seed + two_j + len(edges))
    basis, _ = four_valent_intertwiner_basis(two_j)
    singlet = invariant_singlet_tensor(two_j)
    generators = spin_generators(two_j)
    j = two_j / 2.0
    isolated_qfi = 4.0 * j * (j + 1.0) / 3.0
    probabilities, symmetries, minima, cross_norms, diagonal_errors = [], [], [], [], []
    spectral_minima, spectral_maxima, spectral_ranks, spectral_nullities = [], [], [], []
    path_casimir_errors = []
    normalized_moments = {str(order): [] for order in (1, 2, 3, 4)}
    normalized_quantiles = []
    normalized_variances = []
    rank_fractions = []
    block_spectral_moments = {str(order): [] for order in (1, 2, 3, 4)}
    block_spectral_quantiles = []
    block_spectral_variances = []
    block_rank_fractions = []
    block_symmetry_errors = []
    block_min_eigenvalues = []
    orthonormal_block_moments = {str(order): [] for order in (1, 2, 3, 4)}
    orthonormal_block_quantiles = []
    orthonormal_block_variances = []
    orthonormal_block_rank_fractions = []
    orthonormal_block_min_eigenvalues = []
    orthonormal_block_symmetry_errors = []
    blocks = collective_edge_blocks(name, len(edges))
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
        eigenvalues = np.linalg.eigvalsh(qfi)
        minima.append(float(eigenvalues[0]))
        spectral_scale = max(1.0, float(np.max(np.abs(eigenvalues))))
        spectral_tolerance = 1e-9 * spectral_scale
        spectral_minima.append(float(eigenvalues[0]))
        spectral_maxima.append(float(eigenvalues[-1]))
        spectral_ranks.append(int(np.sum(eigenvalues > spectral_tolerance)))
        spectral_nullities.append(int(np.sum(np.abs(eigenvalues) <= spectral_tolerance)))
        normalized_eigenvalues = eigenvalues / isolated_qfi
        for order in (1, 2, 3, 4):
            normalized_moments[str(order)].append(
                float(np.mean(normalized_eigenvalues ** order))
            )
        normalized_quantiles.append(
            np.quantile(normalized_eigenvalues, [0.0, 0.25, 0.5, 0.75, 1.0]).tolist()
        )
        normalized_variances.append(float(np.var(normalized_eigenvalues)))
        rank_fractions.append(float(spectral_ranks[-1] / (3 * len(edges))))
        effective_qfi, _block_map = collective_edge_pullback_qfi(qfi, blocks)
        effective_symmetry_error = float(np.linalg.norm(effective_qfi - effective_qfi.T))
        effective_eigenvalues = np.linalg.eigvalsh(
            0.5 * (effective_qfi + effective_qfi.T)
        )
        normalized_effective = effective_eigenvalues / isolated_qfi
        for order in (1, 2, 3, 4):
            block_spectral_moments[str(order)].append(
                float(np.mean(normalized_effective ** order))
            )
        block_spectral_quantiles.append(
            np.quantile(normalized_effective, [0.0, 0.25, 0.5, 0.75, 1.0]).tolist()
        )
        block_spectral_variances.append(float(np.var(normalized_effective)))
        block_rank_fractions.append(float(
            np.sum(effective_eigenvalues > 1e-9 * max(1.0, float(np.max(np.abs(effective_eigenvalues)))))
            / effective_qfi.shape[0]
        ))
        block_symmetry_errors.append(effective_symmetry_error)
        block_min_eigenvalues.append(float(effective_eigenvalues[0]))
        unit_qfi, _unit_map = collective_edge_pullback_qfi(
            qfi, blocks, orthonormalize_blocks=True
        )
        unit_symmetry_error = float(np.linalg.norm(unit_qfi - unit_qfi.T))
        unit_eigenvalues = np.linalg.eigvalsh(0.5 * (unit_qfi + unit_qfi.T))
        normalized_unit = unit_eigenvalues / isolated_qfi
        for order in (1, 2, 3, 4):
            orthonormal_block_moments[str(order)].append(
                float(np.mean(normalized_unit ** order))
            )
        orthonormal_block_quantiles.append(
            np.quantile(normalized_unit, [0.0, 0.25, 0.5, 0.75, 1.0]).tolist()
        )
        orthonormal_block_variances.append(float(np.var(normalized_unit)))
        orthonormal_block_rank_fractions.append(float(
            np.sum(unit_eigenvalues > 1e-9 * max(1.0, float(np.max(np.abs(unit_eigenvalues)))))
            / unit_qfi.shape[0]
        ))
        orthonormal_block_min_eigenvalues.append(float(unit_eigenvalues[0]))
        orthonormal_block_symmetry_errors.append(unit_symmetry_error)
        block_diagonal = np.zeros_like(qfi)
        for edge_index in range(len(edges)):
            sl = slice(3 * edge_index, 3 * edge_index + 3)
            block_diagonal[sl, sl] = qfi[sl, sl]
            diagonal_errors.append(float(np.linalg.norm(qfi[sl, sl] - isolated_qfi * np.eye(3))))
        cross_norms.append(float(np.linalg.norm(qfi - block_diagonal)))
        if name == "path_3":
            predicted = -(4.0 / 3.0) * float(
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
        "qfi_spectral_minimum_over_trials": min(spectral_minima),
        "qfi_spectral_maximum_over_trials": max(spectral_maxima),
        "qfi_rank_min_over_trials": min(spectral_ranks),
        "qfi_rank_max_over_trials": max(spectral_ranks),
        "qfi_nullity_min_over_trials": min(spectral_nullities),
        "qfi_nullity_max_over_trials": max(spectral_nullities),
        "qfi_parameter_dimension": 3 * len(edges),
        "qfi_rank_fraction_mean": float(np.mean(rank_fractions)),
        "qfi_spectral_normalization": "eigenvalues divided by isolated-link scale 4*j*(j+1)/3",
        "normalized_spectral_moments_mean": {
            order: float(np.mean(values)) for order, values in normalized_moments.items()
        },
        "normalized_spectral_quantiles_mean": np.mean(
            np.asarray(normalized_quantiles), axis=0
        ).tolist(),
        "normalized_spectral_quantiles_probability_levels": [0.0, 0.25, 0.5, 0.75, 1.0],
        "normalized_spectral_variance_mean": float(np.mean(normalized_variances)),
        "max_diagonal_block_error_from_isolated_spin_qfi": max(diagonal_errors),
        "cross_edge_qfi_frobenius_norm_min": min(cross_norms),
        "cross_edge_qfi_frobenius_norm_max": max(cross_norms),
        "max_path_middle_casimir_formula_error": max(path_casimir_errors) if path_casimir_errors else None,
        "collective_blocking": {
            "edge_partition": [list(block) for block in blocks],
            "parameter_map": "theta_micro[e] = theta_block[b] for e in block b; B has 3x3 identity blocks",
            "effective_parameter_dimension": 3 * len(blocks),
            "pullback": "F_effective = B.T @ F_micro @ B",
            "normalized_spectral_moments_mean": {
                order: float(np.mean(values)) for order, values in block_spectral_moments.items()
            },
            "normalized_spectral_quantiles_mean": np.mean(
                np.asarray(block_spectral_quantiles), axis=0
            ).tolist(),
            "normalized_spectral_variance_mean": float(np.mean(block_spectral_variances)),
            "rank_fraction_mean": float(np.mean(block_rank_fractions)),
            "minimum_eigenvalue_over_trials": min(block_min_eigenvalues),
            "max_symmetry_error": max(block_symmetry_errors),
            "checks": {
                "effective_qfi_is_psd": min(block_min_eigenvalues) > -2e-8,
                "effective_qfi_is_symmetric": max(block_symmetry_errors) < 2e-9,
            },
            "orthonormal_collective_mode": {
                "parameter_map": "theta_micro[e] = theta_unit[b] / sqrt(|b|); columns of B are orthonormal",
                "normalized_spectral_moments_mean": {
                    order: float(np.mean(values)) for order, values in orthonormal_block_moments.items()
                },
                "normalized_spectral_quantiles_mean": np.mean(
                    np.asarray(orthonormal_block_quantiles), axis=0
                ).tolist(),
                "normalized_spectral_variance_mean": float(np.mean(orthonormal_block_variances)),
                "rank_fraction_mean": float(np.mean(orthonormal_block_rank_fractions)),
                "minimum_eigenvalue_over_trials": min(orthonormal_block_min_eigenvalues),
                "max_symmetry_error": max(orthonormal_block_symmetry_errors),
                "checks": {
                    "effective_qfi_is_psd": min(orthonormal_block_min_eigenvalues) > -2e-8,
                    "effective_qfi_is_symmetric": max(orthonormal_block_symmetry_errors) < 2e-9,
                },
            },
        },
        "checks": {
            "qfi_is_symmetric": max(symmetries) < 2e-9,
            "qfi_is_psd": min(minima) > -2e-8,
            "collective_pullback_is_psd": min(block_min_eigenvalues) > -2e-8,
            "collective_pullback_is_symmetric": max(block_symmetry_errors) < 2e-9,
            "qfi_rank_nullity_partition_is_consistent": all(
                rank + nullity == 3 * len(edges)
                for rank, nullity in zip(spectral_ranks, spectral_nullities)
            ),
            "cross_edge_correlations_are_present": max(cross_norms) > 1e-5,
            "path_diagonal_blocks_match_isolated_link": (
                max(diagonal_errors) < 2e-8 if name == "path_3" else True
            ),
            "path_cross_block_matches_middle_casimir": (
                max(path_casimir_errors) < 2e-8 if name == "path_3" else True
            ),
        },
    }


def run(two_j_values=(1, 2), trials=3, graph_names=None):
    available_graphs = graph_fixtures()
    selected_names = tuple(available_graphs) if graph_names is None else tuple(graph_names)
    unknown = set(selected_names) - set(available_graphs)
    if unknown:
        raise ValueError(f"unknown graph fixtures: {sorted(unknown)}")
    rows = {}
    for two_j in two_j_values:
        rows[str(two_j / 2.0)] = {}
        for name in selected_names:
            nodes, edges = available_graphs[name]
            # The path_5 boundary has 12 open legs. Keep its dense boundary
            # vector to j<=1; higher spins remain covered on <=10-leg fixtures.
            if name == "path_5" and two_j > 2:
                continue
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
        "graph_spin_coverage": {
            name: [float(two_j / 2.0) for two_j in two_j_values
                   if name in rows[str(two_j / 2.0)]]
            for name in selected_names
        },
        "interpretation": (
            "All spectra are also reported after division by the exact isolated-link "
            "Casimir scale 4*j*(j+1)/3. Moments, quantiles, rank fraction, and "
            "variance provide dimensionless diagnostics; they are not asserted "
            "to be coarse-graining invariants. Collective blocking is an explicit "
            "parameter restriction, not vertex fusion or boundary tracing. "
            "Path_5 is restricted to j<=1 because "
            "its 12-leg dense boundary vector grows rapidly with representation "
            "dimension. Cycle and branched graph metrics are computed directly, "
            "without path additivity assumptions."
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
