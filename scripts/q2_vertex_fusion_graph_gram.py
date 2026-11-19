#!/usr/bin/env python3
"""Finite graph-block fusion Gram diagnostics for tree and loop contractions.

The input is a product of four-valent SU(2)-invariant tensors. All edges of
the selected graph are contracted inside one block. The Gram matrix C^dagger C
distinguishes tree-like isometric fusion from loop-induced rank loss and
state-dependent normalization. This is a kinematic finite control, not a
physical BQG measure or RG flow.
"""
from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from scripts.q2_relational_higher_spin_graph_qfi_gate import (
    graph_fixtures,
    invariant_singlet_tensor,
    spin_generators,
)
from scripts.q2_relational_higher_spin_network_qfi_gate import (
    four_valent_intertwiner_basis,
)
from scripts.q2_vertex_fusion_coarse_graining import (
    fuse_vertex_blocks,
    fused_block_gauss_residuals,
)


def singlet_multiplicity(two_j: int, leg_count: int) -> int:
    """Multiplicity of total spin zero in leg_count copies of spin two_j/2."""
    multiplicities = {0: 1}
    for _ in range(leg_count):
        updated: dict[int, int] = {}
        for total, multiplicity in multiplicities.items():
            for output in range(abs(total - two_j), total + two_j + 1, 2):
                updated[output] = updated.get(output, 0) + multiplicity
        multiplicities = updated
    return int(multiplicities.get(0, 0))


def graph_cycle_rank(node_count: int, edges) -> int:
    """Return E - V + number_of_connected_components."""
    adjacency = [set() for _ in range(node_count)]
    for u, _lu, v, _lv in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    seen = set()
    components = 0
    for root in range(node_count):
        if root in seen:
            continue
        components += 1
        stack = [root]
        seen.add(root)
        while stack:
            node = stack.pop()
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
    return len(edges) - node_count + components


def _relative_singular_rank(singular_values: np.ndarray, tolerance: float) -> int:
    """Rank estimate invariant under an overall rescaling of the contraction map."""
    values = np.asarray(singular_values, dtype=float)
    if values.size == 0 or values[0] <= 0.0:
        return 0
    return int(np.count_nonzero(values > tolerance * values[0]))


def _positive_spectrum(eigenvalues: np.ndarray, tolerance: float) -> np.ndarray:
    """Select positive eigenvalues using a scale-relative spectral threshold."""
    values = np.asarray(eigenvalues, dtype=float)
    if values.size == 0:
        return values
    scale = float(np.max(np.abs(values)))
    if scale == 0.0:
        return values[:0]
    return values[values > tolerance * scale]


def _relative_threshold(scale: float, tolerance: float) -> float:
    """Return a nonzero relative tolerance without imposing a unit scale."""
    return tolerance * max(abs(float(scale)), np.finfo(float).tiny)


def analyze_graph_block_fusion_gram(
    graph_name: str = "star_4",
    two_j: int = 1,
    rank_tolerance: float = 1e-10,
    scalar_tolerance: float = 2e-9,
) -> dict[str, object]:
    """Build the exact finite contraction matrix and diagnose its Gram spectrum."""
    if graph_name not in graph_fixtures():
        raise ValueError(f"unknown graph fixture: {graph_name}")
    if isinstance(two_j, bool) or not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer")
    if rank_tolerance <= 0 or scalar_tolerance <= 0:
        raise ValueError("tolerances must be positive")

    node_count, edges = graph_fixtures()[graph_name]
    dimension = two_j + 1
    local_basis, local_gauss_residual = four_valent_intertwiner_basis(two_j)
    local_dimension = int(local_basis.shape[1])
    singlet = invariant_singlet_tensor(two_j)
    generators = spin_generators(two_j)
    blocks = (tuple(range(node_count)),)
    columns = []
    max_gauss_residual = 0.0

    for local_indices in product(range(local_dimension), repeat=node_count):
        states = [
            local_basis[:, index].copy()
            for index in local_indices
        ]
        fused, boundaries, coarse_edges = fuse_vertex_blocks(
            states, edges, singlet, blocks
        )
        if len(fused) != 1 or len(boundaries) != 1 or coarse_edges:
            raise AssertionError("single-block fusion did not absorb every internal edge")
        expected_boundary_legs = 4 * node_count - 2 * len(edges)
        if len(boundaries[0]) != expected_boundary_legs:
            raise AssertionError("fused boundary leg count is inconsistent with graph")
        max_gauss_residual = max(
            max_gauss_residual,
            fused_block_gauss_residuals(fused, generators)[0],
        )
        columns.append(np.asarray(fused[0], complex).reshape(-1))

    contraction = np.column_stack(columns)
    gram = contraction.conj().T @ contraction
    gram = 0.5 * (gram + gram.conj().T)
    eigenvalues = np.linalg.eigvalsh(gram)
    singular_values = np.linalg.svd(contraction, compute_uv=False)

    input_dimension = local_dimension ** node_count
    boundary_leg_count = 4 * node_count - 2 * len(edges)
    boundary_invariant_dimension = singlet_multiplicity(two_j, boundary_leg_count)
    rank_threshold = _relative_threshold(
        float(singular_values[0]) if singular_values.size else 0.0,
        rank_tolerance,
    )
    rank = _relative_singular_rank(singular_values, rank_tolerance)
    positive_eigenvalues = _positive_spectrum(eigenvalues, rank_tolerance)
    gram_scalar = float(np.trace(gram).real / input_dimension)
    gram_scalar_error = float(np.max(np.abs(
        gram - gram_scalar * np.eye(input_dimension, dtype=complex)
    )))
    gram_is_scalar = gram_scalar_error <= _relative_threshold(gram_scalar, scalar_tolerance)
    cycle_rank = graph_cycle_rank(node_count, edges)
    tree_expected_scale = (
        float(dimension ** (-2 * len(edges))) if cycle_rank == 0 else None
    )
    tree_expected_gram_error = (
        float(np.max(np.abs(
            gram - tree_expected_scale * np.eye(input_dimension, dtype=complex)
        )))
        if tree_expected_scale is not None else None
    )
    minimum_eigenvalue = float(eigenvalues[0]) if eigenvalues.size else 0.0
    maximum_eigenvalue = float(eigenvalues[-1]) if eigenvalues.size else 0.0
    positive_singular_values = singular_values[singular_values > rank_threshold]
    minimum_positive_singular_value = (
        float(positive_singular_values[-1]) if positive_singular_values.size else 0.0
    )
    maximum_singular_value = (
        float(positive_singular_values[0]) if positive_singular_values.size else 0.0
    )
    condition_number_on_support = (
        maximum_singular_value / minimum_positive_singular_value
        if minimum_positive_singular_value > 0 else None
    )

    passed = (
        rank > 0
        and rank <= min(input_dimension, boundary_invariant_dimension)
        and minimum_eigenvalue >= -_relative_threshold(maximum_eigenvalue, scalar_tolerance)
        and local_gauss_residual < scalar_tolerance
        and max_gauss_residual < scalar_tolerance
        and (
            tree_expected_gram_error is None
            or tree_expected_gram_error < _relative_threshold(tree_expected_scale, scalar_tolerance)
        )
    )
    if cycle_rank == 0:
        interpretation = (
            "The tested tree block has C^dagger C proportional to identity. "
            "Its contraction norm is input-state independent, and the "
            "normalized fusion map is an isometry up to the reported scale."
        )
    else:
        interpretation = (
            "The tested loop block has a non-scalar, rank-deficient Gram map. "
            "Some microscopic intertwiner directions are annihilated and the "
            "norm of the fused tensor depends on the microscopic input. "
            "A normalized pushforward ensemble therefore requires the "
            "state-dependent factor 1/||C psi||^2; the tree covariance formula "
            "cannot be reused."
        )

    return {
        "status": "finite graph-block fusion Gram diagnostic",
        "graph": graph_name,
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "node_count": int(node_count),
        "edge_count": int(len(edges)),
        "cycle_rank": int(cycle_rank),
        "boundary_leg_count": int(boundary_leg_count),
        "boundary_invariant_dimension": int(boundary_invariant_dimension),
        "local_intertwiner_dimension": int(local_dimension),
        "input_dimension": int(input_dimension),
        "output_vector_dimension": int(contraction.shape[0]),
        "fusion_map_rank": int(rank),
        "fusion_map_rank_defect": int(input_dimension - rank),
        "singular_value_max": maximum_singular_value,
        "singular_value_min_positive": minimum_positive_singular_value,
        "singular_value_condition_number_on_support": condition_number_on_support,
        "gram_eigenvalue_min": minimum_eigenvalue,
        "gram_eigenvalue_max": maximum_eigenvalue,
        "gram_scalar_mean": gram_scalar,
        "gram_max_abs_error_from_scalar_identity": gram_scalar_error,
        "gram_is_scalar_identity": bool(gram_is_scalar),
        "tree_expected_gram_scale": tree_expected_scale,
        "tree_expected_gram_max_abs_error": tree_expected_gram_error,
        "local_gauss_residual": float(local_gauss_residual),
        "max_fused_boundary_gauss_relative_residual": float(max_gauss_residual),
        "passed": bool(passed),
        "interpretation": interpretation,
        "claim_boundary": (
            "Finite contraction matrix for the named graph fixture and spin. "
            "No physical history measure, spin sum, boundary trace, effective "
            "action, continuum limit, or physical RG flow is derived."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default="star_4", choices=sorted(graph_fixtures()))
    parser.add_argument("--two-j", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze_graph_block_fusion_gram(args.graph, args.two_j)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
