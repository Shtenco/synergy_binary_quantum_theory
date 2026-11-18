#!/usr/bin/env python3
"""Exact pushforward of independent Haar four-valent intertwiners under vertex fusion.

This is a deliberately specified toy ensemble, not a claim about the physical
BQG history measure. It tests whether independent local Haar measures close
to a Haar measure on the full invariant space of a fused six-valent block.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.q2_relational_higher_spin_graph_qfi_gate import (
    invariant_singlet_tensor,
    spin_generators,
)
from scripts.q2_relational_higher_spin_network_qfi_gate import (
    four_valent_intertwiner_basis,
)
from scripts.q2_vertex_fusion_coarse_graining import fuse_vertex_blocks


def singlet_multiplicity(two_j: int, leg_count: int) -> int:
    """Multiplicity of total spin zero in a tensor product of equal spins."""
    multiplicities = {0: 1}
    for _ in range(leg_count):
        updated: dict[int, int] = {}
        for total, multiplicity in multiplicities.items():
            for output in range(abs(total - two_j), total + two_j + 1, 2):
                updated[output] = updated.get(output, 0) + multiplicity
        multiplicities = updated
    return int(multiplicities.get(0, 0))


def apply_generator_axis(
    tensor: np.ndarray, axis: int, generator: np.ndarray
) -> np.ndarray:
    return np.moveaxis(np.tensordot(generator, tensor, axes=(1, axis)), 0, axis)


def analyze_independent_haar_fusion(two_j: int = 1) -> dict[str, object]:
    if isinstance(two_j, bool) or not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer")
    dimension = two_j + 1
    generators = spin_generators(two_j)
    singlet = invariant_singlet_tensor(two_j)
    local_basis, local_gauss_residual = four_valent_intertwiner_basis(two_j)

    # Two vertices are joined through leg 3 of vertex 0 and leg 0 of vertex 1.
    edge = [(0, 3, 1, 0)]
    columns = []
    for left in range(local_basis.shape[1]):
        for right in range(local_basis.shape[1]):
            fused, boundaries, coarse_edges = fuse_vertex_blocks(
                [local_basis[:, left], local_basis[:, right]],
                edge,
                singlet,
                ((0, 1),),
            )
            if coarse_edges or len(boundaries) != 1 or len(boundaries[0]) != 6:
                raise AssertionError("unexpected two-vertex fusion boundary")
            columns.append(np.asarray(fused[0], complex).reshape(-1))
    contraction_map = np.column_stack(columns)

    gram = contraction_map.conj().T @ contraction_map
    expected_gram = np.eye(gram.shape[0], dtype=complex) / (dimension**2)
    gram_error = float(np.max(np.abs(gram - expected_gram)))
    singular_values = np.linalg.svd(contraction_map, compute_uv=False)
    map_rank = int(np.sum(singular_values > 1e-10))

    # Check the six-leg Gauss law without materializing a dense d^6 by d^6
    # generator matrix. This remains practical through spin j=2.
    max_six_valent_gauss_residual = 0.0
    for column in columns:
        tensor = column.reshape((dimension,) * 6)
        norm = float(np.linalg.norm(tensor))
        for generator in generators:
            total = np.zeros_like(tensor)
            for axis in range(6):
                total += apply_generator_axis(tensor, axis, generator)
            max_six_valent_gauss_residual = max(
                max_six_valent_gauss_residual,
                float(np.linalg.norm(total)) / max(norm, 1e-14),
            )

    six_valent_invariant_dimension = singlet_multiplicity(two_j, 6)
    input_dimension = int(local_basis.shape[1] ** 2)
    pushforward_covariance_eigenvalues = singular_values**2
    positive = pushforward_covariance_eigenvalues[
        pushforward_covariance_eigenvalues > 1e-10
    ]
    covariance_trace = float(np.sum(positive))
    covariance_purity = float(np.sum(positive**2))
    haar_purity = 1.0 / six_valent_invariant_dimension
    non_haar_distance = float(np.sqrt(
        1.0 / input_dimension - 1.0 / six_valent_invariant_dimension
    ))

    return {
        "status": "exact pushforward of independent local Haar intertwiners",
        "ensemble": (
            "two independent Haar-random unit vectors in the four-valent "
            "intertwiner space, contracted through one invariant singlet"
        ),
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "local_intertwiner_dimension": int(local_basis.shape[1]),
        "product_input_dimension": input_dimension,
        "fused_boundary_leg_count": 6,
        "full_six_valent_invariant_dimension": six_valent_invariant_dimension,
        "fusion_map_rank": map_rank,
        "fusion_map_singular_values": singular_values.tolist(),
        "fusion_map_gram_max_abs_error_from_I_over_d_squared": gram_error,
        "local_gauss_residual": float(local_gauss_residual),
        "six_valent_gauss_max_relative_residual": max_six_valent_gauss_residual,
        "pushforward_covariance_trace": covariance_trace,
        "pushforward_covariance_rank": int(positive.size),
        "pushforward_nonzero_covariance_eigenvalues": positive.tolist(),
        "pushforward_covariance_purity": covariance_purity,
        "full_six_valent_haar_covariance_purity": haar_purity,
        "pushforward_vs_full_haar_covariance_frobenius_distance": non_haar_distance,
        "passed": (
            local_basis.shape[1] == dimension
            and map_rank == input_dimension
            and six_valent_invariant_dimension > input_dimension
            and gram_error < 2e-9
            and local_gauss_residual < 2e-9
            and max_six_valent_gauss_residual < 2e-9
            and abs(covariance_trace - 1.0) < 2e-9
            and positive.size == input_dimension
            and np.max(np.abs(positive - 1.0 / input_dimension)) < 2e-9
            and abs(covariance_purity - 1.0 / input_dimension) < 2e-9
            and abs(haar_purity - 1.0 / six_valent_invariant_dimension) < 2e-9
            and non_haar_distance > 1e-3
        ),
        "interpretation": (
            "For this explicitly chosen independent local-Haar ensemble, "
            "two-vertex fusion maps the product input isometrically up to "
            "scale into a proper subspace of the six-valent invariant space. "
            "Its exact pushforward covariance is therefore not the full "
            "six-valent Haar covariance. This demonstrates non-closure of "
            "this toy measure, not a statement about the unspecified physical "
            "BQG measure."
        ),
        "claim_boundary": (
            "Exact finite-dimensional result for the listed spin, two "
            "four-valent vertices, one singlet edge, and independent local "
            "Haar measures only."
        ),
    }


def analyze_spin_family(two_j_values=(1, 2, 3, 4)) -> dict[str, object]:
    controls = {
        str(two_j / 2.0): analyze_independent_haar_fusion(two_j)
        for two_j in two_j_values
    }
    return {
        "status": "higher-spin exact Haar pushforward under two-vertex fusion",
        "tested_spins": [two_j / 2.0 for two_j in two_j_values],
        "passed": all(row["passed"] for row in controls.values()),
        "controls": controls,
        "claim_boundary": (
            "Exact pushforward of independent local Haar intertwiner measures "
            "for the listed finite spins and one two-vertex block only; not "
            "the physical BQG measure or a completed RG flow."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-j", type=int, action="append")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = (
        analyze_independent_haar_fusion(args.two_j[0])
        if args.two_j and len(args.two_j) == 1
        else analyze_spin_family(tuple(args.two_j) if args.two_j else (1, 2, 3, 4))
    )
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
