#!/usr/bin/env python3
"""Exact pushforward of independent Haar four-valent intertwiners under vertex fusion.

This is a deliberately specified toy ensemble, not a claim about the physical
BQG history measure. It tests whether independent local Haar measures close
to a Haar measure on the full invariant space of a fused six-valent block.
"""
from __future__ import annotations

import argparse
import json
from itertools import product
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



def product_vs_global_haar_two_copy_moment(two_j: int, vertex_count: int) -> dict[str, float]:
    """Exact Hilbert-Schmidt comparison of product-Haar and global-Haar 2-copy moments."""
    dimension = two_j + 1
    input_dimension = dimension**vertex_count
    product_purity = (
        2.0 / (dimension * (dimension + 1.0))
    ) ** vertex_count
    global_haar_purity = 2.0 / (
        input_dimension * (input_dimension + 1.0)
    )
    squared_distance = product_purity - global_haar_purity
    if squared_distance < -1e-12:
        raise AssertionError("two-copy moment distance squared became negative")
    return {
        "product_haar_two_copy_moment_purity": product_purity,
        "global_haar_on_input_two_copy_moment_purity": global_haar_purity,
        "two_copy_moment_hilbert_schmidt_distance": float(
            np.sqrt(max(0.0, squared_distance))
        ),
    }


def analyze_independent_haar_fusion(
    two_j: int = 1, vertex_count: int = 2
) -> dict[str, object]:
    if isinstance(two_j, bool) or not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer")
    if isinstance(vertex_count, bool) or not isinstance(vertex_count, int) or vertex_count < 2:
        raise ValueError("vertex_count must be an integer of at least two")
    dimension = two_j + 1
    generators = spin_generators(two_j)
    singlet = invariant_singlet_tensor(two_j)
    local_basis, local_gauss_residual = four_valent_intertwiner_basis(two_j)

    # A path block joins leg 3 of each vertex to leg 0 of its successor.
    edges = [(vertex, 3, vertex + 1, 0) for vertex in range(vertex_count - 1)]
    block = (tuple(range(vertex_count)),)
    columns = []
    for local_indices in product(range(local_basis.shape[1]), repeat=vertex_count):
        states = [local_basis[:, index] for index in local_indices]
        fused, boundaries, coarse_edges = fuse_vertex_blocks(
            states, edges, singlet, block
        )
        expected_boundary_legs = 2 * vertex_count + 2
        if (
            coarse_edges
            or len(boundaries) != 1
            or len(boundaries[0]) != expected_boundary_legs
        ):
            raise AssertionError("unexpected fused chain boundary")
        columns.append(np.asarray(fused[0], complex).reshape(-1))
    contraction_map = np.column_stack(columns)

    input_dimension = int(local_basis.shape[1] ** vertex_count)
    product_haar_projective_support_dimension = vertex_count * (
        local_basis.shape[1] - 1
    )
    full_projective_image_dimension = input_dimension - 1
    segre_codimension_in_image = (
        full_projective_image_dimension - product_haar_projective_support_dimension
    )
    gram = contraction_map.conj().T @ contraction_map
    expected_gram = np.eye(gram.shape[0], dtype=complex) / (
        dimension ** (2 * (vertex_count - 1))
    )
    gram_error = float(np.max(np.abs(gram - expected_gram)))
    singular_values = np.linalg.svd(contraction_map, compute_uv=False)
    map_rank = int(np.sum(singular_values > 1e-10))

    # Check the boundary Gauss law without materializing a dense d^L by d^L
    # generator matrix. This remains practical for the listed chain blocks.
    boundary_leg_count = 2 * vertex_count + 2
    max_boundary_gauss_residual = 0.0
    for column in columns:
        tensor = column.reshape((dimension,) * boundary_leg_count)
        norm = float(np.linalg.norm(tensor))
        for generator in generators:
            total = np.zeros_like(tensor)
            for axis in range(boundary_leg_count):
                total += apply_generator_axis(tensor, axis, generator)
            max_boundary_gauss_residual = max(
                max_boundary_gauss_residual,
                float(np.linalg.norm(total)) / max(norm, 1e-14),
            )

    full_invariant_dimension = singlet_multiplicity(two_j, boundary_leg_count)
    # The normalized map is d^(N-1) C. Pushing forward the product-Haar
    # covariance I/d^N gives P_image/d^N.
    pushforward_covariance_eigenvalues = (
        dimension ** (vertex_count - 2)
    ) * singular_values**2
    positive = pushforward_covariance_eigenvalues[
        pushforward_covariance_eigenvalues > 1e-10
    ]
    covariance_trace = float(np.sum(positive))
    covariance_purity = float(np.sum(positive**2))
    haar_purity = 1.0 / full_invariant_dimension
    non_haar_distance = float(np.sqrt(
        1.0 / input_dimension - 1.0 / full_invariant_dimension
    ))
    two_copy_moment = product_vs_global_haar_two_copy_moment(
        two_j, vertex_count
    )

    return {
        "status": "exact pushforward of independent local Haar intertwiners",
        "ensemble": (
            f"{vertex_count} independent Haar-random unit vectors in the "
            "four-valent intertwiner space, fused along a chain of "
            f"{vertex_count - 1} invariant singlets"
        ),
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "local_intertwiner_dimension": int(local_basis.shape[1]),
        "product_input_dimension": input_dimension,
        "full_projective_image_dimension": full_projective_image_dimension,
        "product_haar_projective_support_dimension": product_haar_projective_support_dimension,
        "segre_codimension_in_image": segre_codimension_in_image,
        "vertex_count": vertex_count,
        "fused_boundary_leg_count": boundary_leg_count,
        "full_boundary_invariant_dimension": full_invariant_dimension,
        "fusion_map_rank": map_rank,
        "fusion_map_singular_values": singular_values.tolist(),
        "fusion_map_gram_max_abs_error_from_expected_scale": gram_error,
        "local_gauss_residual": float(local_gauss_residual),
        "boundary_gauss_max_relative_residual": max_boundary_gauss_residual,
        "pushforward_covariance_trace": covariance_trace,
        "pushforward_covariance_rank": int(positive.size),
        "pushforward_nonzero_covariance_eigenvalues": positive.tolist(),
        "pushforward_covariance_purity": covariance_purity,
        "full_boundary_haar_covariance_purity": haar_purity,
        "pushforward_vs_full_haar_covariance_frobenius_distance": non_haar_distance,
        **two_copy_moment,
        "passed": (
            local_basis.shape[1] == dimension
            and map_rank == input_dimension
            and full_invariant_dimension > input_dimension
            and segre_codimension_in_image > 0
            and gram_error < 2e-9
            and local_gauss_residual < 2e-9
            and max_boundary_gauss_residual < 2e-9
            and abs(covariance_trace - 1.0) < 2e-9
            and positive.size == input_dimension
            and np.max(np.abs(positive - 1.0 / input_dimension)) < 2e-9
            and abs(covariance_purity - 1.0 / input_dimension) < 2e-9
            and abs(haar_purity - 1.0 / full_invariant_dimension) < 2e-9
            and non_haar_distance > 1e-3
            and two_copy_moment["two_copy_moment_hilbert_schmidt_distance"] > 1e-3
        ),
        "interpretation": (
            "For this explicitly chosen independent local-Haar ensemble, "
            "A connected chain block maps the product input isometrically "
            "up to scale into a proper subspace of its boundary invariant "
            "space. Its exact pushforward covariance is therefore not the "
            "full boundary-Haar covariance. The exact two-copy moment also "
            "differs from Haar on the image subspace, so the pushforward is "
            "not Haar even after restricting to its support: it is concentrated "
            "on the lower-dimensional Segre variety of product states in the "
            "induced vertex-factorization. This demonstrates non-closure of this "
            "toy measure, not a statement about the unspecified "
            "physical BQG measure."
        ),
        "claim_boundary": (
            "Exact finite-dimensional result for the listed spin, chain "
            "block size, singlet edge pattern, and independent local Haar "
            "measures only."
        ),
    }


def analyze_block_family(
    specifications=((1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 2), (4, 2))
) -> dict[str, object]:
    if specifications is None:
        specifications = ((1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 2), (4, 2))
    controls = {
        f"j={two_j / 2.0};vertices={vertex_count}": analyze_independent_haar_fusion(
            two_j, vertex_count
        )
        for two_j, vertex_count in specifications
    }
    return {
        "status": "exact Haar pushforward under connected chain-block fusion",
        "tested_blocks": [
            {"spin": two_j / 2.0, "vertex_count": vertex_count}
            for two_j, vertex_count in specifications
        ],
        "passed": all(row["passed"] for row in controls.values()),
        "controls": controls,
        "claim_boundary": (
            "Exact pushforward of independent local Haar intertwiner measures "
            "for the listed finite spins and connected chain blocks only; not "
            "the physical BQG measure or a completed RG flow."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-j", type=int, action="append")
    parser.add_argument("--vertices", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = (
        analyze_independent_haar_fusion(args.two_j[0], args.vertices)
        if args.two_j and len(args.two_j) == 1
        else analyze_block_family(
            tuple((two_j, args.vertices) for two_j in args.two_j)
            if args.two_j else None
        )
    )
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
