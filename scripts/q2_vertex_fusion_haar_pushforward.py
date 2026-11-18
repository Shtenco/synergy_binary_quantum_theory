#!/usr/bin/env python3
"""Exact pushforward of independent Haar four-valent intertwiners under vertex fusion.

This is a deliberately specified toy ensemble, not a claim about the physical
BQG history measure.  It tests whether independent local Haar measures close
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


def total_generators(leg_count: int, generators: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    dimension = generators[0].shape[0]
    identity = np.eye(dimension, dtype=complex)
    totals = []
    for generator in generators:
        total = np.zeros((dimension**leg_count, dimension**leg_count), dtype=complex)
        for leg in range(leg_count):
            factors = [
                generator if axis == leg else identity
                for axis in range(leg_count)
            ]
            term = factors[0]
            for factor in factors[1:]:
                term = np.kron(term, factor)
            total += term
        totals.append(total)
    return tuple(totals)


def invariant_basis(leg_count: int, generators: tuple[np.ndarray, ...]) -> tuple[np.ndarray, float]:
    stacked = np.vstack(total_generators(leg_count, generators))
    _u, singular, vh = np.linalg.svd(stacked, full_matrices=True)
    threshold = (
        max(stacked.shape)
        * np.finfo(float).eps
        * (singular[0] if singular.size else 1.0)
        * 50.0
    )
    rank = int(np.sum(singular > threshold))
    basis = vh[rank:].conj().T
    return basis, float(np.linalg.norm(stacked @ basis))


def analyze_independent_haar_fusion(two_j: int = 1) -> dict[str, object]:
    if two_j != 1:
        raise ValueError("the exact two-vertex pushforward control currently targets two_j=1")
    dimension = two_j + 1
    generators = spin_generators(two_j)
    singlet = invariant_singlet_tensor(two_j)
    local_basis, local_gauss_residual = four_valent_intertwiner_basis(two_j)

    # The two local vertices are joined through leg 3 of vertex 0 and leg 0
    # of vertex 1.  The six remaining half-links form the fused block boundary.
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
    expected_gram = 0.25 * np.eye(gram.shape[0], dtype=complex)
    gram_error = float(np.max(np.abs(gram - expected_gram)))
    singular_values = np.linalg.svd(contraction_map, compute_uv=False)
    map_rank = int(np.linalg.matrix_rank(contraction_map, tol=1e-10))

    six_valent_basis, six_valent_gauss_residual = invariant_basis(6, generators)
    invariant_projector = six_valent_basis @ six_valent_basis.conj().T
    image_projector = contraction_map @ np.linalg.solve(
        gram, contraction_map.conj().T
    )
    invariance_error = float(
        np.linalg.norm((np.eye(dimension**6) - invariant_projector) @ contraction_map)
    )
    image_invariant_overlap_error = float(
        np.linalg.norm(image_projector - invariant_projector @ image_projector @ invariant_projector)
    )

    # A normalized local Haar vector is mapped by 2*M to a normalized fused
    # vector because M^dagger M = I/4.  Independent local Haar second moments
    # are I/2 tensor I/2 = I/4, hence the exact pushforward covariance is M M^dagger.
    pushforward_covariance = contraction_map @ contraction_map.conj().T
    coarse_haar_covariance = invariant_projector / six_valent_basis.shape[1]
    covariance_trace = float(np.trace(pushforward_covariance).real)
    covariance_purity = float(np.trace(pushforward_covariance @ pushforward_covariance).real)
    haar_purity = float(np.trace(coarse_haar_covariance @ coarse_haar_covariance).real)
    non_haar_distance = float(
        np.linalg.norm(pushforward_covariance - coarse_haar_covariance)
    )
    covariance_eigenvalues = np.linalg.eigvalsh(pushforward_covariance)
    positive = covariance_eigenvalues[covariance_eigenvalues > 1e-10]

    return {
        "status": "exact pushforward of independent local Haar intertwiners",
        "ensemble": (
            "two independent Haar-random unit vectors in the four-valent "
            "j=1/2 intertwiner space, contracted through one invariant singlet"
        ),
        "spin": 0.5,
        "local_intertwiner_dimension": int(local_basis.shape[1]),
        "product_input_dimension": int(local_basis.shape[1] ** 2),
        "fused_boundary_leg_count": 6,
        "full_six_valent_invariant_dimension": int(six_valent_basis.shape[1]),
        "fusion_map_rank": map_rank,
        "fusion_map_singular_values": singular_values.tolist(),
        "fusion_map_gram_max_abs_error_from_I_over_4": gram_error,
        "local_gauss_residual": local_gauss_residual,
        "six_valent_gauss_residual": six_valent_gauss_residual,
        "fusion_image_invariant_subspace_residual": invariance_error,
        "image_projector_invariance_residual": image_invariant_overlap_error,
        "pushforward_covariance_trace": covariance_trace,
        "pushforward_covariance_rank": int(positive.size),
        "pushforward_nonzero_covariance_eigenvalues": positive.tolist(),
        "pushforward_covariance_purity": covariance_purity,
        "full_six_valent_haar_covariance_purity": haar_purity,
        "pushforward_vs_full_haar_covariance_frobenius_distance": non_haar_distance,
        "passed": (
            map_rank == 4
            and six_valent_basis.shape[1] == 5
            and gram_error < 2e-10
            and local_gauss_residual < 2e-10
            and six_valent_gauss_residual < 2e-10
            and invariance_error < 2e-10
            and image_invariant_overlap_error < 2e-10
            and abs(covariance_trace - 1.0) < 2e-10
            and positive.size == 4
            and abs(covariance_purity - 0.25) < 2e-10
            and abs(haar_purity - 0.2) < 2e-10
            and non_haar_distance > 0.1
        ),
        "interpretation": (
            "For this explicitly chosen independent local-Haar ensemble, "
            "two-vertex fusion maps the four-dimensional product input "
            "isometrically (up to scale) into a four-dimensional subspace of "
            "the five-dimensional six-valent invariant space.  Its exact "
            "pushforward covariance is therefore not the full six-valent "
            "Haar covariance.  This demonstrates non-closure of this toy "
            "measure under fusion, not a statement about the unspecified "
            "physical BQG measure."
        ),
        "claim_boundary": (
            "Exact finite-dimensional result for j=1/2, two vertices, one "
            "singlet edge, and independent local Haar measures only."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze_independent_haar_fusion()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
