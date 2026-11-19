#!/usr/bin/env python3
"""Global SU(2)^V quotient and two-sided endpoint QFI for finite spin networks.

The base coordinates are all 4V half-link rotations. Vertex Gauss rotations
span a 3V-dimensional subspace. The quotient is represented by its orthogonal
complement; reference-leg sections are checked as alternative coordinates.
Internal edges retain independent rotations at both endpoints.
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
    contract_graph,
)


def vertex_gauge_generators(node_count: int) -> np.ndarray:
    """Return G (12V x 3V), one common SU(2) rotation per vertex."""
    if node_count < 1:
        raise ValueError("node_count must be positive")
    generators = np.zeros((12 * node_count, 3 * node_count), dtype=float)
    for node in range(node_count):
        for axis in range(3):
            for leg in range(4):
                generators[12 * node + 3 * leg + axis, 3 * node + axis] = 1.0
    return generators


def two_sided_edge_parameter_map(
    node_count: int, edges: list[tuple[int, int, int, int]]
) -> np.ndarray:
    """Map 6E endpoint parameters into 12V half-link coordinates.

    Ordering: (edge, endpoint u/v, axis x/y/z). Both endpoint frames are
    independent; no relative-frame convention is silently imposed here.
    """
    if node_count < 1:
        raise ValueError("node_count must be positive")
    mapping = np.zeros((12 * node_count, 6 * len(edges)), dtype=float)
    occupied: set[tuple[int, int]] = set()
    for edge_index, edge in enumerate(edges):
        if len(edge) != 4:
            raise ValueError("each edge must be (u, leg_u, v, leg_v)")
        u, leg_u, v, leg_v = map(int, edge)
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("edge endpoints must be distinct valid nodes")
        if not (0 <= leg_u < 4 and 0 <= leg_v < 4):
            raise ValueError("half-link slots must be in [0,3]")
        for endpoint in ((u, leg_u), (v, leg_v)):
            if endpoint in occupied:
                raise ValueError("a half-link cannot belong to more than one edge")
            occupied.add(endpoint)
        for endpoint, (node, leg) in enumerate(((u, leg_u), (v, leg_v))):
            for axis in range(3):
                column = 6 * edge_index + 3 * endpoint + axis
                mapping[12 * node + 3 * leg + axis, column] = 1.0
    return mapping


def _section_map(node_count: int, reference_leg: int) -> np.ndarray:
    """Gauge section fixing one chosen half-link per vertex to zero."""
    if not 0 <= reference_leg < 4:
        raise ValueError("reference_leg must be in [0,3]")
    columns = []
    for node in range(node_count):
        for leg in range(4):
            if leg == reference_leg:
                continue
            for axis in range(3):
                column = np.zeros(12 * node_count, dtype=float)
                column[12 * node + 3 * leg + axis] = 1.0
                columns.append(column)
    return np.column_stack(columns)


def quotient_diagnostics(
    qfi: np.ndarray, node_count: int, reference_legs: tuple[int, ...] = (0, 3)
) -> dict[str, object]:
    """Construct the orthogonal quotient and compare gauge sections.

    Section spectra are generalized spectra relative to the quotient Euclidean
    metric. Raw principal-submatrix eigenvalues are coordinate-dependent and
    are deliberately not compared as physical invariants.
    """
    expected = 12 * node_count
    matrix = np.asarray(qfi, dtype=float)
    if matrix.shape != (expected, expected):
        raise ValueError("qfi shape must be (12V,12V)")
    if not np.allclose(matrix, matrix.T, atol=1e-10):
        raise ValueError("qfi must be symmetric")
    gauge = vertex_gauge_generators(node_count)
    gauge_residual = float(np.linalg.norm(matrix @ gauge))
    # Columns of gauge are mutually orthogonal, each with norm 2.
    gauge_basis = gauge / 2.0
    complete, _ = np.linalg.qr(gauge_basis, mode="complete")
    quotient_basis = complete[:, 3 * node_count:]
    quotient_metric = quotient_basis.T @ matrix @ quotient_basis
    quotient_metric = 0.5 * (quotient_metric + quotient_metric.T)
    quotient_eigenvalues = np.linalg.eigvalsh(quotient_metric)
    scale = max(1.0, float(np.max(np.abs(quotient_eigenvalues))))
    threshold = 1e-9 * scale
    section_results = {}
    section_spectra = []
    for reference_leg in reference_legs:
        section = _section_map(node_count, reference_leg)
        coordinates = quotient_basis.T @ section
        singular_values = np.linalg.svd(coordinates, compute_uv=False)
        if singular_values[-1] <= 1e-10:
            raise AssertionError("reference-leg section does not span the quotient")
        coordinate_metric = coordinates.T @ quotient_metric @ coordinates
        coordinate_metric = 0.5 * (coordinate_metric + coordinate_metric.T)
        gram = coordinates.T @ coordinates
        gram_values, gram_vectors = np.linalg.eigh(gram)
        if gram_values[0] <= 1e-12:
            raise AssertionError("section quotient Gram matrix is singular")
        gram_inv_sqrt = (
            gram_vectors * (1.0 / np.sqrt(gram_values))
        ) @ gram_vectors.T
        whitened = gram_inv_sqrt @ coordinate_metric @ gram_inv_sqrt
        section_eigenvalues = np.linalg.eigvalsh(0.5 * (whitened + whitened.T))
        section_spectra.append(section_eigenvalues)
        section_results[str(reference_leg)] = {
            "coordinate_count": int(section.shape[1]),
            "coordinate_condition_number": float(singular_values[0] / singular_values[-1]),
            "generalized_minimum_eigenvalue": float(section_eigenvalues[0]),
            "generalized_maximum_eigenvalue": float(section_eigenvalues[-1]),
            "spectrum_matches_orthogonal_quotient": bool(np.allclose(
                section_eigenvalues, quotient_eigenvalues, atol=2e-8, rtol=2e-8
            )),
        }
    max_section_spectrum_error = max(
        float(np.max(np.abs(values - quotient_eigenvalues)))
        for values in section_spectra
    )
    return {
        "full_parameter_count": expected,
        "gauge_generator_count": 3 * node_count,
        "quotient_dimension": 9 * node_count,
        "gauge_annihilation_residual_frobenius": gauge_residual,
        "quotient_minimum_eigenvalue": float(quotient_eigenvalues[0]),
        "quotient_maximum_eigenvalue": float(quotient_eigenvalues[-1]),
        "quotient_numerical_rank": int(np.sum(quotient_eigenvalues > threshold)),
        "quotient_nullity_beyond_gauge": int(
            np.sum(np.abs(quotient_eigenvalues) <= threshold)
        ),
        "quotient_is_psd_up_to_roundoff": bool(quotient_eigenvalues[0] >= -1e-9 * scale),
        # Preserve the complete correlated quotient matrix as the primary
        # finite-network result; spectra alone would discard off-diagonal data.
        "quotient_qfi_matrix": quotient_metric.tolist(),
        "quotient_qfi_offdiagonal_frobenius_norm": float(
            np.linalg.norm(quotient_metric - np.diag(np.diag(quotient_metric)))
        ),
        "reference_sections": section_results,
        "max_reference_section_spectrum_error": max_section_spectrum_error,
        "reference_section_spectra_agree": max_section_spectrum_error < 2e-8,
    }



def projective_tangent_jacobian(
    node_states: list[np.ndarray],
    edges: list[tuple[int, int, int, int]],
) -> tuple[np.ndarray, float]:
    """Real-parameter tangent map of the normalized boundary ray.

    The returned complex matrix has one column per half-link/axis parameter.
    Each column is projected orthogonally to the state, removing both radial
    normalization and physically irrelevant global phase.
    """
    raw, boundary = contract_graph(node_states, edges)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("graph singlet contraction has zero support")
    state = raw / np.sqrt(probability)
    columns = []
    from q2_relational_graph_qfi_gate import PAULI
    for node in range(len(node_states)):
        for leg in range(4):
            for pauli in PAULI:
                derivative_raw, derivative_boundary = contract_graph(
                    node_states, edges,
                    half_link_twist=(node, leg, pauli),
                )
                if derivative_boundary != boundary:
                    raise AssertionError("boundary ordering changed under tangent variation")
                derivative = derivative_raw / np.sqrt(probability)
                # Project out the complete complex state direction, including phase.
                derivative = derivative - state * np.vdot(state, derivative)
                columns.append(derivative)
    return np.column_stack(columns), probability


def stabilizer_diagnostics(
    tangent: np.ndarray, node_count: int, qfi: np.ndarray
) -> dict[str, object]:
    """Identify extra null directions as infinitesimal stabilizers of the boundary ray."""
    if tangent.shape[1] != 12 * node_count:
        raise ValueError("tangent column count must equal 12V")
    gauge = vertex_gauge_generators(node_count)
    gauge_basis = gauge / 2.0
    complete, _ = np.linalg.qr(gauge_basis, mode="complete")
    quotient_basis = complete[:, 3 * node_count:]
    quotient_tangent = tangent @ quotient_basis
    realified = np.concatenate(
        [quotient_tangent.real, quotient_tangent.imag], axis=0
    )
    _u, singular_values, vh = np.linalg.svd(realified, full_matrices=False)
    largest = float(singular_values[0]) if singular_values.size else 0.0
    qfi_quotient = quotient_basis.T @ qfi @ quotient_basis
    qfi_quotient = 0.5 * (qfi_quotient + qfi_quotient.T)
    qfi_eigenvalues = np.linalg.eigvalsh(qfi_quotient)
    qfi_scale = max(1.0, float(np.max(np.abs(qfi_eigenvalues))))
    qfi_tolerance = 1e-9 * qfi_scale
    # Since F = 4 R(D)^T R(D), use the square-root-consistent SVD cutoff.
    tolerance = float(np.sqrt(qfi_tolerance / 4.0))
    rank = int(np.sum(singular_values > tolerance))
    quotient_dimension = 9 * node_count
    nullity = quotient_dimension - rank
    qfi_rank = int(np.sum(qfi_eigenvalues > qfi_tolerance))
    participation = []
    if nullity:
        null_coordinates = vh[-nullity:, :]
        null_full = quotient_basis @ null_coordinates.T
        for node in range(node_count):
            leg_weights = []
            for leg in range(4):
                indices = [12 * node + 3 * leg + axis for axis in range(3)]
                leg_weights.append(float(np.sum(null_full[indices, :] ** 2)))
            participation.append({
                "node": node,
                "half_link_weights": leg_weights,
                "total_weight": float(sum(leg_weights)),
            })
    return {
        "projective_tangent_shape": list(tangent.shape),
        "quotient_tangent_rank": rank,
        "quotient_qfi_rank": qfi_rank,
        "quotient_tangent_nullity": nullity,
        "tangent_qfi_rank_agrees": rank == qfi_rank,
        "smallest_tangent_singular_value": (
            float(singular_values[-1]) if singular_values.size else 0.0
        ),
        "largest_tangent_singular_value": largest,
        "null_mode_half_link_participation": participation,
        "interpretation": (
            "Additional null vectors are infinitesimal stabilizers of the "
            "normalized boundary-state ray in the chosen finite contraction. "
            "Their physical status beyond this model is not established."
        ),
    }

def run(seed: int = 20261118, trials: int = 4) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    results = {}
    for name, (node_count, edges) in graph_fixtures().items():
        gauge_rows = []
        endpoint_rows = []
        for _ in range(trials):
            states = [random_intertwiner(rng) for _ in range(node_count)]
            _probability, full_qfi, symmetry, _boundary = full_half_link_qfi(states, edges)
            gauge = quotient_diagnostics(full_qfi, node_count)
            tangent, _ = projective_tangent_jacobian(states, edges)
            stabilizers = stabilizer_diagnostics(tangent, node_count, full_qfi)
            gauge["stabilizer_diagnostics"] = stabilizers
            mapping = two_sided_edge_parameter_map(node_count, edges)
            endpoint_qfi = mapping.T @ full_qfi @ mapping
            endpoint_qfi = 0.5 * (endpoint_qfi + endpoint_qfi.T)
            endpoint_eigenvalues = np.linalg.eigvalsh(endpoint_qfi)
            endpoint_rows.append({
                "parameter_count": int(mapping.shape[1]),
                "qfi_shape": list(endpoint_qfi.shape),
                "qfi_matrix": endpoint_qfi.tolist(),
                "minimum_eigenvalue": float(endpoint_eigenvalues[0]),
                "maximum_eigenvalue": float(endpoint_eigenvalues[-1]),
                "symmetry_error": float(np.linalg.norm(endpoint_qfi - endpoint_qfi.T)),
                "both_endpoints_independent": bool(mapping.shape[1] == 6 * len(edges)),
            })
            gauge["full_half_link_symmetry_error"] = symmetry
            gauge_rows.append(gauge)
        results[name] = {
            "nodes": node_count,
            "internal_edges": len(edges),
            "two_sided_endpoint_parameter_count": 6 * len(edges),
            "two_sided_endpoint_qfi_shape": [6 * len(edges), 6 * len(edges)],
            "max_gauge_annihilation_residual": max(
                row["gauge_annihilation_residual_frobenius"] for row in gauge_rows
            ),
            "quotient_dimension": gauge_rows[0]["quotient_dimension"],
            "quotient_numerical_rank_min": min(
                row["quotient_numerical_rank"] for row in gauge_rows
            ),
            "quotient_nullity_beyond_gauge_max": max(
                row["quotient_nullity_beyond_gauge"] for row in gauge_rows
            ),
            "quotient_minimum_eigenvalue": min(
                row["quotient_minimum_eigenvalue"] for row in gauge_rows
            ),
            "max_reference_section_spectrum_error": max(
                row["max_reference_section_spectrum_error"] for row in gauge_rows
            ),
            "reference_section_spectra_agree": all(
                row["reference_section_spectra_agree"] for row in gauge_rows
            ),
            "reference_sections": gauge_rows[0]["reference_sections"],
            "quotient_qfi_matrices_by_trial": [
                row["quotient_qfi_matrix"] for row in gauge_rows
            ],
            "quotient_qfi_offdiagonal_frobenius_norm_max": max(
                row["quotient_qfi_offdiagonal_frobenius_norm"] for row in gauge_rows
            ),
            "two_sided_endpoint_qfi_matrices_by_trial": [
                row["qfi_matrix"] for row in endpoint_rows
            ],
            "quotient_tangent_rank_min": min(
                row["stabilizer_diagnostics"]["quotient_tangent_rank"] for row in gauge_rows
            ),
            "projective_stabilizer_nullity_max": max(
                row["stabilizer_diagnostics"]["quotient_tangent_nullity"] for row in gauge_rows
            ),
            "stabilizer_half_link_participation_sample": gauge_rows[0]["stabilizer_diagnostics"]["null_mode_half_link_participation"],
            "stabilizer_rank_consistency": all(
                row["stabilizer_diagnostics"]["tangent_qfi_rank_agrees"] for row in gauge_rows
            ),
            "endpoint_qfi_minimum_eigenvalue": min(
                row["minimum_eigenvalue"] for row in endpoint_rows
            ),
            "endpoint_qfi_max_symmetry_error": max(
                row["symmetry_error"] for row in endpoint_rows
            ),
            "checks": {
                "global_su2v_gauge_annihilated": max(
                    row["gauge_annihilation_residual_frobenius"] for row in gauge_rows
                ) < 1e-9,
                "quotient_psd": all(
                    row["quotient_is_psd_up_to_roundoff"] for row in gauge_rows
                ),
                "reference_section_independence": all(
                    row["reference_section_spectra_agree"] for row in gauge_rows
                ),
                "two_sided_endpoint_parameters_present": all(
                    row["both_endpoints_independent"] for row in endpoint_rows
                ),
                "two_sided_endpoint_qfi_psd": min(
                    row["minimum_eigenvalue"] for row in endpoint_rows
                ) > -1e-9,
                "extra_nulls_are_projective_state_stabilizers": all(
                    row["stabilizer_diagnostics"]["tangent_qfi_rank_agrees"]
                    and row["stabilizer_diagnostics"]["quotient_tangent_nullity"]
                    == row["quotient_nullity_beyond_gauge"]
                    for row in gauge_rows
                ),
            },
        }
        results[name]["passed"] = all(results[name]["checks"].values())
    return {
        "status": "finite spin-network global SU(2)^V quotient control",
        "passed": all(row["passed"] for row in results.values()),
        "seed": seed,
        "trials_per_graph": trials,
        "base_coordinates": "all four half-link SU(2) infinitesimal rotations at every vertex",
        "gauge_action": "one common SU(2) rotation on all four half-links at each invariant vertex",
        "quotient_method": "orthogonal complement of the 3V-dimensional Gauss gauge subspace",
        "section_comparison": "reference half-link 0 versus reference half-link 3; compare generalized spectra using the induced quotient Gram metric",
        "edge_parameterization": "six independent endpoint parameters per internal edge, three at each endpoint",
        "claim_boundary": (
            "This is a finite spin-1/2 path and cycle control. The quotient may "
            "retain additional physical null modes; no arbitrary-graph theorem, "
            "holonomy dynamics, BQG history measure, causal-PL pushforward, "
            "or continuum spectral-dimension claim follows."
        ),
        "graphs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=4)
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
