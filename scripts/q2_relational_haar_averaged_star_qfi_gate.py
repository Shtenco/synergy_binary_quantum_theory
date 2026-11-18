#!/usr/bin/env python3
"""Exact Haar-ensemble average of four-valent spin-j pure-state QFI.

The average is over normalized pure intertwiners with Haar-distributed
coefficients in Inv(V_j^tensor4). It is not the QFI of the Haar-averaged
mixed density operator. Optional spin weights define a classical ensemble
average of sector QFIs; weights are inputs, not derived BQG amplitudes.
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

from q2_relational_higher_spin_network_qfi_gate import (
    four_valent_intertwiner_basis,
    node_generator,
    spin_generators,
)


def haar_averaged_star_qfi(two_j: int) -> dict[str, object]:
    if not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer")
    basis, gauss_residual = four_valent_intertwiner_basis(two_j)
    generators = spin_generators(two_j)
    dimension = two_j + 1
    intertwiner_dimension = basis.shape[1]
    qfi = np.zeros((12, 12), dtype=float)
    for left in range(4):
        for right in range(4):
            for axis, generator_a in enumerate(generators):
                op_left = node_generator(generator_a, left, dimension)
                for other_axis, generator_b in enumerate(generators):
                    op_right = node_generator(generator_b, right, dimension)
                    reduced = basis.conj().T @ (op_left @ op_right) @ basis
                    expectation = np.trace(reduced) / intertwiner_dimension
                    qfi[3 * left + axis, 3 * right + other_axis] = 4.0 * float(expectation.real)
    qfi = 0.5 * (qfi + qfi.T)
    j = two_j / 2.0
    q = 4.0 * j * (j + 1.0) / 3.0
    leg_metric = (4.0 / 3.0) * np.eye(4) - (1.0 / 3.0) * np.ones((4, 4))
    expected = q * np.kron(leg_metric, np.eye(3))
    eigenvalues = np.linalg.eigvalsh(qfi)
    gauss_residuals = []
    for axis in range(3):
        direction = np.zeros(12)
        direction[axis::3] = 1.0
        gauss_residuals.append(float(np.linalg.norm(qfi @ direction)))
    return {
        "two_j": two_j,
        "spin": j,
        "intertwiner_dimension": int(intertwiner_dimension),
        "gauss_nullspace_residual": gauss_residual,
        "expected_single_leg_qfi": q,
        "max_exact_matrix_formula_error": float(np.linalg.norm(qfi - expected)),
        "mean_pair_dot_product": -j * (j + 1.0) / 3.0,
        "mean_cross_link_qfi": -q / 3.0,
        "eigenvalues": eigenvalues.tolist(),
        "zero_eigenvalue_count": int(np.sum(np.abs(eigenvalues) < 1e-9)),
        "positive_eigenvalue_value": 4.0 * q / 3.0,
        "positive_eigenvalue_count": int(np.sum(eigenvalues > 1e-9)),
        "max_gauss_null_residual": max(gauss_residuals),
        "passed": bool(
            intertwiner_dimension == two_j + 1
            and gauss_residual < 2e-10
            and np.linalg.norm(qfi - expected) < 2e-9
            and max(gauss_residuals) < 2e-9
            and eigenvalues[0] > -2e-9
            and int(np.sum(np.abs(eigenvalues) < 1e-9)) == 3
            and int(np.sum(eigenvalues > 1e-9)) == 9
        ),
        "qfi": qfi.tolist(),
    }


def weighted_spin_ensemble(
    two_j_values: tuple[int, ...], weights: tuple[float, ...]
) -> dict[str, object]:
    if not two_j_values or len(two_j_values) != len(weights):
        raise ValueError("spin labels and weights must be nonempty and have equal length")
    if len(set(two_j_values)) != len(two_j_values):
        raise ValueError("spin labels must be unique")
    if any(not isinstance(value, int) or value < 1 for value in two_j_values):
        raise ValueError("each two_j must be a positive integer")
    weight_array = np.asarray(weights, dtype=float)
    if not np.all(np.isfinite(weight_array)) or np.any(weight_array < 0.0):
        raise ValueError("weights must be finite and nonnegative")
    total = float(weight_array.sum())
    if total <= 0.0:
        raise ValueError("at least one spin weight must be positive")
    probabilities = weight_array / total
    sectors = [haar_averaged_star_qfi(two_j) for two_j in two_j_values]
    combined = sum(
        (probabilities[index] * np.asarray(sectors[index]["qfi"], dtype=float)
         for index in range(len(sectors))),
        np.zeros((12, 12), dtype=float),
    )
    q_bar = sum(
        probabilities[index] * float(sectors[index]["expected_single_leg_qfi"])
        for index in range(len(sectors))
    )
    j_metric = q_bar * np.kron(
        (4.0 / 3.0) * np.eye(4) - (1.0 / 3.0) * np.ones((4, 4)),
        np.eye(3),
    )
    eigenvalues = np.linalg.eigvalsh(combined)
    return {
        "spin_labels": [value / 2.0 for value in two_j_values],
        "normalized_weights": probabilities.tolist(),
        "mean_single_leg_qfi": q_bar,
        "mean_cross_link_qfi": -q_bar / 3.0,
        "max_linearity_formula_error": float(np.linalg.norm(combined - j_metric)),
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "zero_eigenvalue_count": int(np.sum(np.abs(eigenvalues) < 1e-9)),
        "positive_eigenvalue_count": int(np.sum(eigenvalues > 1e-9)),
        "passed": bool(
            np.linalg.norm(combined - j_metric) < 2e-9
            and eigenvalues[0] > -2e-9
            and int(np.sum(np.abs(eigenvalues) < 1e-9)) == 3
            and int(np.sum(eigenvalues > 1e-9)) == 9
        ),
        "qfi": combined.tolist(),
    }


def run() -> dict[str, object]:
    spins = [haar_averaged_star_qfi(two_j) for two_j in (1, 2, 3, 4)]
    # Illustrative distribution only. It is deliberately not presented as
    # a BQG-derived spin amplitude or as a fitted physical prior.
    ensemble = weighted_spin_ensemble((1, 2, 3, 4), (0.1, 0.2, 0.3, 0.4))
    return {
        "status": "exact Haar-ensemble averaged four-valent spin QFI",
        "passed": all(item["passed"] for item in spins) and ensemble["passed"],
        "average_type": (
            "classical average of pure-intertwiner QFIs with Haar-distributed "
            "intertwiner labels; not QFI of the averaged mixed state"
        ),
        "exact_formula": {
            "single_leg": "q_j=4*j*(j+1)/3",
            "mean_pair_dot": "-j*(j+1)/3",
            "mean_cross_link": "-q_j/3",
            "leg_metric": "q_j*(4/3*I_4 - 1/3*ones_4)",
            "spectrum": "0 (multiplicity 3), 4*q_j/3 (multiplicity 9)",
        },
        "spin_sector_results": spins,
        "illustrative_spin_weights": ensemble,
        "claim_boundary": (
            "The Haar measure is a mathematical ensemble over pure intertwiners. "
            "The example spin weights are illustrative inputs, not derived from "
            "BQG amplitudes. A physical spin sum needs a specified history measure, "
            "superselection/label-access assumptions, and normalization."
        ),
    }


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2, default=_json_default)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
