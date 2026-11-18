#!/usr/bin/env python3
"""Four-valent SU(2) star: full incident-link QFI and gauge-null audit."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import spatial_qubit_geometry_gate as SQ

PAULI = SQ.PAULI


def leg_generator(pauli: np.ndarray, leg: int) -> np.ndarray:
    return SQ.kron_all([
        0.5 * pauli if index == leg else SQ.I2
        for index in range(4)
    ])


def pair_casimir_matrix(state: np.ndarray) -> np.ndarray:
    """C_ij=(4/3)<J_i dot J_j>; equals the 3x3-axis QFI block scalar."""
    result = np.zeros((4, 4), dtype=float)
    for i, j in itertools.product(range(4), repeat=2):
        pair = sum(
            (leg_generator(PAULI[axis], i) @ leg_generator(PAULI[axis], j)
             for axis in range(3)),
            np.zeros((16, 16), dtype=complex),
        )
        result[i, j] = (4.0 / 3.0) * float(np.vdot(state, pair @ state).real)
    return result


def random_intertwiner(rng: np.random.Generator) -> np.ndarray:
    coefficients = rng.normal(size=2) + 1j * rng.normal(size=2)
    coefficients /= np.linalg.norm(coefficients)
    return SQ.logical_basis() @ coefficients


def run(seed: int = 20261118, trials: int = 32) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    row_errors = []
    diagonal_errors = []
    symmetry_errors = []
    minimum_eigenvalues = []
    relative_minimum_eigenvalues = []
    gauge_residuals = []
    sample_matrices = []

    for _ in range(trials):
        state = random_intertwiner(rng)
        metric = pair_casimir_matrix(state)
        sample_matrices.append(metric)
        row_errors.append(float(np.linalg.norm(metric @ np.ones(4))))
        diagonal_errors.append(float(np.max(np.abs(np.diag(metric) - 1.0))))
        symmetry_errors.append(float(np.linalg.norm(metric - metric.T)))
        minimum_eigenvalues.append(float(np.linalg.eigvalsh(metric)[0]))
        # Gauge-fix the common SU(2) twist by setting leg 3's parameter to zero.
        relative = metric[:3, :3]
        relative_minimum_eigenvalues.append(float(np.linalg.eigvalsh(relative)[0]))
        # The full 12x12 QFI is C tensor I_3 in edge-major parameter ordering.
        full_qfi = np.kron(metric, np.eye(3))
        gauge_vectors = np.tile(np.eye(3), (4, 1))
        gauge_residuals.append(float(np.linalg.norm(full_qfi @ gauge_vectors)))

    mean_metric = np.mean(sample_matrices, axis=0)
    checks = {
        "each_incident_link_has_unit_diagonal_qfi": max(diagonal_errors) < 1e-11,
        "pair_metric_is_symmetric": max(symmetry_errors) < 1e-11,
        "Gauss_law_produces_common_twist_null_mode": max(row_errors) < 1e-11,
        "full_12_parameter_qfi_is_psd": min(minimum_eigenvalues) > -1e-10,
        "gauge_fixed_three_parameter_metric_is_positive": min(
            relative_minimum_eigenvalues
        ) > 1e-10,
        "full_qfi_annihilates_three_global_su2_gauge_directions": max(
            gauge_residuals
        ) < 1e-10,
    }
    return {
        "status": "four-valent star relational QFI and Gauss-null control",
        "passed": all(checks.values()),
        "seed": seed,
        "trials": trials,
        "local_parameter_count": 12,
        "edge_axis_qfi_block": "F_(i,a),(j,b) = C_ij delta_ab",
        "pair_metric_definition": "C_ij=(4/3)<iota|J_i dot J_j|iota>",
        "mean_pair_metric": mean_metric.tolist(),
        "mean_pair_metric_eigenvalues": np.linalg.eigvalsh(mean_metric).tolist(),
        "max_common_twist_null_residual": max(row_errors),
        "max_full_qfi_gauge_null_residual": max(gauge_residuals),
        "minimum_full_pair_metric_eigenvalue": min(minimum_eigenvalues),
        "minimum_gauge_fixed_relative_metric_eigenvalue": min(relative_minimum_eigenvalues),
        "gauge_fix": "set the fourth incident-link twist parameter to zero",
        "gauge_fixed_relative_metric_mean": np.mean(
            np.asarray(sample_matrices)[:, :3, :3], axis=0
        ).tolist(),
        "checks": checks,
        "interpretation": (
            "The four incident link twists do not define four independent "
            "positive edge weights. Their 12-parameter QFI has three exact "
            "null directions from the common SU(2) Gauss rotation. The "
            "remaining edge variables are coupled by the full 4x4 pair "
            "Casimir matrix; after fixing one link as reference, the induced "
            "3x3 relative metric is positive for the sampled generic states."
        ),
        "claim_boundary": (
            "This is the local four-valent invariant-node twist metric. "
            "It does not include a full multi-node mixed-state reduction, "
            "higher-spin sums, the BQG history-amplitude measure, or a "
            "continuum limit."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=32)
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
