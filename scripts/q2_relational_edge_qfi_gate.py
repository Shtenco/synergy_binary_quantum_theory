#!/usr/bin/env python3
"""Relational QFI control for two Gauss-invariant four-valent qubit blocks.

Two normalized four-qubit SU(2) intertwiners are glued across one shared
spin-1/2 face by the invariant epsilon tensor. A relative SU(2) holonomy is
inserted on one half-link before contraction. The six uncontracted boundary
qubits are retained as the output state.
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


def node_generator(pauli: np.ndarray, site: int) -> np.ndarray:
    return SQ.kron_all([
        0.5 * pauli if index == site else SQ.I2
        for index in range(4)
    ])


def contract_shared_face(two_block_state: np.ndarray) -> np.ndarray:
    """Contract A's fourth and B's first qubit with the normalized singlet."""
    tensor = np.asarray(two_block_state, dtype=complex).reshape((2,) * 8)
    return np.tensordot(
        EPSILON, tensor, axes=([0, 1], [3, 4])
    ).reshape(-1)


def relational_qfi_matrix(
    intertwiner_a: np.ndarray,
    intertwiner_b: np.ndarray,
) -> tuple[float, np.ndarray]:
    """QFI of the normalized six-qubit boundary state under B-half-link twists."""
    psi = np.kron(intertwiner_a, intertwiner_b)
    glued = contract_shared_face(psi)
    probability = float(np.vdot(glued, glued).real)
    if probability <= 1e-14:
        raise ValueError("shared-face singlet contraction has zero support")

    state = glued / math.sqrt(probability)
    derivatives = []
    for pauli in PAULI:
        generator = node_generator(pauli, 0)
        derivative_raw = contract_shared_face(
            np.kron(intertwiner_a, -1j * (generator @ intertwiner_b))
        )
        # Differentiate normalization as well as the unnormalized contraction.
        radial = float(np.vdot(glued, derivative_raw).real) / probability
        derivatives.append(derivative_raw / math.sqrt(probability) - radial * state)

    qfi = np.zeros((3, 3), dtype=float)
    overlaps = [np.vdot(state, derivative) for derivative in derivatives]
    for i, j in itertools.product(range(3), repeat=2):
        qfi[i, j] = 4.0 * float(np.real(
            np.vdot(derivatives[i], derivatives[j])
            - np.conjugate(overlaps[i]) * overlaps[j]
        ))
    return probability, qfi


def random_intertwiner(rng: np.random.Generator) -> np.ndarray:
    coefficients = rng.normal(size=2) + 1j * rng.normal(size=2)
    coefficients /= np.linalg.norm(coefficients)
    return SQ.logical_basis() @ coefficients


def run(seed: int = 20261118, trials: int = 16) -> dict[str, object]:
    if trials < 1:
        raise ValueError("trials must be positive")
    rng = np.random.default_rng(seed)
    probabilities = []
    matrices = []
    gauss_residuals = []
    for _ in range(trials):
        intertwiner_a = random_intertwiner(rng)
        intertwiner_b = random_intertwiner(rng)
        probability, qfi = relational_qfi_matrix(intertwiner_a, intertwiner_b)
        probabilities.append(probability)
        matrices.append(qfi)

        for node_state in (intertwiner_a, intertwiner_b):
            for axis in range(3):
                total = sum(
                    (0.5 * SQ.kron_all([
                        SQ.PAULI[axis] if site == target else SQ.I2
                        for site in range(4)
                    ]) for target in range(4)),
                    np.zeros((16, 16), dtype=complex),
                )
                gauss_residuals.append(float(np.linalg.norm(total @ node_state)))

    mean_qfi = np.mean(np.asarray(matrices), axis=0)
    max_qfi_error = float(max(np.linalg.norm(matrix - np.eye(3)) for matrix in matrices))
    checks = {
        "singlet_gluing_probability_is_one_quarter": max(
            abs(value - 0.25) for value in probabilities
        ) < 1e-12,
        "each_node_satisfies_gauss_constraint": max(gauss_residuals) < 1e-12,
        "relational_qfi_is_identity_for_every_sample": max_qfi_error < 1e-11,
        "mean_relational_qfi_is_identity": float(
            np.linalg.norm(mean_qfi - np.eye(3))
        ) < 1e-12,
    }
    return {
        "status": "two-block spin-1/2 singlet-gluing relational QFI control",
        "passed": all(checks.values()),
        "seed": seed,
        "trials": trials,
        "shared_face": "spin-1/2 half-links contracted with normalized SU(2) epsilon singlet",
        "boundary_state": "normalized pure state on the six uncontracted qubits",
        "relative_generator": "J^a on the B-side shared half-link before invariant contraction",
        "gluing_probability_min": min(probabilities),
        "gluing_probability_max": max(probabilities),
        "relational_qfi_mean": mean_qfi.tolist(),
        "relational_qfi_max_frobenius_error_from_identity": max_qfi_error,
        "normalized_qfi_metric_F_over_4": (mean_qfi / 4.0).tolist(),
        "max_node_gauss_residual": max(gauss_residuals),
        "checks": checks,
        "claim_boundary": (
            "This derives a relational QFI for the specified two-intertwiner "
            "spin-network gluing channel and retained six-qubit boundary state. "
            "It is not yet a unique BQG dynamical edge weight: boundary tracing, "
            "holonomy gauge reduction, spin sums, and coarse-graining over the "
            "causal-to-PL fibres remain to be derived."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=16)
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
