#!/usr/bin/env python3
"""Three-node, two-link relational-QFI test for network non-additivity.

Three four-valent SU(2)-invariant spin-1/2 intertwiners form a path A-B-C.
The A3-B0 and B3-C0 half-links are contracted with normalized epsilon
singlets. Independent relative twists are inserted on B0 and B3 before
contraction. The resulting 8-qubit boundary state defines a 6x6 multiparameter
pure-state QFI matrix. This is a finite network control, not a full BQG
history-amplitude or continuum metric.
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


def path_boundary_amplitude(
    intertwiner_a: np.ndarray,
    intertwiner_b: np.ndarray,
    intertwiner_c: np.ndarray,
) -> np.ndarray:
    """Contract A3-B0 and B3-C0; retain A0..2, B1..2, C1..3."""
    return np.einsum(
        np.asarray(intertwiner_a).reshape((2,) * 4), [0, 1, 2, 3],
        np.asarray(intertwiner_b).reshape((2,) * 4), [4, 5, 6, 7],
        np.asarray(intertwiner_c).reshape((2,) * 4), [8, 9, 10, 11],
        EPSILON, [3, 4],
        EPSILON, [7, 8],
        [0, 1, 2, 5, 6, 9, 10, 11],
        optimize=True,
    ).reshape(-1)


def _apply_boundary_twist(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, site: int, pauli: np.ndarray
) -> np.ndarray:
    generator = node_generator(pauli, site)
    return path_boundary_amplitude(a, -1j * (generator @ b), c)


def network_relational_qfi(
    intertwiner_a: np.ndarray,
    intertwiner_b: np.ndarray,
    intertwiner_c: np.ndarray,
) -> tuple[float, np.ndarray, float]:
    """Return path gluing probability, 6x6 QFI, and Gauss pullthrough residual."""
    a, b, c = intertwiner_a, intertwiner_b, intertwiner_c
    raw = path_boundary_amplitude(a, b, c)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("two-link path contraction has zero support")
    state = raw / math.sqrt(probability)

    derivatives = []
    residuals = []
    # Two independent link twists, both represented on B's shared half-links.
    for shared_site, other_sites in ((0, (1, 2, 3)), (3, (0, 1, 2))):
        for pauli in PAULI:
            derivative_raw = _apply_boundary_twist(a, b, c, shared_site, pauli)
            derivative = derivative_raw / math.sqrt(probability)
            # B-node Gauss law: J_shared = -sum of the other three half-link J's.
            pulled_generator = sum(
                (node_generator(pauli, site) for site in other_sites),
                np.zeros((16, 16), dtype=complex),
            )
            pulled_raw = path_boundary_amplitude(
                a, -1j * (pulled_generator @ b), c
            )
            pulled_derivative = pulled_raw / math.sqrt(probability)
            residuals.append(float(np.linalg.norm(derivative - pulled_derivative)))
            derivatives.append(derivative)

    overlaps = [np.vdot(state, derivative) for derivative in derivatives]
    qfi = np.zeros((6, 6), dtype=float)
    for i, j in itertools.product(range(6), repeat=2):
        qfi[i, j] = 4.0 * float(np.real(
            np.vdot(derivatives[i], derivatives[j])
            - np.conjugate(overlaps[i]) * overlaps[j]
        ))
    return probability, qfi, max(residuals)


def random_intertwiner(rng: np.random.Generator) -> np.ndarray:
    coefficients = rng.normal(size=2) + 1j * rng.normal(size=2)
    coefficients /= np.linalg.norm(coefficients)
    return SQ.logical_basis() @ coefficients


def run(seed: int = 20261118, trials: int = 32) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required to test correlation variability")
    rng = np.random.default_rng(seed)
    probabilities = []
    matrices = []
    pullthrough = []
    cross_scalars = []
    diagonal_errors = []
    cross_isotropic_errors = []
    minimum_eigenvalues = []

    for _ in range(trials):
        a, b, c = (random_intertwiner(rng) for _ in range(3))
        probability, qfi, residual = network_relational_qfi(a, b, c)
        probabilities.append(probability)
        matrices.append(qfi)
        pullthrough.append(residual)

        diagonal_errors.append(max(
            np.linalg.norm(qfi[:3, :3] - np.eye(3)),
            np.linalg.norm(qfi[3:, 3:] - np.eye(3)),
        ))
        cross = qfi[:3, 3:]
        scalar = float(np.trace(cross) / 3.0)
        cross_scalars.append(scalar)
        cross_isotropic_errors.append(float(np.linalg.norm(cross - scalar * np.eye(3))))
        minimum_eigenvalues.append(float(np.linalg.eigvalsh(qfi)[0]))

    cross_abs = [abs(value) for value in cross_scalars]
    checks = {
        "two_singlet_gluing_probability_is_one_sixteenth": max(
            abs(value - 1.0 / 16.0) for value in probabilities
        ) < 1e-12,
        "each_link_diagonal_qfi_block_is_identity": max(diagonal_errors) < 1e-10,
        "cross_block_is_su2_isotropic_for_path_intertwiners": max(
            cross_isotropic_errors
        ) < 1e-10,
        "cross_link_correlations_are_not_identically_zero": max(cross_abs) > 1e-3,
        "cross_link_correlations_depend_on_intermediate_intertwiner": (
            max(cross_scalars) - min(cross_scalars) > 1e-2
        ),
        "full_qfi_is_positive_semidefinite": min(minimum_eigenvalues) > -1e-10,
        "gauss_law_twists_pull_through_each_shared_face": max(pullthrough) < 1e-10,
    }
    return {
        "status": "three-node two-link relational network QFI non-additivity control",
        "passed": all(checks.values()),
        "seed": seed,
        "trials": trials,
        "network": "path A-B-C; one normalized spin-1/2 epsilon singlet on A3-B0 and B3-C0",
        "boundary_state": "normalized pure state on eight retained half-links",
        "relative_twist_parameters": "three SU(2) generators on B0 and three on B3",
        "qfi_shape": [6, 6],
        "gluing_probability_min": min(probabilities),
        "gluing_probability_max": max(probabilities),
        "diagonal_link_qfi_max_error_from_identity": max(diagonal_errors),
        "cross_block_form": "c I_3 within numerical tolerance for the tested SU(2)-invariant path family",
        "cross_block_scalar_min": min(cross_scalars),
        "cross_block_scalar_max": max(cross_scalars),
        "cross_block_absolute_max": max(cross_abs),
        "cross_block_isotropy_max_error": max(cross_isotropic_errors),
        "full_qfi_minimum_eigenvalue": min(minimum_eigenvalues),
        "max_gauss_pullthrough_residual": max(pullthrough),
        "checks": checks,
        "interpretation": (
            "The individual 3x3 diagonal link blocks equal I_3, but the two "
            "links have state-dependent nonzero cross-QFI. Therefore summing "
            "independent two-node edge QFIs is not the multiparameter network "
            "QFI in general. A coarse graph metric must retain or justify "
            "integrating out these correlations."
        ),
        "claim_boundary": (
            "Finite three-node spin-1/2 path with pure invariant intertwiners "
            "and two singlet contractions only. No boundary gauge reduction, "
            "mixed-state tracing, spin sums, BQG history-amplitude measure, "
            "or continuum-limit result is claimed."
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
