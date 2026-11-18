#!/usr/bin/env python3
"""Higher-spin four-valent SU(2) intertwiner path and relational QFI gate.

Constructs Inv(V_j^tensor4) by an explicit nullspace of total generators,
glues a three-node path with normalized spin-j invariant bilinear tensors,
and computes the full 6x6 pure-boundary-state QFI. No external symbolic
or angular-momentum package is required.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from functools import lru_cache

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

PAIRS = ((0, 1), (1, 2), (2, 0))


def spin_generators(two_j: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer (j >= 1/2)")
    j = two_j / 2.0
    dimension = two_j + 1
    m_values = np.array([j - k for k in range(dimension)], dtype=float)
    raising = np.zeros((dimension, dimension), dtype=complex)
    for column, m in enumerate(m_values):
        if column > 0:
            raising[column - 1, column] = math.sqrt((j - m) * (j + m + 1.0))
    lowering = raising.conj().T
    return (
        0.5 * (raising + lowering),
        (raising - lowering) / (2.0j),
        np.diag(m_values).astype(complex),
    )


def invariant_singlet_tensor(two_j: int) -> np.ndarray:
    if not isinstance(two_j, int) or two_j < 1:
        raise ValueError("two_j must be a positive integer (j >= 1/2)")
    j = two_j / 2.0
    dimension = two_j + 1
    tensor = np.zeros((dimension, dimension), dtype=complex)
    for row in range(dimension):
        tensor[row, two_j - row] = (-1.0) ** row / math.sqrt(dimension)
    return tensor


@lru_cache(maxsize=None)
def four_valent_intertwiner_basis(two_j: int) -> tuple[np.ndarray, float]:
    """Orthonormal columns span the kernel of all three total node generators."""
    generators = spin_generators(two_j)
    dimension = two_j + 1
    identity = np.eye(dimension, dtype=complex)
    total = []
    for generator in generators:
        terms = []
        for leg in range(4):
            factors = [generator if k == leg else identity for k in range(4)]
            term = factors[0]
            for factor in factors[1:]:
                term = np.kron(term, factor)
            terms.append(term)
        total.append(sum(terms, np.zeros((dimension**4, dimension**4), complex)))
    stacked = np.vstack(total)
    _u, singular, vh = np.linalg.svd(stacked, full_matrices=True)
    threshold = max(stacked.shape) * np.finfo(float).eps * (singular[0] if singular.size else 1.0) * 50.0
    rank = int(np.sum(singular > threshold))
    basis = vh[rank:].conj().T
    residual = float(np.linalg.norm(stacked @ basis))
    expected_dimension = two_j + 1
    if basis.shape != (dimension**4, expected_dimension):
        raise AssertionError(
            f"unexpected intertwiner dimension {basis.shape[1]}, expected {expected_dimension}"
        )
    return basis, residual


def random_intertwiner(
    rng: np.random.Generator, basis: np.ndarray
) -> np.ndarray:
    coefficients = rng.normal(size=basis.shape[1]) + 1j * rng.normal(size=basis.shape[1])
    coefficients /= np.linalg.norm(coefficients)
    return basis @ coefficients


def node_generator(generator: np.ndarray, leg: int, dimension: int) -> np.ndarray:
    factors = [generator if k == leg else np.eye(dimension, dtype=complex) for k in range(4)]
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor)
    return result


def contract_path(
    states: list[np.ndarray],
    singlet: np.ndarray,
    twist: tuple[int, int, np.ndarray] | None = None,
) -> np.ndarray:
    dimension = singlet.shape[0]
    tensors = [np.asarray(state, complex).reshape((dimension,) * 4).copy() for state in states]
    if twist is not None:
        node, leg, generator = twist
        tensors[node] = (
            -1j * (node_generator(generator, leg, dimension) @ tensors[node].reshape(-1))
        ).reshape((dimension,) * 4)
    operands = [
        tensors[0], [0, 1, 2, 3],
        tensors[1], [4, 5, 6, 7],
        tensors[2], [8, 9, 10, 11],
        singlet, [3, 4],
        singlet, [7, 8],
        [0, 1, 2, 5, 6, 9, 10, 11],
    ]
    return np.einsum(*operands, optimize=True).reshape(-1)


def path_qfi(
    states: list[np.ndarray], singlet: np.ndarray, generators: tuple[np.ndarray, ...]
) -> tuple[float, np.ndarray, float]:
    raw = contract_path(states, singlet)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("higher-spin path contraction has zero support")
    state = raw / math.sqrt(probability)
    derivatives = []
    for leg in (0, 3):
        for generator in generators:
            derivative_raw = contract_path(states, singlet, (1, leg, generator))
            derivative = derivative_raw / math.sqrt(probability)
            derivative -= state * np.vdot(state, derivative)
            derivatives.append(derivative)
    overlaps = [np.vdot(state, derivative) for derivative in derivatives]
    qfi = np.zeros((6, 6), dtype=float)
    for i in range(6):
        for k in range(6):
            qfi[i, k] = 4.0 * float(np.real(
                np.vdot(derivatives[i], derivatives[k])
                - np.conjugate(overlaps[i]) * overlaps[k]
            ))
    return probability, qfi, float(np.linalg.norm(qfi - qfi.T))


def pair_dot_operator(generators: tuple[np.ndarray, ...]) -> np.ndarray:
    dimension = generators[0].shape[0]
    return sum(
        (node_generator(g, 0, dimension) @ node_generator(g, 3, dimension)
         for g in generators),
        np.zeros((dimension**4, dimension**4), dtype=complex),
    )


def analyze_spin(two_j: int, seed: int = 20261118, trials: int = 8) -> dict[str, object]:
    rng = np.random.default_rng(seed + two_j)
    generators = spin_generators(two_j)
    singlet = invariant_singlet_tensor(two_j)
    basis, intertwiner_residual = four_valent_intertwiner_basis(two_j)
    dimension = two_j + 1
    j = two_j / 2.0
    expected_diagonal = 4.0 * j * (j + 1.0) / 3.0
    pair_operator = pair_dot_operator(generators)
    projected_pair = basis.conj().T @ pair_operator @ basis
    pair_spectrum = np.linalg.eigvalsh(projected_pair)
    probabilities, diagonal_errors, cross_errors, pullthrough_errors = [], [], [], []
    min_qfi_eigenvalues, cross_values, analytic_cross_errors = [], [], []
    for _ in range(trials):
        states = [random_intertwiner(rng, basis) for _ in range(3)]
        probability, qfi, symmetry = path_qfi(states, singlet, generators)
        probabilities.append(probability)
        diagonal_errors.append(max(
            np.linalg.norm(qfi[:3, :3] - expected_diagonal * np.eye(3)),
            np.linalg.norm(qfi[3:, 3:] - expected_diagonal * np.eye(3)),
        ))
        cross = qfi[:3, 3:]
        scalar = float(np.trace(cross) / 3.0)
        cross_values.append(scalar)
        predicted_cross = (4.0 / 3.0) * float(
            np.vdot(states[1], pair_operator @ states[1]).real
        )
        analytic_cross_errors.append(abs(scalar - predicted_cross))
        cross_errors.append(float(np.linalg.norm(cross - scalar * np.eye(3))))
        min_qfi_eigenvalues.append(float(np.linalg.eigvalsh(qfi)[0]))
        # The shared-node Gauss law moves the twist to the other three legs.
        for leg in (0, 3):
            generator = generators[0]
            direct = contract_path(states, singlet, (1, leg, generator))
            other = -sum(
                (node_generator(generator, other_leg, dimension) @ states[1]
                 for other_leg in range(4) if other_leg != leg),
                np.zeros(dimension**4, dtype=complex),
            )
            pulled_states = [states[0], -1j * other, states[2]]
            pulled = contract_path(pulled_states, singlet)
            pullthrough_errors.append(float(np.linalg.norm(direct - pulled)))
        if symmetry > 1e-10:
            raise AssertionError("higher-spin path QFI is not symmetric")
    expected_pair_values = np.array([
        0.5 * (k * (k + 1.0) - 2.0 * j * (j + 1.0))
        for k in range(two_j + 1)
    ])
    return {
        "spin": j,
        "dimension": dimension,
        "intertwiner_dimension": int(basis.shape[1]),
        "intertwiner_nullspace_residual": intertwiner_residual,
        "singlet_normalization_error": abs(float(np.linalg.norm(singlet)) - 1.0),
        "singlet_lie_invariance_residual": max(
            float(np.linalg.norm(g @ singlet + singlet @ g.T)) for g in generators
        ),
        "gluing_probability_min": min(probabilities),
        "gluing_probability_max": max(probabilities),
        "expected_single_link_qfi_diagonal": expected_diagonal,
        "max_diagonal_block_error": max(diagonal_errors),
        "cross_block_isotropy_max_error": max(cross_errors),
        "cross_block_scalar_min": min(cross_values),
        "cross_block_scalar_max": max(cross_values),
        "cross_block_casimir_formula_max_error": max(analytic_cross_errors),
        "middle_pair_casimir_eigenvalues": pair_spectrum.tolist(),
        "expected_pair_casimir_eigenvalues": expected_pair_values.tolist(),
        "pair_casimir_spectrum_error": float(np.linalg.norm(pair_spectrum - expected_pair_values)),
        "max_gauss_pullthrough_residual": max(pullthrough_errors),
        "minimum_network_qfi_eigenvalue": min(min_qfi_eigenvalues),
        "cross_spectrum_interval_violation": max(analytic_cross_errors),
        "checks": {
            "intertwiner_dimension_is_2j_plus_1": bool(basis.shape[1] == two_j + 1),
            "intertwiner_nullspace_is_accurate": bool(intertwiner_residual < 2e-10),
            "singlet_is_normalized_and_invariant": bool(
                abs(np.linalg.norm(singlet) - 1.0) < 2e-12
                and max(np.linalg.norm(g @ singlet + singlet @ g.T) for g in generators) < 2e-10
            ),
            "single_link_qfi_is_isotropic_and_exact": bool(max(diagonal_errors) < 2e-9),
            "network_cross_block_is_isotropic": bool(max(cross_errors) < 2e-9),
            "network_cross_block_matches_middle_node_casimir": bool(max(analytic_cross_errors) < 2e-9),
            "pair_casimir_spectrum_is_exact": bool(np.linalg.norm(pair_spectrum - expected_pair_values) < 2e-9),
            "gauss_pullthrough_holds": bool(max(pullthrough_errors) < 2e-9),
            "network_qfi_is_psd": bool(min(min_qfi_eigenvalues) > -2e-9),
            "cross_block_respects_pair_casimir_bounds": bool(
                min(cross_values) >= (4.0 / 3.0) * pair_spectrum[0] - 2e-9
                and max(cross_values) <= (4.0 / 3.0) * pair_spectrum[-1] + 2e-9
            ),
        },
    }


def run(two_j_values: tuple[int, ...] = (1, 2, 3, 4), trials: int = 8) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rows = [analyze_spin(two_j, trials=trials) for two_j in two_j_values]
    return {
        "status": "higher-spin four-valent path relational QFI gate",
        "passed": all(all(row["checks"].values()) for row in rows),
        "tested_spins": [row["spin"] for row in rows],
        "trials_per_spin": trials,
        "exact_predictions": {
            "intertwiner_dimension": "2j+1",
            "single_link_diagonal_qfi": "4*j*(j+1)/3 times I3",
            "network_cross_block": "(4/3)*<J_middle,left dot J_middle,right> times I3",
            "pair_casimir_eigenvalues": "0.5*(k*(k+1)-2*j*(j+1)), k=0,...,2j",
        },
        "claim_boundary": (
            "Finite three-node path, pure four-valent SU(2) intertwiners and "
            "normalized invariant spin-j singlet contractions. No spin sum, "
            "mixed boundary reduction, BQG history amplitude, or continuum "
            "conductance is claimed."
        ),
        "representations": rows,
    }


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=8)
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
