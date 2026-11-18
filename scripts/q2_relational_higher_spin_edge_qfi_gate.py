#!/usr/bin/env python3
"""Higher-spin two-node singlet gluing and relational QFI.

Build random normalized four-valent SU(2) invariant node states by computing
the common kernel of the four-leg total generators. Glue one leg from each
node with the normalized spin-j singlet and compute the normalized boundary
ray QFI. For an invariant node, Schur's lemma gives rho_leg=I/(2j+1), hence
F_ab=4 j(j+1)/3 delta_ab and gluing probability=(2j+1)^(-2).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from q2_su2_irrep_singlet_invariance_gate import (
    invariant_singlet_tensor,
    spin_generators,
)


def _kron_on_leg(operator: np.ndarray, leg: int, dimension: int, valence: int = 4) -> np.ndarray:
    factors = [np.eye(dimension, dtype=complex) for _ in range(valence)]
    factors[leg] = operator
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor)
    return result


def invariant_node_basis(two_j: int, tolerance: float = 1e-10) -> tuple[np.ndarray, dict[str, float]]:
    generators = spin_generators(two_j)
    dimension = two_j + 1
    total_generators = [
        sum(
            (_kron_on_leg(generator, leg, dimension) for leg in range(4)),
            np.zeros((dimension**4, dimension**4), dtype=complex),
        )
        for generator in generators
    ]
    stacked = np.vstack(total_generators)
    _u, singular, vh = np.linalg.svd(stacked, full_matrices=True)
    scale = float(singular[0]) if singular.size else 0.0
    cutoff = tolerance * max(1.0, scale)
    rank = int(np.sum(singular > cutoff))
    basis = vh.conj().T[:, rank:]
    expected_dimension = two_j + 1
    if basis.shape != (dimension**4, expected_dimension):
        raise AssertionError(
            f"unexpected invariant multiplicity {basis.shape[1]} != {expected_dimension}"
        )
    residual = max(float(np.linalg.norm(generator @ basis)) for generator in total_generators)
    orthogonality = float(np.linalg.norm(basis.conj().T @ basis - np.eye(expected_dimension)))
    return basis, {
        "total_generator_kernel_rank": dimension**4 - rank,
        "max_gauss_residual": residual,
        "basis_orthogonality_error": orthogonality,
        "smallest_nonzero_singular_value": float(singular[rank - 1]) if rank else 0.0,
        "largest_null_singular_value": float(singular[rank]) if rank < len(singular) else 0.0,
    }


def random_invariant_node(
    rng: np.random.Generator, basis: np.ndarray
) -> np.ndarray:
    coefficients = rng.normal(size=basis.shape[1]) + 1j * rng.normal(size=basis.shape[1])
    coefficients /= np.linalg.norm(coefficients)
    state = basis @ coefficients
    return state / np.linalg.norm(state)


def apply_generator_to_leg(
    state: np.ndarray, generator: np.ndarray, leg: int, dimension: int
) -> np.ndarray:
    tensor = np.asarray(state, dtype=complex).reshape((dimension,) * 4)
    moved = np.moveaxis(tensor, leg, 0)
    acted = np.tensordot(generator, moved, axes=(1, 0))
    return np.moveaxis(acted, 0, leg).reshape(-1)


def glue_two_nodes(
    node_a: np.ndarray,
    node_b: np.ndarray,
    singlet: np.ndarray,
    dimension: int,
) -> np.ndarray:
    a = np.asarray(node_a, dtype=complex).reshape((dimension,) * 4)
    b = np.asarray(node_b, dtype=complex).reshape((dimension,) * 4)
    return np.einsum("abcx,ydef,xy->abcdef", a, b, singlet, optimize=True).reshape(-1)


def projective_derivative(raw: np.ndarray, derivative_raw: np.ndarray, probability: float) -> np.ndarray:
    state = raw / math.sqrt(probability)
    derivative = derivative_raw / math.sqrt(probability)
    return derivative - state * np.vdot(state, derivative)


def analyze_sample(
    two_j: int,
    node_a: np.ndarray,
    node_b: np.ndarray,
    generators: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict[str, object]:
    dimension = two_j + 1
    singlet = invariant_singlet_tensor(two_j)
    raw = glue_two_nodes(node_a, node_b, singlet, dimension)
    probability = float(np.vdot(raw, raw).real)
    if probability <= 1e-14:
        raise ValueError("zero-support higher-spin singlet gluing")
    derivatives_a, derivatives_b = [], []
    for generator in generators:
        a_acted = apply_generator_to_leg(node_a, generator, 3, dimension)
        b_acted = apply_generator_to_leg(node_b, generator, 0, dimension)
        derivatives_a.append(projective_derivative(
            raw, glue_two_nodes(a_acted, node_b, singlet, dimension), probability
        ))
        derivatives_b.append(projective_derivative(
            raw, glue_two_nodes(node_a, b_acted, singlet, dimension), probability
        ))
    da, db = np.column_stack(derivatives_a), np.column_stack(derivatives_b)
    qfi = 4.0 * np.real(da.conj().T @ da)
    qfi = 0.5 * (qfi + qfi.T)
    qfi_b = 4.0 * np.real(db.conj().T @ db)
    qfi_b = 0.5 * (qfi_b + qfi_b.T)
    expected = (4.0 * (0.5 * two_j) * (0.5 * two_j + 1.0) / 3.0) * np.eye(3)
    tangent_sum_residual = float(np.linalg.norm(da + db))
    qfi_eigenvalues = np.linalg.eigvalsh(qfi)
    checks = {
        "gluing_probability_is_inverse_dimension_squared": bool(
            np.isclose(probability, 1.0 / dimension**2, atol=3e-10, rtol=3e-10)
        ),
        "first_endpoint_qfi_matches_casimir_formula": bool(
            np.allclose(qfi, expected, atol=3e-9, rtol=3e-9)
        ),
        "second_endpoint_qfi_matches_first": bool(
            np.allclose(qfi_b, qfi, atol=3e-9, rtol=3e-9)
        ),
        "singlet_endpoint_derivatives_are_opposite": bool(tangent_sum_residual < 3e-9),
        "qfi_is_positive_definite": bool(qfi_eigenvalues[0] > 0.0),
    }
    return {
        "spin": 0.5 * two_j,
        "dimension": dimension,
        "gluing_probability": probability,
        "expected_gluing_probability": 1.0 / dimension**2,
        "qfi_matrix": qfi.tolist(),
        "expected_qfi_diagonal": float(expected[0, 0]),
        "qfi_minimum_eigenvalue": float(qfi_eigenvalues[0]),
        "qfi_maximum_eigenvalue": float(qfi_eigenvalues[-1]),
        "endpoint_tangent_sum_residual": tangent_sum_residual,
        "qfi_second_endpoint_error": float(np.linalg.norm(qfi_b - qfi)),
        "checks": checks,
        "passed": all(checks.values()),
    }


def run(seed: int = 20261118, trials: int = 3, two_j_values: tuple[int, ...] = (1, 2, 3)) -> dict[str, object]:
    if trials < 2:
        raise ValueError("at least two trials are required")
    rng = np.random.default_rng(seed)
    results = {}
    for two_j in two_j_values:
        basis, basis_diagnostics = invariant_node_basis(two_j)
        generators = spin_generators(two_j)
        rows = []
        for _ in range(trials):
            node_a = random_invariant_node(rng, basis)
            node_b = random_invariant_node(rng, basis)
            rows.append(analyze_sample(two_j, node_a, node_b, generators))
        results[str(0.5 * two_j)] = {
            "spin": 0.5 * two_j,
            "invariant_space_dimension": int(basis.shape[1]),
            "basis_diagnostics": basis_diagnostics,
            "trials": trials,
            "max_gluing_probability_error": max(
                abs(row["gluing_probability"] - row["expected_gluing_probability"])
                for row in rows
            ),
            "max_qfi_casimir_formula_error": max(
                float(np.linalg.norm(
                    np.asarray(row["qfi_matrix"]) -
                    row["expected_qfi_diagonal"] * np.eye(3)
                )) for row in rows
            ),
            "max_endpoint_tangent_sum_residual": max(
                row["endpoint_tangent_sum_residual"] for row in rows
            ),
            "checks": {
                key: all(row["checks"][key] for row in rows)
                for key in rows[0]["checks"]
            },
            "passed": all(row["passed"] for row in rows)
                and basis_diagnostics["max_gauss_residual"] < 3e-9
                and basis_diagnostics["basis_orthogonality_error"] < 3e-9,
            "samples": rows,
        }
    return {
        "status": "higher-spin two-node relational QFI control",
        "passed": all(row["passed"] for row in results.values()),
        "seed": seed,
        "trials_per_spin": trials,
        "analytic_prediction": {
            "gluing_probability": "1/(2j+1)^2",
            "relational_qfi": "4*j*(j+1)/3 * I_3",
            "endpoint_tangent_relation": "D_v=-D_u",
        },
        "claim_boundary": (
            "Higher-spin two-node gluing is computed with explicit four-valent "
            "invariant node states for the listed spins. This is not yet a "
            "higher-spin multi-edge graph contraction, spin foam sum, or "
            "history-weighted BQG conductance."
        ),
        "spins": results,
    }


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    result = run(seed=args.seed, trials=args.trials)
    rendered = json.dumps(result, indent=2, default=_json_default)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
