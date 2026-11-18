#!/usr/bin/env python3
"""SU(2) invariant singlet bilinear form for arbitrary finite spin j.

For the spin-j irrep, construct the invariant two-index tensor C_j and
verify J_a C_j + C_j J_a^T = 0 and U_a C_j U_a^T = C_j. This proves the
local common-endpoint cancellation for integer and half-integer spins, but
does not construct higher-spin node intertwiners or network QFI.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def spin_generators(two_j: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not isinstance(two_j, int) or two_j < 0:
        raise ValueError("two_j must be a nonnegative integer")
    j = 0.5 * two_j
    dimension = two_j + 1
    m_values = np.array([j - index for index in range(dimension)], dtype=float)
    raising = np.zeros((dimension, dimension), dtype=complex)
    for column, m in enumerate(m_values):
        if column > 0:
            raising[column - 1, column] = math.sqrt((j - m) * (j + m + 1.0))
    lowering = raising.conj().T
    jx = 0.5 * (raising + lowering)
    jy = (raising - lowering) / (2.0j)
    jz = np.diag(m_values).astype(complex)
    return jx, jy, jz


def invariant_singlet_tensor(two_j: int) -> np.ndarray:
    if not isinstance(two_j, int) or two_j < 0:
        raise ValueError("two_j must be a nonnegative integer")
    j = 0.5 * two_j
    dimension = two_j + 1
    tensor = np.zeros((dimension, dimension), dtype=complex)
    for row in range(dimension):
        m = j - row
        column = two_j - row
        phase_exponent = int(round(j - m))
        tensor[row, column] = (-1.0) ** phase_exponent / math.sqrt(dimension)
    return tensor


def unitary_from_generator(generator: np.ndarray, angle: float) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(generator)
    return (eigenvectors * np.exp(-1.0j * angle * eigenvalues)) @ eigenvectors.conj().T


def analyze_spin(two_j: int) -> dict[str, object]:
    generators = spin_generators(two_j)
    singlet = invariant_singlet_tensor(two_j)
    dimension = two_j + 1
    lie_residuals = [
        float(np.linalg.norm(generator @ singlet + singlet @ generator.T))
        for generator in generators
    ]
    finite_residuals = []
    for axis, generator in enumerate(generators):
        unitary = unitary_from_generator(generator, 0.371 + 0.113 * axis)
        finite_residuals.append(float(np.linalg.norm(
            unitary @ singlet @ unitary.T - singlet
        )))
    norm_error = abs(float(np.linalg.norm(singlet)) - 1.0)
    transpose_parity_error = float(np.linalg.norm(
        singlet.T - ((-1.0) ** two_j) * singlet
    ))
    return {
        "two_j": two_j,
        "spin": 0.5 * two_j,
        "dimension": dimension,
        "singlet_norm_error": norm_error,
        "max_lie_algebra_invariance_residual": max(lie_residuals),
        "max_finite_group_invariance_residual": max(finite_residuals),
        "singlet_transpose_parity": "antisymmetric" if two_j % 2 else "symmetric",
        "transpose_parity_residual": transpose_parity_error,
        "passed": (
            norm_error < 2e-12
            and max(lie_residuals) < 2e-12
            and max(finite_residuals) < 2e-12
            and transpose_parity_error < 2e-12
        ),
    }


def run(two_j_values: tuple[int, ...] = (1, 2, 3, 4, 5, 6)) -> dict[str, object]:
    rows = [analyze_spin(two_j) for two_j in two_j_values]
    return {
        "status": "arbitrary-irrep local SU(2) singlet invariance control",
        "passed": all(row["passed"] for row in rows),
        "tested_spins": [row["spin"] for row in rows],
        "identity": "J_a C_j + C_j J_a^T = 0; U_j C_j U_j^T = C_j",
        "interpretation": (
            "The common-endpoint infinitesimal cancellation is representation "
            "independent for the invariant two-leg singlet tensor, including "
            "integer and half-integer SU(2) irreps. This is a local representation "
            "identity, not a higher-spin spin-network or QFI computation."
        ),
        "representations": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
