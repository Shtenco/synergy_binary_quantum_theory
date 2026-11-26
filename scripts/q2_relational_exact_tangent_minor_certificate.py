#!/usr/bin/env python3
"""Exact nonzero-minor certificates for the finite T75 graph fixtures.

The certificate uses fixed rows/columns selected during the numerical audit,
but reconstructs the logical intertwiners, edge contractions, projective
tangent, and determinants exactly in QQ(sqrt(2), sqrt(3), i). A nonzero minor
plus the exact Gauss/edge-singlet kernel count proves generic full quotient
rank on each fixed graph fixture. It does not prove a refinement-uniform gap.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from q2_relational_graph_qfi_gate import graph_fixtures


I = sp.I
SQRT = sp.sqrt
PAULI = (
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -I], [I, 0]]),
    sp.Matrix([[1, 0], [0, -1]]),
)

CERTIFICATES = {
    "path_3": {
        "coefficients": [
            [(1, 0), (-1, 2)],
            [(-1, 2), (2, -2)],
            [(-1, 0), (1, -1)],
        ],
        "rows": [307, 85, 52, 45, 203, 476, 492, 75, 459, 311, 362, 59,
                 277, 55, 179, 396, 54, 42, 213, 317, 314],
        "columns": [0, 1, 2, 30, 31, 32, 18, 19, 20, 4, 3, 5, 33, 34, 35,
                    6, 13, 8, 23, 22, 15],
        "quotient_dimension": 21,
        "independent_kernel_dimension": 15,
    },
    "cycle_4": {
        "coefficients": [
            [(0, 1), (0, -2)],
            [(1, 0), (0, 2)],
            [(2, 0), (-1, -1)],
            [(1, 0), (1, 1)],
        ],
        "rows": [102, 331, 347, 43, 294, 23, 103, 152, 356, 411, 286, 70,
                 145, 401, 326, 480, 283, 155, 11, 26, 357, 281, 282, 89],
        "columns": [3, 4, 5, 13, 12, 14, 34, 33, 35, 9, 10, 11, 45, 47, 46,
                    24, 25, 26, 37, 21, 23, 15, 16, 17],
        "quotient_dimension": 24,
        "independent_kernel_dimension": 24,
    },
}


def _logical_basis() -> sp.Matrix:
    zero, one = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    kron = sp.kronecker_product
    singlet = (kron(zero, one) - kron(one, zero)) / SQRT(2)
    l0 = kron(singlet, singlet)
    triplet_plus = kron(zero, zero)
    triplet_zero = (kron(zero, one) + kron(one, zero)) / SQRT(2)
    triplet_minus = kron(one, one)
    l1 = (
        kron(triplet_plus, triplet_minus)
        - kron(triplet_zero, triplet_zero)
        + kron(triplet_minus, triplet_plus)
    ) / SQRT(3)
    return sp.Matrix.hstack(l0, l1)


def _coefficient(pair: tuple[int, int]) -> sp.Expr:
    return sp.Integer(pair[0]) + I * sp.Integer(pair[1])


def _contract(states: list[sp.Matrix], edges: list[tuple[int, int, int, int]]) -> np.ndarray:
    tensors = [np.array(state, dtype=object).reshape((2,) * 4) for state in states]
    labels: list[list[int | None]] = [[None] * 4 for _ in states]
    next_label = 0
    epsilon = np.array([[0, 1], [-1, 0]], dtype=object) / SQRT(2)
    operands: list[object] = []
    for u, leg_u, v, leg_v in edges:
        label_u, label_v = next_label, next_label + 1
        next_label += 2
        labels[u][leg_u] = label_u
        labels[v][leg_v] = label_v
        operands.extend((epsilon, [label_u, label_v]))

    boundary: list[tuple[int, int]] = []
    for node in range(len(states)):
        for leg in range(4):
            if labels[node][leg] is None:
                labels[node][leg] = next_label
                boundary.append((node, leg))
                next_label += 1
    for node, tensor in enumerate(tensors):
        operands.extend((tensor, labels[node]))
    operands.append([labels[node][leg] for node, leg in boundary])
    return np.einsum(*operands, optimize="greedy").reshape(-1)


def exact_minor(name: str) -> dict[str, object]:
    spec = CERTIFICATES[name]
    node_count, edges = graph_fixtures()[name]
    logical = _logical_basis()
    states = [
        logical * sp.Matrix([_coefficient(pair) for pair in node_coefficients])
        for node_coefficients in spec["coefficients"]
    ]
    raw = _contract(states, edges)
    probability = sum(sp.conjugate(value) * value for value in raw)

    tangent_columns: list[list[sp.Expr]] = []
    for node in range(node_count):
        for leg in range(4):
            for pauli in PAULI:
                local_ops = [sp.eye(2)] * 4
                local_ops[leg] = pauli
                generator = sp.kronecker_product(*local_ops)
                varied = list(states)
                varied[node] = (-I / 2) * generator * states[node]
                derivative = _contract(varied, edges)
                overlap = sum(
                    sp.conjugate(raw[k]) * derivative[k]
                    for k in range(len(raw))
                )
                tangent_columns.append([
                    derivative[k] - raw[k] * overlap / probability
                    for k in range(len(raw))
                ])

    rows = spec["rows"]
    columns = spec["columns"]
    boundary_dimension = len(raw)
    minor = sp.Matrix([
        [
            sp.re(tangent_columns[column][row])
            if row < boundary_dimension
            else sp.im(tangent_columns[column][row - boundary_dimension])
            for column in columns
        ]
        for row in rows
    ])
    field = sp.QQ.algebraic_field(SQRT(2), SQRT(3), I)
    determinant = DomainMatrix.from_Matrix(minor).convert_to(field).det()
    if determinant.is_zero:
        raise ArithmeticError(f"the exact {name} minor is zero")

    return {
        "fixture": name,
        "minor_order": len(rows),
        "quotient_dimension": spec["quotient_dimension"],
        "independent_kernel_dimension": spec["independent_kernel_dimension"],
        "rows_zero_based": rows,
        "columns_zero_based": columns,
        "determinant_exact": str(field.to_sympy(determinant)),
        "determinant_decimal": str(sp.N(field.to_sympy(determinant), 30)),
        "exactly_nonzero": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", choices=["all", *CERTIFICATES], default="all")
    args = parser.parse_args()
    names = list(CERTIFICATES) if args.fixture == "all" else [args.fixture]
    print(json.dumps({"certificates": [exact_minor(name) for name in names]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

