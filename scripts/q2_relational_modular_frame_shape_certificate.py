#!/usr/bin/env python3
"""Modular exact nonzero-minor certificates for the finite frame-shape bridge.

A nonzero determinant after reduction in a valid residue field proves the
corresponding algebraic determinant over QQ(sqrt(2),sqrt(3),i) is nonzero.
The input minors are the frozen witnesses from the companion response gate.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import sympy as sp

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
        "coefficients": [[(1, 0), (-1, 2)], [(-1, 2), (2, -2)], [(-1, 0), (1, -1)]],
        "rows": [307,460,309,458,54,195,405,299,468,308,42,362,201,277,35,490,49,452,305,85,461,203,82,205,462,314,469],
        "columns": [38,39,40,36,37,41,6,7,8,30,31,32,18,19,20,14,12,13,0,4,2,34,33,35,24,16,29],
        "kernel_dimension": 15,
    },
    "cycle_4": {
        "coefficients": [[(1, 1), (2, 0)], [(1, 0), (0, 2)], [(2, 0), (-1, -1)], [(1, 0), (1, 1)]],
        "rows": [283,484,102,153,357,412,436,23,355,106,43,347,101,294,422,103,152,281,486,70,410,145,401,326,480,359,331,11,356,86,229,100],
        "columns": [54,55,52,53,50,51,48,49,9,11,13,36,38,4,14,12,34,33,35,10,45,47,46,25,24,26,37,5,21,8,16,15],
        "kernel_dimension": 24,
    },
}

# The prime 97 splits all generators: 14^2=2, 10^2=3, 22^2=-1 (mod 97).
PRIME = 97
ROOTS = {SQRT(2): 14, SQRT(3): 10, SQRT(6): 43, I: 22}


def _logical_basis() -> sp.Matrix:
    kron = sp.kronecker_product
    z, o = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    singlet = (kron(z, o) - kron(o, z)) / SQRT(2)
    l0 = kron(singlet, singlet)
    tp, t0, tm = kron(z, z), (kron(z, o) + kron(o, z)) / SQRT(2), kron(o, o)
    l1 = (kron(tp, tm) - kron(t0, t0) + kron(tm, tp)) / SQRT(3)
    return sp.Matrix.hstack(l0, l1)


def _c(pair: tuple[int, int]) -> sp.Expr:
    return sp.Integer(pair[0]) + I * sp.Integer(pair[1])


def _contract(states: list[sp.Matrix], edges: list[tuple[int, int, int, int]]) -> np.ndarray:
    tensors = [np.array(state, dtype=object).reshape((2,) * 4) for state in states]
    labels: list[list[int | None]] = [[None] * 4 for _ in states]
    next_label = 0
    epsilon = np.array([[0, 1], [-1, 0]], dtype=object) / SQRT(2)
    operands: list[object] = []
    for u, lu, v, lv in edges:
        a, b = next_label, next_label + 1
        next_label += 2
        labels[u][lu], labels[v][lv] = a, b
        operands.extend((epsilon, [a, b]))
    boundary = []
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


def _mod_rational(expr: sp.Expr) -> int:
    value = sp.cancel(sp.expand(expr).subs(ROOTS))
    numerator, denominator = sp.fraction(value)
    numerator, denominator = sp.Integer(numerator), sp.Integer(denominator)
    n, d = int(numerator) % PRIME, int(denominator) % PRIME
    if d == 0:
        raise ZeroDivisionError("the selected prime divides an algebraic denominator")
    return (n * pow(d, -1, PRIME)) % PRIME


def _mod_det(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    determinant = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] % PRIME), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            determinant = -determinant
        value = a[col][col] % PRIME
        determinant = (determinant * value) % PRIME
        inverse = pow(value, -1, PRIME)
        for row in range(col + 1, n):
            factor = (a[row][col] * inverse) % PRIME
            if factor:
                for k in range(col + 1, n):
                    a[row][k] = (a[row][k] - factor * a[col][k]) % PRIME
                a[row][col] = 0
    return determinant % PRIME


def certify(name: str) -> dict[str, object]:
    spec = CERTIFICATES[name]
    node_count, edges = graph_fixtures()[name]
    logical = _logical_basis()
    states = [logical * sp.Matrix([_c(x) for x in row]) for row in spec["coefficients"]]
    raw = _contract(states, edges)
    norm2 = sum(sp.conjugate(x) * x for x in raw)
    derivatives: list[list[sp.Expr]] = []

    # Frame coordinates use the T75 ordering: node, leg, Pauli axis.
    for node in range(node_count):
        for leg in range(4):
            for pauli in PAULI:
                ops = [sp.eye(2)] * 4
                ops[leg] = pauli
                varied = list(states)
                varied[node] = (-I / 2) * sp.kronecker_product(*ops) * states[node]
                d = _contract(varied, edges)
                overlap = sum(sp.conjugate(raw[k]) * d[k] for k in range(len(raw)))
                derivatives.append([d[k] - raw[k] * overlap / norm2 for k in range(len(raw))])

    # Two exact real tangent directions for each logical ray.
    for node, (a_pair, b_pair) in enumerate(spec["coefficients"]):
        a, b = _c(a_pair), _c(b_pair)
        u = sp.Matrix([-sp.conjugate(b), sp.conjugate(a)])
        for direction in (u, I * u):
            varied = list(states)
            varied[node] = logical * direction
            d = _contract(varied, edges)
            overlap = sum(sp.conjugate(raw[k]) * d[k] for k in range(len(raw)))
            derivatives.append([d[k] - raw[k] * overlap / norm2 for k in range(len(raw))])

    boundary_dimension = len(raw)
    residues = []
    for row in spec["rows"]:
        residue_row = []
        for col in spec["columns"]:
            value = (sp.re(derivatives[col][row]) if row < boundary_dimension
                     else sp.im(derivatives[col][row - boundary_dimension]))
            residue_row.append(_mod_rational(value))
        residues.append(residue_row)
    residue_det = _mod_det(residues)
    if residue_det == 0:
        raise ArithmeticError(f"the selected residue of the {name} minor vanished")
    quotient_dim = 12 * node_count - spec["kernel_dimension"]
    return {
        "fixture": name,
        "prime": PRIME,
        "sqrt2_residue": ROOTS[SQRT(2)],
        "sqrt3_residue": ROOTS[SQRT(3)],
        "i_residue": ROOTS[I],
        "minor_order": len(spec["rows"]),
        "frame_quotient_dimension": quotient_dim,
        "node_shape_dimension": 2 * node_count,
        "combined_dimension": quotient_dim + 2 * node_count,
        "rows_zero_based": spec["rows"],
        "columns_zero_based": spec["columns"],
        "determinant_mod_prime": residue_det,
        "exact_determinant_is_nonzero": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", choices=["all", *CERTIFICATES], default="all")
    args = parser.parse_args()
    names = list(CERTIFICATES) if args.fixture == "all" else [args.fixture]
    print(json.dumps({"certificates": [certify(name) for name in names]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

