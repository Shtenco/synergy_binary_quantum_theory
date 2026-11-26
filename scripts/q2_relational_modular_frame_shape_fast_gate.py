#!/usr/bin/env python3
"""Fast finite-field certificates for combined frame and node-shape rank.

The exact algebraic Jacobian is reduced directly modulo 97. The prime splits
sqrt(2), sqrt(3), and i; projective tangent contractions and determinants are
then evaluated in the resulting residue algebra. A nonzero real determinant
residue proves the exact algebraic minor is nonzero over characteristic zero.
"""
from __future__ import annotations

import json

import numpy as np
import sympy as sp

from q2_relational_graph_qfi_gate import graph_fixtures

P = 97
RT2, RT3, IMOD = 14, 10, 22
INV2, INVRT2, INVRT3 = pow(2, -1, P), pow(RT2, -1, P), pow(RT3, -1, P)


class C:
    """Element a+i*b in F_97[i], stored as a pair to retain conjugation."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        if isinstance(a, C):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = int(a) % P, int(b) % P

    @staticmethod
    def coerce(x):
        return x if isinstance(x, C) else C(x)

    def __add__(self, other):
        o = self.coerce(other)
        return C(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        o = self.coerce(other)
        return C(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inverse(self):
        norm = (self.a * self.a + self.b * self.b) % P
        if norm == 0:
            raise ZeroDivisionError("noninvertible residue in F_97[i]")
        ni = pow(norm, -1, P)
        return C(self.a * ni, -self.b * ni)

    def __truediv__(self, other):
        return self * self.coerce(other).inverse()

    def conjugate(self):
        return C(self.a, -self.b)

    def __repr__(self):
        return f"C({self.a},{self.b})"


ZERO, ONE, IMAG = C(0), C(1), C(0, 1)
PAULI = (
    np.array([[ZERO, ONE], [ONE, ZERO]], dtype=object),
    np.array([[ZERO, -IMAG], [IMAG, ZERO]], dtype=object),
    np.array([[ONE, ZERO], [ZERO, C(-1)]], dtype=object),
)

CERTIFICATES = {
    "path_3": {
        "coeffs": [[(1, 0), (-1, 2)], [(-1, 2), (2, -2)], [(-1, 0), (1, -1)]],
        "rows": [307,460,309,458,54,195,405,299,468,308,42,362,201,277,35,490,49,452,305,85,461,203,82,205,462,314,469],
        "cols": [38,39,40,36,37,41,6,7,8,30,31,32,18,19,20,14,12,13,0,4,2,34,33,35,24,16,29],
        "kernel": 15,
    },
    "cycle_4": {
        "coeffs": [[(1, 1), (2, 0)], [(1, 0), (0, 2)], [(2, 0), (-1, -1)], [(1, 0), (1, 1)]],
        "rows": [283,484,102,153,357,412,436,23,355,106,43,347,101,294,422,103,152,281,486,70,410,145,401,326,480,359,331,11,356,86,229,100],
        "cols": [54,55,52,53,50,51,48,49,9,11,13,36,38,4,14,12,34,33,35,10,45,47,46,25,24,26,37,5,21,8,16,15],
        "kernel": 24,
    },
}


def _rat_mod(x) -> int:
    x = sp.cancel(x)
    n, d = sp.fraction(x)
    n, d = int(n), int(d)
    if d % P == 0:
        raise ZeroDivisionError("prime 97 divides an exact denominator")
    return (n * pow(d % P, -1, P)) % P


def _logical_basis_mod():
    # Exact SU(2)-invariant four-spin basis, then reduction of its real entries.
    kron = sp.kronecker_product
    z, o = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    sing = (kron(z, o) - kron(o, z)) / sp.sqrt(2)
    l0 = kron(sing, sing)
    tp, t0, tm = kron(z, z), (kron(z, o) + kron(o, z)) / sp.sqrt(2), kron(o, o)
    l1 = (kron(tp, tm) - kron(t0, t0) + kron(tm, tp)) / sp.sqrt(3)
    basis = sp.Matrix.hstack(l0, l1)
    subst = {sp.sqrt(2): RT2, sp.sqrt(3): RT3, sp.sqrt(6): (RT2 * RT3) % P}
    return np.array([[_rat_mod(sp.re(basis[r, c]).subs(subst)) for c in range(2)]
                     for r in range(16)], dtype=object)


def _contract(states, edges):
    tensors = [np.asarray(x, dtype=object).reshape((2,) * 4) for x in states]
    labels = [[None] * 4 for _ in states]
    next_label = 0
    eps = np.array([[ZERO, C(INVRT2)], [C(-INVRT2), ZERO]], dtype=object)
    operands = []
    for u, lu, v, lv in edges:
        a, b = next_label, next_label + 1
        next_label += 2
        labels[u][lu], labels[v][lv] = a, b
        operands.extend((eps, [a, b]))
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


def _dot_conj(x, y):
    total = ZERO
    for a, b in zip(x, y):
        total += a.conjugate() * b
    return total


def _project(raw, derivative, norm2):
    overlap = _dot_conj(raw, derivative)
    return [d - r * (overlap / norm2) for r, d in zip(raw, derivative)]


def _det_mod(a):
    a = [row[:] for row in a]
    n = len(a)
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] % P), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        q = a[col][col] % P
        det = (det * q) % P
        qi = pow(q, -1, P)
        for row in range(col + 1, n):
            factor = (a[row][col] * qi) % P
            for k in range(col + 1, n):
                a[row][k] = (a[row][k] - factor * a[col][k]) % P
            a[row][col] = 0
    return det % P


def certify(name):
    spec = CERTIFICATES[name]
    n, edges = graph_fixtures()[name]
    basis = _logical_basis_mod()
    states = []
    for coeffs in spec["coeffs"]:
        c = [C(a, b) for a, b in coeffs]
        states.append(basis @ np.array(c, dtype=object))
    raw = _contract(states, edges)
    norm2 = _dot_conj(raw, raw)
    if norm2.b != 0 or norm2.a == 0:
        raise ZeroDivisionError(f"postselection norm is not invertible for {name}")
    tangent_cols = []

    # Frame variations, ordered node/leg/axis as in the T75 tangent gate.
    for node in range(n):
        for leg in range(4):
            for pauli in PAULI:
                ops = [np.eye(2, dtype=object) for _ in range(4)]
                ops[leg] = pauli
                gen = np.array([[C(INV2) * x for x in row]
                                for row in np.kron(np.kron(np.kron(ops[0], ops[1]), ops[2]), ops[3])], dtype=object)
                varied = list(states)
                varied[node] = np.array([(-IMAG) * x for x in gen @ states[node]], dtype=object)
                tangent_cols.append(_project(raw, _contract(varied, edges), norm2))

    # Logical ray directions u=(-conj(b),conj(a)) and i*u per node.
    for node, coeffs in enumerate(spec["coeffs"]):
        a, b = (C(*coeffs[0]), C(*coeffs[1]))
        u = [ -b.conjugate(), a.conjugate() ]
        for direction in (u, [IMAG * x for x in u]):
            varied = list(states)
            varied[node] = basis @ np.array(direction, dtype=object)
            tangent_cols.append(_project(raw, _contract(varied, edges), norm2))

    boundary_dim = len(raw)
    mod_matrix = []
    for row in spec["rows"]:
        rr = []
        for col in spec["cols"]:
            entry = tangent_cols[col][row] if row < boundary_dim else tangent_cols[col][row - boundary_dim]
            rr.append(entry.a if row < boundary_dim else entry.b)
        mod_matrix.append(rr)
    residue = _det_mod(mod_matrix)
    if residue == 0:
        raise ArithmeticError(f"zero determinant residue for {name}")
    quotient_dim = 12 * n - spec["kernel"]
    return {
        "fixture": name,
        "prime": P,
        "residues": {"sqrt2": RT2, "sqrt3": RT3, "i": IMOD},
        "minor_order": len(spec["rows"]),
        "frame_quotient_dimension": quotient_dim,
        "node_shape_dimension": 2 * n,
        "combined_dimension": quotient_dim + 2 * n,
        "determinant_mod_prime": residue,
        "exact_minor_nonzero": True,
    }


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixture", choices=["all", *CERTIFICATES], default="all")
    a = p.parse_args()
    names = list(CERTIFICATES) if a.fixture == "all" else [a.fixture]
    print(json.dumps({"certificates": [certify(name) for name in names]}, indent=2))


if __name__ == "__main__":
    main()

