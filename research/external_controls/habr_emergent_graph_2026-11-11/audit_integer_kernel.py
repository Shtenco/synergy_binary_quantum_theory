#!/usr/bin/env python3
"""Exact integer-kernel audit of draft/stress_kernel.py.

Requires SymPy. It compares:
(1) the script's rational-nullspace -> primitive-row scaling procedure;
(2) a saturated Z-kernel basis obtained from Smith normal decomposition;
(3) LLL reduction of both lattices.

The matrix is transcribed from the source blob
0fa851abb181a121da11ae309fdc6f0de2f88f21.
"""
from __future__ import annotations
from functools import reduce
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp, smith_normal_form
from sympy.polys.domains import ZZ

NAMES = [
    "hbar", "h", "c", "lP", "tP", "EP", "G", "mP", "TP", "kB",
    "alpha", "me", "mp", "eps0", "mu0", "qe", "Rinf", "a0", "Z0",
    "Phi0", "lambda_e", "lambda_p", "vH", "muon", "tau", "pion_pm",
    "pion0", "Lambda", "kappa", "mH_mW"
]

M = sp.Matrix([
[0,0,0,0,3,sp.Rational(-1,3)],
[2,0,1,0,3,sp.Rational(-1,3)],
[0,0,1,-1,4,0],
[4,0,0,1,2,sp.Rational(-1,3)],
[4,0,-1,2,-2,sp.Rational(-1,3)],
[-8,-6,1,-2,5,0],
[-6,-10,3,-1,13,sp.Rational(-1,3)],
[0,2,-1,0,-3,0],
[6,0,1,0,-4,sp.Rational(1,3)],
[-6,0,-2,0,8,0],
[2,0,-1,2,-1,0],
[0,0,1,0,4,sp.Rational(-1,3)],
[0,0,sp.Rational(1,2),0,6,sp.Rational(-1,3)],
[-6,0,-3,-1,-20,sp.Rational(1,3)],
[6,0,1,3,12,sp.Rational(-1,3)],
[0,0,-1,0,-7,0],
[4,0,-1,3,3,0],
[-6,0,-1,-1,-4,0],
[6,0,2,2,16,sp.Rational(-1,3)],
[0,0,2,0,10,sp.Rational(-1,3)],
[-2,0,-1,1,-5,0],
[2,0,sp.Rational(-1,2),1,-7,0],
[5,0,sp.Rational(3,2),0,6,sp.Rational(-1,3)],
[4,-1,2,0,5,sp.Rational(-1,3)],
[0,0,sp.Rational(1,2),0,5,sp.Rational(-1,3)],
[-1,0,-2,0,6,sp.Rational(-1,3)],
[2,0,1,0,4,sp.Rational(-1,3)],
[0,0,sp.Rational(-1,2),0,12,sp.Rational(-2,3)],
[8,6,0,3,-3,sp.Rational(-1,3)],
[2,0,-1,0,0,0],
])


def primitive_integer_vector(vec):
    den = sp.ilcm(*[sp.denom(x) for x in vec])
    ints = [int(x*den) for x in vec]
    g = reduce(gcd, (abs(x) for x in ints if x), 0)
    return [x//g for x in ints]


def main():
    assert M.shape == (30, 6)
    assert len(NAMES) == M.rows
    assert M.rank() == 6

    # Reproduce the source's row-wise denominator clearing.
    naive_rows = [primitive_integer_vector(v) for v in M.T.nullspace()]
    K_naive = sp.Matrix(naive_rows)
    assert K_naive.shape == (24, 30)
    assert (M.T * K_naive.T).is_zero_matrix

    # Index of the source-generated sublattice in its saturation.
    S_naive = smith_normal_form(K_naive, domain=ZZ)
    inv_naive = [abs(int(S_naive[i,i])) for i in range(24)]
    index = 1
    for d in inv_naive:
        index *= d

    # A=6M^T is integer. Smith decomposition gives U*A*V=D,
    # with V unimodular. The last 24 columns of V span the full Z-kernel.
    A = (6*M.T).applyfunc(int)
    D, U, V = smith_normal_decomp(A, domain=ZZ)
    K_sat = V[:, 6:].T
    assert K_sat.shape == (24, 30)
    assert (A * K_sat.T).is_zero_matrix
    assert abs(int(V.det())) == 1

    # Express the naive basis in the saturated basis and verify index.
    _, pivots = K_sat.rref()
    cols = list(pivots)
    C = K_naive[:, cols] * K_sat[:, cols].inv()
    assert all(x.q == 1 for x in C)
    assert C*K_sat == K_naive
    assert abs(int(C.det())) == index

    # LLL both lattices under the same delta.
    L_naive = K_naive.lll(delta=sp.Rational(3,4))
    L_sat = K_sat.lll(delta=sp.Rational(3,4))

    def stats(B):
        rows = [[int(x) for x in row] for row in B.tolist()]
        norms = [sum(x*x for x in row) for row in rows]
        maxabs = [max(abs(x) for x in row) for row in rows]
        return {
            "basis_rows": len(rows),
            "max_abs_entry": max(maxabs),
            "rows_max_abs_le_10": sum(x <= 10 for x in maxabs),
            "norm_sq_min": min(norms),
            "norm_sq_max": max(norms),
            "all_exact_kernel": all((A*sp.Matrix(row)) == sp.zeros(6,1) for row in rows),
        }

    print("matrix_shape=", M.shape, "rank=", M.rank(), "rational_nullity=", 30-M.rank())
    print("SNF_diagonal=", [int(D[i,i]) for i in range(6)])
    print("naive_kernel_basis_shape=", K_naive.shape)
    print("naive_lattice_index_in_full_integer_kernel=", index)
    print("naive_LLL=", stats(L_naive))
    print("saturated_LLL=", stats(L_sat))
    print("saturated_LLL_basis_vectors:")
    for i, row in enumerate(L_sat.tolist(), 1):
        terms = [(NAMES[j], int(c)) for j,c in enumerate(row) if c]
        print(i, "norm2=", sum(int(c)**2 for c in row), terms)


if __name__ == "__main__":
    main()
