#!/usr/bin/env python3
"""Saturated integer-kernel audit of the first 30 rows in
physics_constants/constant_matrix_analysis.py.

This is the newer 30-observable matrix, expressed in the basis
(sqrt(2), sqrt(3), pi, lnK, lnN, N^(1/3)).
Requires SymPy with smith_normal_decomp support.
"""
from __future__ import annotations
from functools import reduce
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp, smith_normal_form
from sympy.polys.domains import ZZ

NAMES = [
"hbar","h","c","lP","tP","EP","G","mP","TP","kB","alpha","me","mp",
"eps0","mu0","qe","Rinf","a0","Z0","Phi0","lambda_e","lambda_p","vH",
"muon","tau","pion_pm","pion0","Lambda","kappa","mH_mW"
]
M = sp.Matrix([
[-2,-2,0,0,3,-1],
[0,-2,1,0,3,-1],
[-4,-4,1,-1,4,0],
[4,0,0,1,2,-1],
[8,4,-1,2,-2,-1],
[-10,-6,1,-2,5,0],
[-2,-10,3,-1,13,-1],
[-2,2,-1,0,-3,0],
[6,0,1,0,-4,1],
[-6,0,-2,0,8,-1],
[2,0,-1,2,-1,0],
[3,-1,1,0,4,-1],
[-3,-3,sp.Rational(1,2),0,6,-1],
[-6,0,-3,-1,-20,1],
[14,8,1,3,12,-1],
[-3,-3,-1,0,-7,0],
[1,-3,-1,3,3,0],
[-3,3,-1,-1,-4,0],
[10,4,2,2,16,-1],
[1,1,2,0,10,-1],
[1,3,-1,1,-5,0],
[7,5,sp.Rational(-1,2),1,-7,0],
[5,0,sp.Rational(3,2),0,6,-1],
[2,-3,2,0,5,-1],
[4,4,sp.Rational(1,2),0,5,-1],
[-5,0,-2,0,6,-1],
[8,6,1,0,4,-1],
[0,0,sp.Rational(-1,2),0,12,-2],
[20,6,0,3,-3,-1],
[3,1,-1,0,0,0],
])


def primitive(vec):
    den=sp.ilcm(*[sp.denom(x) for x in vec])
    vals=[int(x*den) for x in vec]
    g=reduce(gcd,(abs(x) for x in vals if x),0)
    return [x//g for x in vals]


def main():
    assert M.shape==(30,6) and len(NAMES)==30 and M.rank()==6
    naive=sp.Matrix([primitive(v) for v in M.T.nullspace()])
    assert naive.shape==(24,30)
    assert (M.T*naive.T).is_zero_matrix

    # Denominator 2 clears all rational entries.
    A=(2*M.T).applyfunc(int)
    D,U,V=smith_normal_decomp(A,domain=ZZ)
    saturated=V[:,6:].T
    assert (A*saturated.T).is_zero_matrix
    assert abs(int(V.det()))==1

    S=smith_normal_form(naive,domain=ZZ)
    inv=[abs(int(S[i,i])) for i in range(24)]
    index=1
    for d in inv: index*=d

    _, piv=saturated.rref()
    cols=list(piv)
    C=naive[:,cols]*saturated[:,cols].inv()
    assert all(x.q==1 for x in C)
    assert C*saturated==naive
    assert abs(int(C.det()))==index

    lll_naive=naive.lll(delta=sp.Rational(3,4))
    lll_full=saturated.lll(delta=sp.Rational(3,4))

    def report(B):
        rows=[[int(x) for x in row] for row in B.tolist()]
        return (max(max(abs(x) for x in row) for row in rows),
                sum(max(abs(x) for x in row)<=10 for row in rows),
                all((A*sp.Matrix(row))==sp.zeros(6,1) for row in rows))
    print("matrix_shape",M.shape,"rank",M.rank(),"nullity",24)
    print("SNF_diagonal",[int(D[i,i]) for i in range(6)])
    print("naive_lattice_index",index)
    print("naive_LLL_max_abs,rows_le_10,exact",report(lll_naive))
    print("saturated_LLL_max_abs,rows_le_10,exact",report(lll_full))
    print("saturated_LLL_vectors:")
    for i,row in enumerate(lll_full.tolist(),1):
        v=[int(x) for x in row]
        print(i,"norm2",sum(x*x for x in v),
              [(NAMES[j],x) for j,x in enumerate(v) if x])


if __name__=="__main__":
    main()
