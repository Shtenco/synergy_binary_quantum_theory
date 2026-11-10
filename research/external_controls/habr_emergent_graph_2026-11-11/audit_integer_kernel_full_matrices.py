#!/usr/bin/env python3
"""Exact saturated integer-kernel + LLL control for the full 69x6 and 70x6 matrices.

Frozen from homoastricus/emergent_graph_theory:
- physics_constants/constant_matrix_analysis.py (79b5e570ab990c6143ace748881ba20203de049c)
- matrix_basis/matrix_analysis.py (e58dd577e804b07d07e6c0df1a01041639ad9100)

Requires SymPy with smith_normal_decomp support.
"""
from __future__ import annotations
import json
from functools import reduce
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp
from sympy.polys.domains import ZZ

NAMES69 = ["hbar","h","c","lP","tP","EP","G","mP","TP","kB","alpha","me","mp","eps0","mu0","qe","Rinf","a0","Z0","Phi0","lambda_e","lambda_p","vH","muon","tau","pion_pm","pion0","Lambda","kappa","mH_mW","m_qu_u","m_qu_d","m_qu_s","m_qu_c","m_qu_b","m_qu_t","m_proton_to_m_electron","m_tau_m_electron","m_Λ_barion","m_z_bozon","m_w_bozon","m_D0","m_J_ψ","m_eta","r_electron","D0_lifetime","Λ_b_lifetime","B_plus_lifetime","kaon_lifetime","mu_lifetime","tau_lifetime","pion_lifetime","D_plus_lifetime","m_Higgs","m_neitrino","Sigma_plus","Ksi_0","Omega_minus","Ksi_plus","Omega0_c","lambda_B0","Ksi_minus","Sigma_minus","phi_meson","omega_meson","eta_shtrih","rho_meson","neutron_lifetime","h2_connection_energy"]
M69 = sp.Matrix([
    [-2, -2, 0, 0, 3, -1],
    [0, -2, 1, 0, 3, -1],
    [-4, -4, 1, -1, 4, 0],
    [4, 0, 0, 1, 2, -1],
    [8, 4, -1, 2, -2, -1],
    [-10, -6, 1, -2, 5, 0],
    [-2, -10, 3, -1, 13, -1],
    [-2, 2, -1, 0, -3, 0],
    [6, 0, 1, 0, -4, 1],
    [-6, 0, -2, 0, 8, -1],
    [2, 0, -1, 2, -1, 0],
    [3, -1, 1, 0, 4, -1],
    [-3, -3, sp.Rational(1, 2), 0, 6, -1],
    [-6, 0, -3, -1, -20, 1],
    [14, 8, 1, 3, 12, -1],
    [-3, -3, -1, 0, -7, 0],
    [1, -3, -1, 3, 3, 0],
    [-3, 3, -1, -1, -4, 0],
    [10, 4, 2, 2, 16, -1],
    [1, 1, 2, 0, 10, -1],
    [1, 3, -1, 1, -5, 0],
    [7, 5, -sp.Rational(1, 2), 1, -7, 0],
    [5, 0, sp.Rational(3, 2), 0, 6, -1],
    [2, -3, 2, 0, 5, -1],
    [4, 4, sp.Rational(1, 2), 0, 5, -1],
    [-5, 0, -2, 0, 6, -1],
    [8, 6, 1, 0, 4, -1],
    [0, 0, -sp.Rational(1, 2), 0, 12, -2],
    [20, 6, 0, 3, -3, -1],
    [3, 1, -1, 0, 0, 0],
    [-4, 1, -2, 0, 5, -1],
    [-2, -3, 0, 0, 5, -1],
    [2, 2, -sp.Rational(1, 2), 0, 5, -1],
    [-4, -6, 2, 0, 6, -1],
    [-2, -3, 1, 0, 6, -1],
    [6, 6, -2, 0, 6, -1],
    [-6, -2, -sp.Rational(1, 2), 0, 2, 0],
    [1, 5, -sp.Rational(1, 2), 0, 1, 0],
    [1, 0, -2, 0, 6, -1],
    [2, -2, sp.Rational(5, 2), 0, 6, -1],
    [0, -2, 3, 0, 6, -1],
    [-1, -3, sp.Rational(1, 2), 0, 6, -1],
    [7, 0, 2, 0, 5, -1],
    [2, 0, 2, 0, 5, -1],
    [1, 3, -3, 3, -6, 0],
    [-6, -4, -2, 1, -4, 0],
    [1, 0, 0, 1, -5, 0],
    [-2, 0, 1, 1, -5, 0],
    [1, -3, 0, 0, -3, 0],
    [-2, -3, 0, 1, -2, 0],
    [-2, 0, 0, 0, -5, 0],
    [5, 4, 1, 0, -4, 0],
    [-5, -5, -sp.Rational(1, 2), 0, -4, 0],
    [3, -1, 2, 0, 6, -1],
    [1, 0, 0, -1, 2, -1],
    [-2, 2, -2, 0, 6, -1],
    [-2, -3, 0.5, 0, 6, -1],
    [-3, -3, 1, 0, 6, -1],
    [0, 0, -1, 0, 6, -1],
    [5, 1, 2, 0, 5, -1],
    [-1, -1, 0.5, 0, 6, -1],
    [-2, -3, 0.5, 0, 6, -1],
    [-2, 2, -2, 0, 6, -1],
    [4, 3, 0.5, 0, 5, -1],
    [3, 0, 2, 0, 5, -1],
    [4, 6, -1, 0, 5, -1],
    [0, 1, 2.5, 0, 5, -1],
    [0, 0, 1, 0, 1, 0],
    [-3, -9, 1, 2, 10, -1],
])
NAMES70 = ["hbar","h","c","lP","tP","EP","G","mP","TP","kB","alpha","me","mp","eps0","mu0","qe","Rinf","a0","Z0","Phi0","lambda_e","lambda_p","vH","muon","tau","pion_pm","pion0","Lambda","kappa","mH_mW","m_qu_u","m_qu_d","m_qu_s","m_qu_c","m_qu_b","m_qu_t","m_proton_to_m_electron","m_tau_m_electron","m_Λ_barion","m_z_bozon","m_w_bozon","m_D0","m_J_ψ","m_eta","r_electron","D0_lifetime","Λ_b_lifetime","B_plus_lifetime","kaon_lifetime","mu_lifetime","tau_lifetime","pion_lifetime","D_plus_lifetime","m_Higgs","m_neitrino","Sigma_plus","Ksi_0","Omega_minus","Ksi_plus","Omega0_c","lambda_B0","Ksi_minus","Sigma_minus","phi_meson","omega_meson","eta_shtrih","rho_meson","K_star_meson","neutron_lifetime","h2_connection_energy"]
M70 = sp.Matrix([
    [-2, -2, 0, 0, 3, -1],
    [0, -2, 1, 0, 3, -1],
    [-4, -4, 1, -1, 4, 0],
    [4, 0, 0, 1, 2, -1],
    [8, 4, -1, 2, -2, -1],
    [-10, -6, 1, -2, 5, 0],
    [-2, -10, 3, -1, 13, -1],
    [-2, 2, -1, 0, -3, 0],
    [6, 0, 1, 0, -4, 1],
    [-6, 0, -2, 0, 8, -1],
    [2, 0, -1, 2, -1, 0],
    [3, -1, 1, 0, 4, -1],
    [-3, -3, sp.Rational(1, 2), 0, 6, -1],
    [-6, 0, -3, -1, -20, 1],
    [14, 8, 1, 3, 12, -1],
    [-3, -3, -1, 0, -7, 0],
    [1, -3, -1, 3, 3, 0],
    [-3, 3, -1, -1, -4, 0],
    [10, 4, 2, 2, 16, -1],
    [1, 1, 2, 0, 10, -1],
    [1, 3, -1, 1, -5, 0],
    [7, 5, -sp.Rational(1, 2), 1, -7, 0],
    [5, 0, sp.Rational(3, 2), 0, 6, -1],
    [2, -3, 2, 0, 5, -1],
    [4, 4, sp.Rational(1, 2), 0, 5, -1],
    [-5, 0, -2, 0, 6, -1],
    [8, 6, 1, 0, 4, -1],
    [0, 0, -sp.Rational(1, 2), 0, 12, -2],
    [20, 6, 0, 3, -3, -1],
    [3, 1, -1, 0, 0, 0],
    [-4, 1, -2, 0, 5, -1],
    [-2, -3, 0, 0, 5, -1],
    [2, 2, -sp.Rational(1, 2), 0, 5, -1],
    [-4, -6, 2, 0, 6, -1],
    [-2, -3, 1, 0, 6, -1],
    [6, 6, -2, 0, 6, -1],
    [-6, -2, -sp.Rational(1, 2), 0, 2, 0],
    [1, 5, -sp.Rational(1, 2), 0, 1, 0],
    [1, 0, -2, 0, 6, -1],
    [2, -2, sp.Rational(5, 2), 0, 6, -1],
    [0, -2, 3, 0, 6, -1],
    [-1, -3, sp.Rational(1, 2), 0, 6, -1],
    [7, 0, 2, 0, 5, -1],
    [2, 0, 2, 0, 5, -1],
    [1, 3, -3, 3, -6, 0],
    [-6, -4, -2, 1, -4, 0],
    [1, 0, 0, 1, -5, 0],
    [-2, 0, 1, 1, -5, 0],
    [1, -3, 0, 0, -3, 0],
    [-2, -3, 0, 1, -2, 0],
    [-2, 0, 0, 0, -5, 0],
    [5, 4, 1, 0, -4, 0],
    [-5, -5, -sp.Rational(1, 2), 0, -4, 0],
    [3, -1, 2, 0, 6, -1],
    [1, 0, 0, -1, 2, -1],
    [-2, 2, -2, 0, 6, -1],
    [-2, -3, 0.5, 0, 6, -1],
    [-3, -3, 1, 0, 6, -1],
    [0, 0, -1, 0, 6, -1],
    [5, 1, 2, 0, 5, -1],
    [-1, -1, 0.5, 0, 6, -1],
    [-2, -3, 0.5, 0, 6, -1],
    [-2, 2, -2, 0, 6, -1],
    [4, 3, 0.5, 0, 5, -1],
    [3, 0, 2, 0, 5, -1],
    [4, 6, -1, 0, 5, -1],
    [0, 1, 2.5, 0, 5, -1],
    [2, 0, -2.5, 0, 6, -1],
    [0, 0, 1, 0, 1, 0],
    [-3, -9, 1, 2, 10, -1],
])

def primitive(v):
    den = sp.ilcm(*[sp.denom(x) for x in v])
    vals = [int(x*den) for x in v]
    g = reduce(gcd, (abs(x) for x in vals if x), 0)
    return [x//g for x in vals]

def audit(label, M, names):
    assert M.cols == 6 and M.rows == len(names)
    assert M.rank() == 6
    A = (2*M.T).applyfunc(int)
    D, U, V = smith_normal_decomp(A, domain=ZZ)
    diag = [int(D[i,i]) for i in range(6)]
    assert abs(int(V.det())) == 1

    saturated = V[:,6:].T
    assert saturated.shape == (M.rows-6, M.rows)
    assert A*saturated.T == sp.zeros(6, M.rows-6)
    reduced = saturated.lll(delta=sp.Rational(3,4))
    assert A*reduced.T == sp.zeros(6, M.rows-6)
    sat_rows = [[int(x) for x in row] for row in reduced.tolist()]
    for i, row in enumerate(sat_rows):
        assert A*sp.Matrix(row) == sp.zeros(6,1), (label, "saturated row", i)
    sat_max = max(abs(x) for row in sat_rows for x in row)

    naive = sp.Matrix([primitive(v) for v in M.T.nullspace()])
    assert naive.shape == saturated.shape
    assert A*naive.T == sp.zeros(6, M.rows-6)
    _, pivots = saturated.rref()
    cols = list(pivots)
    C = naive[:,cols] * saturated[:,cols].inv()
    assert all(x.q == 1 for x in C)
    assert C*saturated == naive
    index = abs(int(C.det()))

    naive_reduced = naive.lll(delta=sp.Rational(3,4))
    naive_rows = [[int(x) for x in row] for row in naive_reduced.tolist()]
    for i, row in enumerate(naive_rows):
        assert A*sp.Matrix(row) == sp.zeros(6,1), (label, "naive LLL row", i)
    naive_max = max(abs(x) for row in naive_rows for x in row)

    payload = {
        "source_matrix": label,
        "shape": list(M.shape),
        "rank": int(M.rank()),
        "nullity": M.rows-6,
        "smith_diagonal_of_2Mt": diag,
        "unimodular_transform_det_abs": abs(int(V.det())),
        "saturated_kernel_index": 1,
        "saturated_lll_max_abs_coefficient": sat_max,
        "saturated_lll_vectors": sat_rows,
        "all_saturated_rows_directly_verified": True,
        "rowwise_primitive_basis_index_in_saturated_kernel": str(index),
        "rowwise_index_power_of_two": index.bit_length()-1 if index & (index-1)==0 else None,
        "rowwise_lll_max_abs_coefficient": naive_max,
        "rowwise_lll_vectors": naive_rows,
        "all_rowwise_rows_directly_verified": True,
        "names": names,
    }
    with open(f"integer_kernel_{label.lower()}_audit.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(label, "shape", M.shape, "rank", M.rank(), "nullity", M.rows-6)
    print(label, "SNF(2*M.T)", diag)
    print(label, "saturated LLL maxabs", sat_max,
          "vectors", len(sat_rows), "direct_product_checks", len(sat_rows), "PASS")
    print(label, "rowwise basis index", index,
          "power_of_two", payload["rowwise_index_power_of_two"])
    print(label, "rowwise LLL maxabs", naive_max,
          "vectors", len(naive_rows), "direct_product_checks", len(naive_rows), "PASS")
    for i, row in enumerate(sat_rows, 1):
        norm2 = sum(x*x for x in row)
        support = [(names[j], x) for j,x in enumerate(row) if x]
        print(label, "SAT", i, "norm2", norm2, support)

def main():
    audit("M69", M69, NAMES69)
    audit("M70", M70, NAMES70)

if __name__ == "__main__":
    main()
