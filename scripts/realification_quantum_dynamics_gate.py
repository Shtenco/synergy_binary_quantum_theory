#!/usr/bin/env python3
"""Exact realification gate for finite-dimensional complex quantum dynamics.

A complex matrix H=A+iB is represented by the real block
R(H)=[[A,-B],[B,A]]. For Hermitian H, A is symmetric and B skew-symmetric,
so R(H) is real symmetric. With J=[[0,-I],[I,0]], the Schrodinger equation
    i dpsi/dt = H psi
is exactly
    dv/dt = -J R(H) v
on ordinary real coordinates v=(Re psi, Im psi).
The real generator is skew-symmetric and preserves the Euclidean norm.

This is an exact representation theorem, not a derivation of the physical
Hamiltonian or time variable from the microscopic information graph.
"""

import argparse
import json
from fractions import Fraction


def zeros(n, m):
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def eye(n):
    A = zeros(n, n)
    for i in range(n):
        A[i][i] = Fraction(1)
    return A


def add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def neg(A):
    return [[-x for x in row] for row in A]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def eq(A, B):
    return A == B


def block2(A, B, C, D):
    n = len(A)
    m = len(A[0])
    out = []
    for i in range(n):
        out.append(A[i] + B[i])
    for i in range(len(C)):
        out.append(C[i] + D[i])
    return out


def realify(Z):
    A, B = Z
    return block2(A, neg(B), B, A)


def cmat_mul(Z, W):
    A, B = Z
    C, D = W
    return (add(mul(A, C), neg(mul(B, D))), add(mul(A, D), mul(B, C)))


def cmat_adjoint(Z):
    A, B = Z
    return (transpose(A), neg(transpose(B)))


def Jn(n):
    I = eye(n)
    O = zeros(n, n)
    return block2(O, neg(I), I, O)


def scale(A, s):
    return [[s * x for x in row] for row in A]


def is_symmetric(A):
    return A == transpose(A)


def is_skew(A):
    return A == neg(transpose(A))


def matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def dot(v, w):
    return sum(a * b for a, b in zip(v, w))


def complex_apply(Z, psi):
    A, B = Z
    x, y = psi
    # (A+iB)(x+iy)=(Ax-By)+i(Bx+Ay)
    return (add_vec(matvec(A, x), neg_vec(matvec(B, y))), add_vec(matvec(B, x), matvec(A, y)))


def add_vec(x, y):
    return [a + b for a, b in zip(x, y)]


def neg_vec(x):
    return [-a for a in x]


def flatten_complex_vec(psi):
    return list(psi[0]) + list(psi[1])


def schrodinger_rhs_complex(H, psi):
    # dpsi = -i H psi. If Hpsi=(u,v), -i(u+iv)=v-iu.
    u, v = complex_apply(H, psi)
    return (v, neg_vec(u))


def examples():
    F = Fraction
    return [
        (
            [[F(1), F(2)], [F(2), F(3)]],
            [[F(0), F(0)], [F(0), F(0)]],
        ),
        (
            [[F(2), F(1, 2)], [F(1, 2), F(-1)]],
            [[F(0), F(3, 4)], [F(-3, 4), F(0)]],
        ),
        (
            [[F(1, 3), F(-2, 5)], [F(-2, 5), F(7, 6)]],
            [[F(0), F(-5, 7)], [F(5, 7), F(0)]],
        ),
    ]


def arbitrary_complex_examples():
    F = Fraction
    return [
        (
            [[F(1), F(2)], [F(0), F(-1)]],
            [[F(3), F(0)], [F(1), F(2)]],
        ),
        (
            [[F(2, 3), F(-1, 2)], [F(4, 5), F(3, 7)]],
            [[F(0), F(1, 3)], [F(-2, 3), F(5, 4)]],
        ),
    ]


def check_algebra_homomorphism():
    Zs = arbitrary_complex_examples() + examples()
    checked = 0
    ok = True
    for Z in Zs:
        if realify(cmat_adjoint(Z)) != transpose(realify(Z)):
            ok = False
            break
        for W in Zs:
            if realify(cmat_mul(Z, W)) != mul(realify(Z), realify(W)):
                ok = False
                break
            checked += 1
        if not ok:
            break
    return {"matrix_pairs_checked": checked, "adjoint_to_transpose": ok, "pass": ok}


def check_hermitian_dynamics():
    rows = []
    ok = True
    test_psis = [
        ([Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]),
        ([Fraction(2, 3), Fraction(-1, 4)], [Fraction(5, 7), Fraction(3, 2)]),
    ]
    for idx, H in enumerate(examples()):
        A, B = H
        RH = realify(H)
        J = Jn(len(A))
        K = scale(mul(J, RH), Fraction(-1))
        hermitian = cmat_adjoint(H) == H
        symmetric = is_symmetric(RH)
        commutes = mul(J, RH) == mul(RH, J)
        skew_generator = is_skew(K)
        j2 = mul(J, J) == scale(eye(2 * len(A)), Fraction(-1))
        vector_ok = True
        norm_derivatives = []
        for psi in test_psis:
            v = flatten_complex_vec(psi)
            rhs_complex = flatten_complex_vec(schrodinger_rhs_complex(H, psi))
            rhs_real = matvec(K, v)
            if rhs_complex != rhs_real:
                vector_ok = False
            deriv = dot(v, rhs_real) * 2
            norm_derivatives.append(str(deriv))
            if deriv != 0:
                vector_ok = False
        case_ok = hermitian and symmetric and commutes and skew_generator and j2 and vector_ok
        ok = ok and case_ok
        rows.append(
            {
                "case": idx,
                "complex_H_hermitian": hermitian,
                "realification_symmetric": symmetric,
                "R_H_commutes_with_J": commutes,
                "J_squared_minus_I": j2,
                "real_generator_skew": skew_generator,
                "complex_and_real_rhs_equal": vector_ok,
                "norm_derivatives_2v_dot_Kv": norm_derivatives,
                "pass": case_ok,
            }
        )
    return {"cases": rows, "pass": ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path")
    args = ap.parse_args()

    report = {
        "realification_star_algebra": check_algebra_homomorphism(),
        "hermitian_schrodinger_real_flow": check_hermitian_dynamics(),
    }
    report["overall_pass"] = all(v["pass"] for v in report.values())

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    if not report["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
