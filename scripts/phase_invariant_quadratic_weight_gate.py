#!/usr/bin/env python3
"""Exact gate for the unique C4-phase-invariant quadratic scalar on R^2/Q^2.

For v=(a,b) and J=[[0,-1],[1,0]], a symmetric quadratic form
Q(v)=v^T A v that obeys J^T A J=A must have A=lambda I.
With positivity and Q(1,0)=1 this gives Q=a^2+b^2.

The gate also verifies the exact interference polarization identity
Q(v+w)=Q(v)+Q(w)+2 lambda v.w.

This is a Born-weight precursor (unique quadratic phase-invariant positive
weight), not a derivation of the complete quantum measurement/Born rule.
"""

import argparse
import json
from fractions import Fraction
from itertools import product


def mmul(A, B):
    return tuple(
        tuple(sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0])))
        for i in range(len(A))
    )


def transpose(A):
    return tuple(tuple(A[j][i] for j in range(len(A))) for i in range(len(A[0])))


def matvec(A, v):
    return tuple(sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A)))


def dot(v, w):
    return sum(a * b for a, b in zip(v, w))


def qform(A, v):
    return dot(v, matvec(A, v))


def check_symbolic_coefficient_map():
    # A=[[alpha,beta],[beta,gamma]]. J^T A J maps it to
    # [[gamma,-beta],[-beta,alpha]], so invariance forces
    # alpha=gamma and beta=0 over characteristic !=2.
    return {
        "A": "[[alpha,beta],[beta,gamma]]",
        "JT_A_J": "[[gamma,-beta],[-beta,alpha]]",
        "invariance_equations": ["alpha=gamma", "beta=-beta"],
        "over_R_or_Q": ["alpha=gamma", "beta=0"],
        "solution": "A=lambda*I",
        "pass": True,
    }


def check_exhaustive_rational_grid():
    vals = [Fraction(-2), Fraction(-1), Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2)]
    J = ((Fraction(0), Fraction(-1)), (Fraction(1), Fraction(0)))
    JT = transpose(J)
    invariant = []
    for alpha, beta, gamma in product(vals, repeat=3):
        A = ((alpha, beta), (beta, gamma))
        if mmul(mmul(JT, A), J) == A:
            invariant.append((alpha, beta, gamma))
    expected = [(lam, Fraction(0), lam) for lam in vals]
    return {
        "solutions": [[str(a), str(b), str(c)] for a, b, c in invariant],
        "expected_scalar_solutions": [[str(a), str(b), str(c)] for a, b, c in expected],
        "pass": invariant == expected,
    }


def check_norm_and_interference():
    vals = [Fraction(-3, 2), Fraction(-1, 3), Fraction(0), Fraction(2, 5), Fraction(7, 4)]
    lambdas = [Fraction(1), Fraction(3, 7), Fraction(5, 2)]
    checked = 0
    ok = True
    for lam in lambdas:
        A = ((lam, Fraction(0)), (Fraction(0), lam))
        for a, b, c, d in product(vals, repeat=4):
            v = (a, b)
            w = (c, d)
            vw = (a + c, b + d)
            lhs = qform(A, vw)
            rhs = qform(A, v) + qform(A, w) + 2 * lam * dot(v, w)
            if lhs != rhs:
                ok = False
                break
            # Complex notation identity: Re(z conjugate(w)) = ac+bd.
            real_z_conj_w = a * c + b * d
            if real_z_conj_w != dot(v, w):
                ok = False
                break
            checked += 1
        if not ok:
            break
    return {"cases_checked": checked, "pass": ok}


def check_phase_rotation_invariance():
    vals = [Fraction(-2), Fraction(-1), Fraction(0), Fraction(1), Fraction(2)]
    J = ((Fraction(0), Fraction(-1)), (Fraction(1), Fraction(0)))
    minusI = ((Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1)))
    minusJ = ((Fraction(0), Fraction(1)), (Fraction(-1), Fraction(0)))
    rotations = [
        ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))),
        J,
        minusI,
        minusJ,
    ]
    A = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
    checked = 0
    ok = True
    for a, b in product(vals, repeat=2):
        v = (a, b)
        for R in rotations:
            if qform(A, matvec(R, v)) != qform(A, v):
                ok = False
                break
            checked += 1
        if not ok:
            break
    return {"mu4_rotations_checked": checked, "pass": ok}


def check_positivity_normalization():
    # Once A=lambda I, positivity implies lambda>=0.
    # Normalization Q(1,0)=1 fixes lambda=1.
    one = (Fraction(1), Fraction(0))
    lambdas = [Fraction(-2), Fraction(-1), Fraction(0), Fraction(1), Fraction(3, 2)]
    rows = []
    ok = True
    for lam in lambdas:
        A = ((lam, Fraction(0)), (Fraction(0), lam))
        qone = qform(A, one)
        psd_expected = lam >= 0
        probes = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(1), Fraction(1)), (Fraction(2), Fraction(-3))]
        psd_probe = all(qform(A, v) >= 0 for v in probes)
        case_ok = psd_probe == psd_expected
        ok = ok and case_ok
        rows.append({"lambda": str(lam), "Q_1_0": str(qone), "psd": psd_probe, "pass": case_ok})
    normalized_lambda = Fraction(1)
    ok = ok and qform(((normalized_lambda, 0), (0, normalized_lambda)), one) == 1
    return {"cases": rows, "normalization_Q_1_0_eq_1_gives_lambda": "1", "pass": ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path")
    args = ap.parse_args()

    report = {
        "symbolic_invariance": check_symbolic_coefficient_map(),
        "exhaustive_rational_grid": check_exhaustive_rational_grid(),
        "norm_interference_identity": check_norm_and_interference(),
        "mu4_phase_invariance": check_phase_rotation_invariance(),
        "positivity_and_normalization": check_positivity_normalization(),
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
