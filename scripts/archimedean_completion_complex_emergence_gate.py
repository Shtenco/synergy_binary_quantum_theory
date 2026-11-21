#!/usr/bin/env python3
"""Exact controls for the arithmetic ladder Q -> R and Q[J] -> C.

The script does not prove Ostrowski's theorem (used as a standard theorem in
accompanying documentation). It supplies exact finite controls showing:
- the Archimedean completion can be approached by nested rational intervals;
- p-adic and Archimedean notions of closeness are genuinely different;
- the rational product formula balances the infinite and finite places;
- the real 2x2 operator J with J^2=-I implements ordinary complex arithmetic
  already over Q, so after Archimedean completion the same structure is C.

No floating-point arithmetic is used for theorem checks.
"""

import argparse
import json
from fractions import Fraction
from math import isqrt


def factor_integer(n):
    n = abs(n)
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def vp_fraction(x, p):
    if x == 0:
        raise ValueError("v_p(0) not finite")
    num = abs(x.numerator)
    den = x.denominator
    v = 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v


def p_adic_abs_fraction(x, p):
    v = vp_fraction(x, p)
    if v >= 0:
        return Fraction(1, p**v)
    return Fraction(p ** (-v), 1)


def check_sqrt2_archimedean_completion():
    rows = []
    ok = True
    prev_lo = None
    prev_hi = None
    for n in range(1, 25):
        D = 2**n
        # floor(sqrt(2)*D) = floor(sqrt(2*D^2)) exactly.
        a = isqrt(2 * D * D)
        lo = Fraction(a, D)
        hi = Fraction(a + 1, D)
        bracket = lo * lo <= 2 < hi * hi
        nested = True
        if prev_lo is not None:
            nested = lo >= prev_lo and hi <= prev_hi
        width = hi - lo
        width_ok = width == Fraction(1, D)
        ok = ok and bracket and nested and width_ok
        rows.append(
            {
                "n": n,
                "lo": f"{lo.numerator}/{lo.denominator}",
                "hi": f"{hi.numerator}/{hi.denominator}",
                "width": f"{width.numerator}/{width.denominator}",
                "bracket": bracket,
                "nested": nested,
            }
        )
        prev_lo, prev_hi = lo, hi
    return {
        "nested_rational_intervals_for_sqrt2": rows,
        "final_width": rows[-1]["width"],
        "pass": ok,
    }


def check_padic_vs_archimedean_closeness():
    result = {}
    ok = True
    for p in [2, 3, 5, 7]:
        rows = []
        x = 0
        for n in range(1, 9):
            old = x
            x += p ** (n - 1)
            diff = Fraction(x - old, 1)
            v = vp_fraction(diff, p)
            padic = p_adic_abs_fraction(diff, p)
            real_abs = abs(diff)
            expected_v = n - 1
            case_ok = v == expected_v and padic == Fraction(1, p ** (n - 1))
            ok = ok and case_ok
            rows.append(
                {
                    "n": n,
                    "x_n": x,
                    "increment": int(diff),
                    "v_p_increment": v,
                    "p_adic_abs_increment": f"{padic.numerator}/{padic.denominator}",
                    "archimedean_abs_increment": int(real_abs),
                    "pass": case_ok,
                }
            )
        # x_n = 1+p+...+p^(n-1) tends p-adically to 1/(1-p),
        # while it diverges in the ordinary absolute value.
        target = Fraction(1, 1 - p)
        err = Fraction(x, 1) - target
        target_check = vp_fraction(err, p) == 8
        ok = ok and target_check
        result[str(p)] = {
            "rows": rows,
            "p_adic_limit": f"{target.numerator}/{target.denominator}",
            "v_p_final_error": vp_fraction(err, p),
            "pass": all(r["pass"] for r in rows) and target_check,
        }
    return {"cases": result, "pass": ok}


def check_product_formula():
    cases = [
        Fraction(-84, 275),
        Fraction(17, 29),
        Fraction(2**7 * 3**2, 5**3 * 11),
        Fraction(-13 * 17, 2**4 * 7**2),
    ]
    rows = []
    ok = True
    for x in cases:
        primes = sorted(set(factor_integer(x.numerator)) | set(factor_integer(x.denominator)))
        finite_product = Fraction(1, 1)
        vals = {}
        for p in primes:
            ap = p_adic_abs_fraction(x, p)
            vals[str(p)] = f"{ap.numerator}/{ap.denominator}"
            finite_product *= ap
        arch = abs(x)
        total = arch * finite_product
        case_ok = total == 1
        ok = ok and case_ok
        rows.append(
            {
                "x": f"{x.numerator}/{x.denominator}",
                "archimedean_abs": f"{arch.numerator}/{arch.denominator}",
                "finite_p_adic_abs": vals,
                "finite_product": f"{finite_product.numerator}/{finite_product.denominator}",
                "total_product": f"{total.numerator}/{total.denominator}",
                "pass": case_ok,
            }
        )
    return {"cases": rows, "pass": ok}


def pair_add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def pair_mul(z, w):
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c)


def matrix_of_pair(z):
    a, b = z
    return ((a, -b), (b, a))


def mmul2(A, B):
    return tuple(
        tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def madd2(A, B):
    return tuple(tuple(A[i][j] + B[i][j] for j in range(2)) for i in range(2))


def check_rational_complex_structure():
    vals = [
        Fraction(-3, 2),
        Fraction(-1, 3),
        Fraction(0, 1),
        Fraction(2, 5),
        Fraction(7, 4),
    ]
    J = ((Fraction(0), Fraction(-1)), (Fraction(1), Fraction(0)))
    minusI = ((Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1)))
    ok = mmul2(J, J) == minusI
    checked = 0
    for a in vals:
        for b in vals:
            for c in vals:
                for d in vals:
                    z = (a, b)
                    w = (c, d)
                    if matrix_of_pair(pair_add(z, w)) != madd2(matrix_of_pair(z), matrix_of_pair(w)):
                        ok = False
                        break
                    if matrix_of_pair(pair_mul(z, w)) != mmul2(matrix_of_pair(z), matrix_of_pair(w)):
                        ok = False
                        break
                    checked += 1
                if not ok:
                    break
            if not ok:
                break
        if not ok:
            break
    return {
        "J_squared_eq_minus_I": mmul2(J, J) == minusI,
        "rational_pair_cases_checked": checked,
        "QJ_is_Qi_on_checked_cases": ok,
        "pass": ok,
    }


def check_one_common_scale_additivity():
    # Finite exhaustive control: any additive map f: Z -> Q with f(1)=s
    # must equal s*n. We verify the recurrence on a symmetric finite window.
    scales = [Fraction(1, 1), Fraction(3, 7), Fraction(-5, 2)]
    out = {}
    ok = True
    for s in scales:
        f = {0: Fraction(0), 1: s}
        for n in range(2, 21):
            f[n] = f[n - 1] + f[1]
        for n in range(1, 21):
            f[-n] = -f[n]
        case_ok = all(f[n] == s * n for n in range(-20, 21))
        # Exhaustively check additivity wherever sum remains in window.
        for a in range(-10, 11):
            for b in range(-10, 11):
                if f[a + b] != f[a] + f[b]:
                    case_ok = False
        ok = ok and case_ok
        out[f"{s.numerator}/{s.denominator}"] = case_ok
    return {"scales": out, "pass": ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path")
    args = ap.parse_args()

    report = {
        "archimedean_completion_sqrt2_control": check_sqrt2_archimedean_completion(),
        "padic_vs_archimedean_closeness": check_padic_vs_archimedean_closeness(),
        "rational_product_formula": check_product_formula(),
        "rational_complex_structure_QJ": check_rational_complex_structure(),
        "one_common_scale_additivity": check_one_common_scale_additivity(),
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
