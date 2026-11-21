#!/usr/bin/env python3
"""Exact gate for modular-complex -> ordinary arithmetic bridges.

No floating point and no external dependencies are used.

The gate checks five independent facts:
1. Gaussian modular arithmetic a+bi mod N embeds exactly in ordinary 2x2
   modular integer matrices.
2. For primes p == 1 mod 4, the Gaussian ring splits exactly into two
   ordinary scalar modular channels F_p x F_p.
3. For primes p == 3 mod 4, x^2+1 is irreducible and the nonzero Gaussian
   pairs form the quadratic field F_{p^2}; p=2 is the ramified exception.
4. The oriented C4 shift contains the real 2x2 quarter-turn block J with
   J^2=-I. Reversing the cycle sends J -> -J, i.e. complex conjugation.
5. Several coprime modular channels can be reconstructed by CRT into ordinary
   Gaussian integers whenever a centered bound prevents wrap ambiguity.

This is an arithmetic/representation theorem. It is not a claim that physical
real numbers or the Archimedean order are already derived from one finite ring.
"""

from itertools import product
from math import prod
import argparse
import json


def gadd(z, w, n):
    return ((z[0] + w[0]) % n, (z[1] + w[1]) % n)


def gmul(z, w, n):
    a, b = z
    c, d = w
    return ((a * c - b * d) % n, (a * d + b * c) % n)


def phi_matrix(z, n):
    a, b = z
    return ((a % n, (-b) % n), (b % n, a % n))


def madd(A, B, n):
    return tuple(
        tuple((A[i][j] + B[i][j]) % n for j in range(2)) for i in range(2)
    )


def mmul(A, B, n):
    return tuple(
        tuple(
            sum(A[i][k] * B[k][j] for k in range(2)) % n
            for j in range(2)
        )
        for i in range(2)
    )


def mmul_int(A, B):
    return tuple(
        tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def mtranspose(A):
    return ((A[0][0], A[1][0]), (A[0][1], A[1][1]))


def det2(A, n):
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % n


def trace2(A, n):
    return (A[0][0] + A[1][1]) % n


def roots_minus_one(p):
    return [r for r in range(p) if (r * r + 1) % p == 0]


def split_pair(z, p, r):
    a, b = z
    return ((a + r * b) % p, (a - r * b) % p)


def unsplit_pair(uv, p, r):
    u, v = uv
    inv2 = pow(2, -1, p)
    inv2r = pow((2 * r) % p, -1, p)
    return (((u + v) * inv2) % p, ((u - v) * inv2r) % p)


def matvec(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M)))


def crt(residues, moduli):
    M = prod(moduli)
    x = 0
    for a, m in zip(residues, moduli):
        Mi = M // m
        x = (x + a * Mi * pow(Mi, -1, m)) % M
    return x


def centered(x, M):
    x %= M
    return x if x <= M // 2 else x - M


def check_matrix_decomplexification():
    checked_moduli = [2, 3, 5, 7, 8, 13]
    result = {}
    for n in checked_moduli:
        ok = True
        for a, b, c, d in product(range(n), repeat=4):
            z = (a, b)
            w = (c, d)
            Pz = phi_matrix(z, n)
            Pw = phi_matrix(w, n)
            if phi_matrix(gadd(z, w, n), n) != madd(Pz, Pw, n):
                ok = False
                break
            if phi_matrix(gmul(z, w, n), n) != mmul(Pz, Pw, n):
                ok = False
                break
            if phi_matrix((a, (-b) % n), n) != mtranspose(Pz):
                ok = False
                break
            if det2(Pz, n) != (a * a + b * b) % n:
                ok = False
                break
            if trace2(Pz, n) != (2 * a) % n:
                ok = False
                break
        result[str(n)] = ok
    return result


def check_split_primes():
    result = {}
    for p in [5, 13, 17]:
        roots = roots_minus_one(p)
        r = roots[0]
        ok = len(roots) == 2
        seen = set()
        for a, b in product(range(p), repeat=2):
            z = (a, b)
            uv = split_pair(z, p, r)
            seen.add(uv)
            if unsplit_pair(uv, p, r) != z:
                ok = False
                break
        if ok:
            for a, b, c, d in product(range(p), repeat=4):
                z = (a, b)
                w = (c, d)
                sz = split_pair(z, p, r)
                sw = split_pair(w, p, r)
                if split_pair(gadd(z, w, p), p, r) != (
                    (sz[0] + sw[0]) % p,
                    (sz[1] + sw[1]) % p,
                ):
                    ok = False
                    break
                if split_pair(gmul(z, w, p), p, r) != (
                    (sz[0] * sw[0]) % p,
                    (sz[1] * sw[1]) % p,
                ):
                    ok = False
                    break
        result[str(p)] = {
            "root_r": r,
            "roots_of_minus_one": roots,
            "bijective": len(seen) == p * p,
            "pass": ok and len(seen) == p * p,
        }
    return result


def check_nonsplit_primes():
    result = {}
    for p in [3, 7, 11, 19]:
        ok = roots_minus_one(p) == []
        for a, b in product(range(p), repeat=2):
            if a == 0 and b == 0:
                continue
            norm = (a * a + b * b) % p
            if norm == 0:
                ok = False
                break
            invnorm = pow(norm, -1, p)
            inv = ((a * invnorm) % p, (-b * invnorm) % p)
            if gmul((a, b), inv, p) != (1, 0):
                ok = False
                break
        result[str(p)] = ok
    return result


def check_ramified_two():
    eps = (1, 1)  # 1+i mod 2
    return {
        "epsilon_nonzero": eps != (0, 0),
        "epsilon_squared": gmul(eps, eps, 2),
        "pass": eps != (0, 0) and gmul(eps, eps, 2) == (0, 0),
    }


def check_c4_complex_block():
    # S|k> = |k+1 mod 4> in the ordinary binary/integer basis.
    S = (
        (0, 0, 0, 1),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
    )
    u0 = (1, 1, 1, 1)
    u2 = (1, -1, 1, -1)
    uc = (1, 0, -1, 0)
    us = (0, 1, 0, -1)
    J = ((0, -1), (1, 0))
    minus_I = ((-1, 0), (0, -1))
    checks = {
        "S_u0_eq_u0": matvec(S, u0) == u0,
        "S_u2_eq_minus_u2": matvec(S, u2) == tuple(-x for x in u2),
        "S_uc_eq_us": matvec(S, uc) == us,
        "S_us_eq_minus_uc": matvec(S, us) == tuple(-x for x in uc),
        "J_squared_eq_minus_I": mmul_int(J, J) == minus_I,
        "orientation_reversal_J_to_minus_J": ((0, 1), (-1, 0))
        == tuple(tuple(-J[i][j] for j in range(2)) for i in range(2)),
    }
    checks["pass"] = all(checks.values())
    return checks


def check_mu4_to_mod5():
    # 2^2 = -1 mod 5, so powers of 2 model 1,i,-1,-i multiplicatively.
    images = [pow(2, k, 5) for k in range(4)]
    hom = all(
        pow(2, (a + b) % 4, 5) == (pow(2, a, 5) * pow(2, b, 5)) % 5
        for a, b in product(range(4), repeat=2)
    )
    return {
        "images_1_i_minus1_minusi": images,
        "expected": [1, 2, 4, 3],
        "pass": images == [1, 2, 4, 3] and hom,
    }


def check_crt_lift():
    moduli = [5, 7, 11]
    M = prod(moduli)
    z = (7, 11)
    w = (-5, 6)
    targets = {
        "sum": (z[0] + w[0], z[1] + w[1]),
        "product": (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]),
    }
    recovered = {}
    ok = True
    for name, target in targets.items():
        re = crt([target[0] % m for m in moduli], moduli)
        im = crt([target[1] % m for m in moduli], moduli)
        rec = (centered(re, M), centered(im, M))
        recovered[name] = {"target": target, "recovered": rec, "pass": rec == target}
        ok = ok and rec == target
    return {
        "moduli": moduli,
        "product_modulus": M,
        "centered_unique_window": [-M // 2, M // 2],
        "cases": recovered,
        "pass": ok,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path", default=None)
    args = ap.parse_args()

    report = {
        "matrix_decomplexification": check_matrix_decomplexification(),
        "split_primes_p_eq_1_mod_4": check_split_primes(),
        "nonsplit_primes_p_eq_3_mod_4": check_nonsplit_primes(),
        "ramified_p_2": check_ramified_two(),
        "q2_oriented_C4_complex_block": check_c4_complex_block(),
        "mu4_to_mod5": check_mu4_to_mod5(),
        "crt_bounded_integer_lift": check_crt_lift(),
    }

    def all_true(obj):
        if isinstance(obj, bool):
            return obj
        if isinstance(obj, dict):
            if "pass" in obj and isinstance(obj["pass"], bool):
                return obj["pass"]
            return all(all_true(v) for v in obj.values())
        return True

    report["overall_pass"] = all(
        [
            all(report["matrix_decomplexification"].values()),
            all(v["pass"] for v in report["split_primes_p_eq_1_mod_4"].values()),
            all(report["nonsplit_primes_p_eq_3_mod_4"].values()),
            report["ramified_p_2"]["pass"],
            report["q2_oriented_C4_complex_block"]["pass"],
            report["mu4_to_mod5"]["pass"],
            report["crt_bounded_integer_lift"]["pass"],
        ]
    )

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    if not report["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
