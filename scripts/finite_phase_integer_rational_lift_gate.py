#!/usr/bin/env python3
"""Exact finite-phase -> ordinary integer/rational lift controls.

Checks:
- primitive M-th roots in F_p when M divides p-1;
- exact modular Fourier diagonalization of an integer cyclic shift;
- modular phase addition + winding reconstructs ordinary integer addition;
- rational reconstruction from one sufficiently large modulus and bounds;
- nested dyadic phase groups provide a controlled dense U(1) target.

The final conversion to a dimensionful physical observable still requires one
common scale factor; this script does not fit or derive that factor.
"""

import argparse
import json
from itertools import product
from math import gcd, pi


def divisors(n):
    out = []
    for d in range(1, n + 1):
        if n % d == 0:
            out.append(d)
    return out


def multiplicative_order(a, p):
    if a % p == 0:
        return 0
    x = 1
    for k in range(1, p):
        x = (x * a) % p
        if x == 1:
            return k
    raise RuntimeError("order not found")


def primitive_root_of_order(M, p):
    for a in range(2, p):
        if multiplicative_order(a, p) == M:
            return a
    return None


def shift_apply(v):
    # S|j> = |j+1>; coefficient at output index l is old l-1.
    return (v[-1],) + tuple(v[:-1])


def vec_scale(a, v, p):
    return tuple((a * x) % p for x in v)


def rank_mod_p(A, p):
    A = [list(row) for row in A]
    m = len(A)
    n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if A[i][c] % p != 0), None)
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        inv = pow(A[r][c] % p, -1, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(m):
            if i == r:
                continue
            f = A[i][c] % p
            if f:
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def check_modular_fourier():
    cases = [(4, 5), (8, 17), (12, 13), (16, 17)]
    result = {}
    for M, p in cases:
        root = primitive_root_of_order(M, p)
        ok = root is not None and (p - 1) % M == 0
        columns = []
        if ok:
            for k in range(M):
                # v_k[j] = root^(-j*k), so S v_k = root^k v_k.
                v = tuple(pow(root, (-j * k) % M, p) for j in range(M))
                columns.append(v)
                if shift_apply(v) != vec_scale(pow(root, k, p), v, p):
                    ok = False
                    break
        if ok:
            F = [[columns[c][r] for c in range(M)] for r in range(M)]
            ok = rank_mod_p(F, p) == M
        result[f"M{M}_p{p}"] = {
            "root": root,
            "order": multiplicative_order(root, p) if root else None,
            "fourier_rank": M if ok else None,
            "pass": ok,
        }
    return result


def check_root_existence_criterion():
    # Exhaustively verify small primes/M that an element of exact order M exists
    # exactly when M divides p-1 (for M>1).
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    mismatches = []
    checked = 0
    for p in primes:
        for M in range(2, min(16, p)):
            found = primitive_root_of_order(M, p) is not None
            expected = (p - 1) % M == 0
            checked += 1
            if found != expected:
                mismatches.append([p, M, found, expected])
    return {"checked": checked, "mismatches": mismatches, "pass": not mismatches}


def check_winding_lift():
    results = {}
    for M in [4, 5, 8, 16]:
        pair_ok = True
        for a, b in product(range(M), repeat=2):
            total = a + b
            q, r = divmod(total, M)
            if total != r + M * q:
                pair_ok = False
                break
        seq = [M - 1, M - 2, 3, 1, M - 1]
        ordinary = sum(seq)
        winding, residue = divmod(ordinary, M)
        results[str(M)] = {
            "pair_exhaustive": pair_ok,
            "sequence": seq,
            "ordinary_sum": ordinary,
            "residue": residue,
            "winding": winding,
            "reconstructed": residue + M * winding,
            "pass": pair_ok and ordinary == residue + M * winding,
        }
    return results


def centered(x, M):
    x %= M
    return x if x <= M // 2 else x - M


def rational_reconstruct_bruteforce(x, M, A, B):
    candidates = []
    for b in range(1, B + 1):
        if gcd(b, M) != 1:
            continue
        a = centered((x * b) % M, M)
        if abs(a) <= A and gcd(abs(a), b) == 1:
            candidates.append((a, b))
    return sorted(set(candidates))


def check_rational_reconstruction():
    M = 5 * 7 * 11 * 13  # 5005
    A = 20
    B = 30
    cases = [(17, 29), (-19, 23)]
    result = {}
    ok = 2 * A * B < M
    for a, b in cases:
        x = (a * pow(b, -1, M)) % M
        candidates = rational_reconstruct_bruteforce(x, M, A, B)
        good = candidates == [(a, b)]
        result[f"{a}/{b}"] = {
            "residue": x,
            "candidates": candidates,
            "unique": good,
        }
        ok = ok and good
    return {
        "M": M,
        "A": A,
        "B": B,
        "uniqueness_bound_2AB_lt_M": 2 * A * B < M,
        "cases": result,
        "pass": ok,
    }


def check_dyadic_phase_refinement():
    # Purely arithmetic nesting: mu_(2^g) embeds in mu_(2^(g+1)) by exponent doubling.
    nesting_ok = True
    for g in range(2, 10):
        M = 2**g
        Mp = 2 ** (g + 1)
        # exponent k in mu_M maps to 2k in mu_2M; addition is respected.
        for a, b in product(range(M), repeat=2):
            lhs = (2 * ((a + b) % M)) % Mp
            rhs = ((2 * a) + (2 * b)) % Mp
            if lhs != rhs:
                nesting_ok = False
                break
        if not nesting_ok:
            break

    # Deterministic angular approximation control for several irrational-looking targets.
    targets = [0.123456789, 1.0, pi / 7, 2.345678901]
    errors = {}
    monotone_ok = True
    for theta in targets:
        seq = []
        prev_bound = None
        for g in range(2, 15):
            M = 2**g
            k = round(M * theta / (2 * pi))
            approx = 2 * pi * k / M
            err = abs(theta - approx)
            bound = pi / M + 1e-15
            if err > bound:
                monotone_ok = False
            seq.append({"M": M, "error": err, "bound_pi_over_M": pi / M})
            prev_bound = bound
        errors[str(theta)] = seq
    return {"nesting_pass": nesting_ok, "density_bound_pass": monotone_ok, "errors": errors, "pass": nesting_ok and monotone_ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_path")
    args = ap.parse_args()

    report = {
        "root_existence_criterion": check_root_existence_criterion(),
        "modular_fourier_diagonalization": check_modular_fourier(),
        "phase_plus_winding_to_integer": check_winding_lift(),
        "rational_reconstruction": check_rational_reconstruction(),
        "dyadic_phase_refinement": check_dyadic_phase_refinement(),
    }
    report["overall_pass"] = all(
        [
            report["root_existence_criterion"]["pass"],
            all(v["pass"] for v in report["modular_fourier_diagonalization"].values()),
            all(v["pass"] for v in report["phase_plus_winding_to_integer"].values()),
            report["rational_reconstruction"]["pass"],
            report["dyadic_phase_refinement"]["pass"],
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
