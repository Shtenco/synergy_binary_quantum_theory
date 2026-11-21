#!/usr/bin/env python3
"""Independent dense regression for the parity block-Lanczos B3 formula.

Construct a finite Hermitian block-Jacobi operator with known positive hopping
blocks B1, B2, B3 and zero diagonal blocks. Recover B3^dag B3 using only the
seed moments H2, H4 and H6, exactly as in the Peter-Weyl B3 gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

TOL = 5e-11


def herm(M):
    return (M + M.conj().T) / 2


def positive_block(rng, m, shift):
    X = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
    M = X.conj().T @ X + shift * np.eye(m)
    w, U = np.linalg.eigh(herm(M))
    return herm((U * np.sqrt(w)) @ U.conj().T)


def sqrt_psd(M):
    w, U = np.linalg.eigh(herm(M))
    if np.min(w) <= 0:
        raise RuntimeError("expected strictly positive test block")
    return herm((U * np.sqrt(w)) @ U.conj().T)


def run(seed=260906, m=4):
    rng = np.random.default_rng(seed)
    B1_true = positive_block(rng, m, 2.0)
    B2_true = positive_block(rng, m, 2.5)
    B3_true = positive_block(rng, m, 3.0)

    H = np.zeros((4 * m, 4 * m), dtype=complex)
    blocks = (B1_true, B2_true, B3_true)
    for n, B in enumerate(blocks):
        a = slice(n * m, (n + 1) * m)
        b = slice((n + 1) * m, (n + 2) * m)
        H[b, a] = B
        H[a, b] = B.conj().T
    H = herm(H)

    Q0 = np.zeros((4 * m, m), dtype=complex)
    Q0[:m, :] = np.eye(m)
    H2 = H @ H
    H3Q0 = H @ H @ H @ Q0
    K = herm(Q0.conj().T @ H2 @ Q0)
    H4 = herm(Q0.conj().T @ H2 @ H2 @ Q0)
    H6 = herm(H3Q0.conj().T @ H3Q0)

    B1 = sqrt_psd(K)
    B1inv = np.linalg.solve(B1, np.eye(m))
    Lambda = herm(B1inv @ (H4 - K @ K) @ B1inv)
    B2 = sqrt_psd(Lambda)

    Q1 = H @ Q0 @ B1inv
    q1_orth = float(np.linalg.norm(Q1.conj().T @ Q1 - np.eye(m)))
    C1 = B1inv @ H4
    recurrence = B1 @ K + Lambda @ B1
    recurrence_err = float(np.linalg.norm(C1 - recurrence))

    R_direct = H3Q0 - Q1 @ C1
    Rg_direct = herm(R_direct.conj().T @ R_direct)
    Rg_moment = herm(H6 - C1.conj().T @ C1)
    residual_cross = float(np.linalg.norm(Rg_direct - Rg_moment))

    P21 = B2 @ B1
    P21inv = np.linalg.solve(P21, np.eye(m))
    recovered = herm(P21inv.conj().T @ Rg_moment @ P21inv)
    target = herm(B3_true.conj().T @ B3_true)

    b1_err = float(np.linalg.norm(B1 - B1_true) / np.linalg.norm(B1_true))
    b2_err = float(np.linalg.norm(B2 - B2_true) / np.linalg.norm(B2_true))
    b3gram_err = float(np.linalg.norm(recovered - target) / np.linalg.norm(target))
    min_recovered = float(np.min(np.linalg.eigvalsh(recovered)))

    checks = {
        "Q1_orthonormal": q1_orth < TOL,
        "C1_recurrence_identity": recurrence_err < TOL,
        "direct_residual_equals_moment_residual": residual_cross < TOL,
        "B1_recovered": b1_err < TOL,
        "B2_recovered": b2_err < TOL,
        "B3dagB3_recovered": b3gram_err < TOL,
        "B3dagB3_positive": min_recovered > 0,
    }
    return {
        "status": "independent dense parity block-Lanczos B3 formula regression",
        "passed": bool(all(checks.values())),
        "seed": seed,
        "block_dimension": m,
        "errors": {
            "Q1_orthogonality": q1_orth,
            "C1_recurrence": recurrence_err,
            "direct_vs_moment_R3gram": residual_cross,
            "B1_relative": b1_err,
            "B2_relative": b2_err,
            "B3dagB3_relative": b3gram_err,
        },
        "recovered_B3dagB3_min_eigenvalue": min_recovered,
        "checks": checks,
        "scope_note": "Algebraic regression only; no BQG microscopic coefficient is supplied by this test.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
