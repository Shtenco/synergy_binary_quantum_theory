#!/usr/bin/env python3
"""Exact finite positive control for physical boundary source dressing.

The gate checks the operator ordering

    M -> P0 -> boundary Gram -> projected/whitened observables -> Z[J] -> W[J]

in a model where the boundary-compressed master is strictly positive even
though the boundary has nonzero overlap with the enlarged physical zero sector.

This is an operator-algebra/source-dressing control only. It is not a BQG
physical projector, cosmological vacuum, dark-sector prediction or propagator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

TOL = 2.0e-11


def invsqrt_on_support(M: np.ndarray, rel_tol: float = 1.0e-12):
    H = (M + M.conj().T) / 2.0
    w, U = np.linalg.eigh(H)
    scale = max(float(np.max(np.abs(w))), 1.0)
    keep = w > rel_tol * scale
    if not np.any(keep):
        raise ValueError("empty support")
    Q = U[:, keep]
    X = (Q * (1.0 / np.sqrt(w[keep]))) @ Q.conj().T
    P = Q @ Q.conj().T
    return X, P, w, keep


def opnorm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A, 2))


def run() -> dict:
    e = np.eye(4)
    p1 = (e[:, 0] + e[:, 2]) / np.sqrt(2.0)
    p2 = (e[:, 1] + e[:, 3]) / np.sqrt(2.0)
    q1 = (e[:, 0] - e[:, 2]) / np.sqrt(2.0)
    q2 = (e[:, 1] - e[:, 3]) / np.sqrt(2.0)

    P0 = np.outer(p1, p1) + np.outer(p2, p2)
    M = 2.0 * np.outer(q1, q1) + 5.0 * np.outer(q2, q2)
    B = e[:, :2]

    OX = np.outer(p1, p2) + np.outer(p2, p1)
    OZ = np.outer(p1, p1) - np.outer(p2, p2)

    evals = np.linalg.eigvalsh((M + M.T) / 2.0)
    gap = float(np.min(evals[evals > 1.0e-12]))
    MB = B.T @ M @ B
    MB_evals = np.linalg.eigvalsh(MB)
    G0 = B.T @ P0 @ B

    G0_i2, G0_supp, _, _ = invsqrt_on_support(G0)
    AX0 = B.T @ P0 @ OX @ P0 @ B
    AZ0 = B.T @ P0 @ OZ @ P0 @ B
    OXbar0 = G0_i2 @ AX0 @ G0_i2
    OZbar0 = G0_i2 @ AZ0 @ G0_i2

    X2 = np.array([[0.0, 1.0], [1.0, 0.0]])
    Z2 = np.array([[1.0, 0.0], [0.0, -1.0]])
    I2 = np.eye(2)

    exact_checks = {
        "master_positive_semidefinite": bool(float(np.min(evals)) > -TOL),
        "full_zero_rank_two": bool(np.sum(np.abs(evals) < 1.0e-12) == 2),
        "boundary_compressed_master_positive_definite": bool(float(np.min(MB_evals)) > 1.0e-12),
        "boundary_zero_projector_is_half_identity": bool(np.linalg.norm(G0 - 0.5 * I2) < TOL),
        "whitened_X_is_Pauli_X": bool(np.linalg.norm(OXbar0 - X2) < TOL),
        "whitened_Z_is_Pauli_Z": bool(np.linalg.norm(OZbar0 - Z2) < TOL),
        "physical_support_is_full_boundary_rank": bool(np.linalg.norm(G0_supp - I2) < TOL),
    }

    rows = []
    for tau in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        K = expm(-tau * M)
        Kh = expm(-0.5 * tau * M)
        Gt = B.T @ K @ B
        Gt_i2, _, _, _ = invsqrt_on_support(Gt)
        AXt = B.T @ Kh @ OX @ Kh @ B
        AZt = B.T @ Kh @ OZ @ Kh @ B
        Xbar = Gt_i2 @ AXt @ Gt_i2
        Zbar = Gt_i2 @ AZt @ Gt_i2

        # Zero-source Hessian of log[(1/2) Tr exp(jx Xbar + jz Zbar)]
        # at j=0. For the maximally mixed support this is the symmetrized
        # covariance matrix.
        ops = (Xbar, Zbar)
        means = np.array([np.trace(O).real / 2.0 for O in ops])
        H = np.zeros((2, 2), float)
        for a in range(2):
            for b in range(2):
                sym = 0.5 * (ops[a] @ ops[b] + ops[b] @ ops[a])
                H[a, b] = float(np.trace(sym).real / 2.0 - means[a] * means[b])

        wrong = expm(-tau * MB)
        heat_err = opnorm(Gt - G0)
        heat_bound = float(np.exp(-tau * gap))
        source_err = opnorm(H - I2)

        rows.append({
            "tau": tau,
            "boundary_heat_Gram": Gt.tolist(),
            "boundary_heat_error_to_P0": heat_err,
            "gap_bound_exp_minus_tau_Delta": heat_bound,
            "compressed_master_heat_kernel": wrong.tolist(),
            "compressed_master_heat_norm": opnorm(wrong),
            "whitened_X": Xbar.tolist(),
            "whitened_Z": Zbar.tolist(),
            "source_Hessian": H.tolist(),
            "source_Hessian_error_to_identity": source_err,
        })

    last = rows[-1]

    # Independent finite-source check at a nonzero source point.
    jx, jz = 0.37, -0.51
    R = np.sqrt(jx * jx + jz * jz)
    exact_Z = float(np.cosh(R))
    exact_W = float(np.log(exact_Z))
    source_matrix = jx * OXbar0 + jz * OZbar0
    numerical_Z = float(np.trace(expm(source_matrix)).real / 2.0)
    numerical_W = float(np.log(numerical_Z))

    # Wrong compressed-master projector is empty: exp(-tau B^T M B) -> 0.
    # Full-space projection retains G0 = I/2.
    checks = dict(exact_checks)
    checks.update({
        "heat_boundary_Gram_respects_gap_bound": bool(all(r["boundary_heat_error_to_P0"] <= r["gap_bound_exp_minus_tau_Delta"] + 5e-12 for r in rows)),
        "heat_boundary_Gram_converges": bool(last["boundary_heat_error_to_P0"] < 1.0e-7),
        "whitened_source_Hessian_converges": bool(last["source_Hessian_error_to_identity"] < 2.0e-6),
        "compressed_master_heat_dies_while_full_projector_overlap_survives": bool(last["compressed_master_heat_norm"] < 5.0e-4 and opnorm(G0) > 0.49),
        "finite_source_matches_cosh": bool(abs(numerical_Z - exact_Z) < TOL and abs(numerical_W - exact_W) < TOL),
    })

    return {
        "status": "exact boundary-projector source-dressing positive control",
        "passed": bool(all(checks.values())),
        "science_status": "FINITE_OPERATOR_SOURCE_DRESSING_CONTROL",
        "master_eigenvalues": [float(x) for x in evals],
        "master_gap": gap,
        "boundary_compressed_master": MB.tolist(),
        "boundary_compressed_master_eigenvalues": [float(x) for x in MB_evals],
        "boundary_physical_Gram_G0": G0.tolist(),
        "exact_whitened_observables": {
            "X": OXbar0.tolist(),
            "Z": OZbar0.tolist(),
        },
        "heat_kernel_staircase": rows,
        "finite_source_check": {
            "jx": jx,
            "jz": jz,
            "radius": float(R),
            "Z_expected_cosh": exact_Z,
            "Z_numerical": numerical_Z,
            "W_expected_log_cosh": exact_W,
            "W_numerical": numerical_W,
        },
        "checks": checks,
        "interpretation": (
            "A strictly positive boundary-compressed master can coexist with nonzero boundary overlap with the full zero projector. "
            "The legal source order is full master -> zero projector/heat filter -> boundary Gram and projected observable -> support whitening -> Z[J]."
        ),
        "production_target": (
            "Replace the 4D control by the 32-column q=2 boundary block and the parity-complete Peter-Weyl master; evaluate G_tau and source matrices with block Krylov rather than a dense ambient diagonalization."
        ),
        "claim_boundary": (
            "No BQG physical zero mode, physical clock, cosmological vacuum, dark matter, dark energy, lensing signal or physical frequency is established by this control."
        ),
    }


def main() -> int:
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
