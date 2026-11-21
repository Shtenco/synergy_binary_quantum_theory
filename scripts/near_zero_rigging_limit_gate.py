#!/usr/bin/env python3
"""Deterministic positive/negative controls for an asymptotic master zero sector.

Checks the non-arbitrary heat schedule
    tau = 1/sqrt(lambda_r lambda_{r+1})
when lambda_r/lambda_{r+1}->0, plus two negative controls:
(1) low/high eigenvalues collapse without scale separation;
(2) the low eigenspace rotates and the boundary projector fails to converge.

This is a mathematical limit control, not a BQG physical-sector result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def projector(cols: np.ndarray) -> np.ndarray:
    return cols @ cols.conj().T


def positive_sequence() -> dict:
    e = np.eye(4)
    p1 = (e[:, 0] + e[:, 2]) / np.sqrt(2.0)
    p2 = (e[:, 1] + e[:, 3]) / np.sqrt(2.0)
    q1 = (e[:, 0] - e[:, 2]) / np.sqrt(2.0)
    q2 = (e[:, 1] - e[:, 3]) / np.sqrt(2.0)
    U = np.column_stack([p1, p2, q1, q2])
    Plow = projector(U[:, :2])
    B = e[:, :2]
    Gtarget = B.T @ Plow @ B

    rows = []
    for n in range(1, 8):
        eps = 2.0 ** (-n)
        vals = np.array([eps**2, 2.0 * eps**2, 1.0, 3.0])
        M = (U * vals) @ U.T
        lr, lnext = vals[1], vals[2]
        ratio = lr / lnext
        tau = 1.0 / np.sqrt(lr * lnext)
        H = expm(-tau * M)
        heat_err = float(np.linalg.norm(H - Plow, 2))
        bound = max(1.0 - np.exp(-tau * lr), np.exp(-tau * lnext))
        Gheat = B.T @ H @ B
        Gerr = float(np.linalg.norm(Gheat - Gtarget, 2))
        rows.append({
            "level": n,
            "epsilon": eps,
            "lambda_r": float(lr),
            "lambda_r_plus_1": float(lnext),
            "ratio": float(ratio),
            "tau": float(tau),
            "tau_lambda_r": float(tau * lr),
            "tau_lambda_r_plus_1": float(tau * lnext),
            "heat_to_low_projector_error": heat_err,
            "analytic_error_bound": float(bound),
            "boundary_Gram_error": Gerr,
        })

    ratios = np.array([r["ratio"] for r in rows])
    errs = np.array([r["heat_to_low_projector_error"] for r in rows])
    gerrs = np.array([r["boundary_Gram_error"] for r in rows])
    return {
        "target_boundary_Gram": Gtarget.tolist(),
        "rows": rows,
        "ratio_monotone_to_zero": bool(np.all(np.diff(ratios) < 0) and ratios[-1] < 1e-3),
        "heat_error_monotone": bool(np.all(np.diff(errs) < 0)),
        "boundary_error_monotone": bool(np.all(np.diff(gerrs) < 0)),
        "final_heat_error": float(errs[-1]),
        "final_boundary_error": float(gerrs[-1]),
        "all_errors_inside_bound": bool(all(r["heat_to_low_projector_error"] <= r["analytic_error_bound"] + 2e-12 for r in rows)),
    }


def no_separation_negative() -> dict:
    # Both candidate low and next eigenvalue collapse with fixed ratio 1/2.
    # The prescribed tau cannot become identity on the low mode while killing
    # the next mode: both exponents remain O(1).
    rows = []
    for n in range(1, 8):
        eps = 2.0 ** (-n)
        lr = eps
        lnext = 2.0 * eps
        ratio = lr / lnext
        tau = 1.0 / np.sqrt(lr * lnext)
        low_deviation = 1.0 - np.exp(-tau * lr)
        high_survival = np.exp(-tau * lnext)
        err = max(low_deviation, high_survival)
        rows.append({
            "level": n,
            "epsilon": eps,
            "ratio": ratio,
            "tau_lambda_low": float(tau * lr),
            "tau_lambda_next": float(tau * lnext),
            "separation_error": float(err),
        })
    errors = np.array([r["separation_error"] for r in rows])
    return {
        "rows": rows,
        "ratio_is_constant": bool(max(abs(r["ratio"] - 0.5) for r in rows) < 1e-15),
        "error_fails_to_converge_to_zero": bool(errors[-1] > 0.2 and abs(errors[-1] - errors[0]) < 1e-12),
    }


def rotating_subspace_negative() -> dict:
    e = np.eye(2)
    B = e[:, [0]]
    rows = []
    overlaps = []
    for n in range(1, 9):
        eps = 2.0 ** (-n)
        theta = 0.0 if n % 2 == 0 else np.pi / 3.0
        u = np.array([np.cos(theta), np.sin(theta)])
        v = np.array([-np.sin(theta), np.cos(theta)])
        P = np.outer(u, u)
        M = eps**2 * np.outer(u, u) + np.outer(v, v)
        # eigenvalue separation is excellent, but the subspace itself rotates.
        overlap = float((B.T @ P @ B)[0, 0])
        overlaps.append(overlap)
        rows.append({
            "level": n,
            "epsilon": eps,
            "lambda_ratio": eps**2,
            "theta": float(theta),
            "boundary_low_projector_overlap": overlap,
            "master_min_eigenvalue": float(np.min(np.linalg.eigvalsh(M))),
        })
    return {
        "rows": rows,
        "overlap_min": float(min(overlaps)),
        "overlap_max": float(max(overlaps)),
        "boundary_projector_fails_to_converge": bool(max(overlaps) - min(overlaps) > 0.5),
    }


def run() -> dict:
    pos = positive_sequence()
    neg1 = no_separation_negative()
    neg2 = rotating_subspace_negative()
    checks = {
        "positive_ratio_separates": pos["ratio_monotone_to_zero"],
        "positive_heat_projector_converges": pos["heat_error_monotone"] and pos["final_heat_error"] < 0.02,
        "positive_boundary_Gram_converges": pos["boundary_error_monotone"] and pos["final_boundary_error"] < 0.01,
        "positive_respects_analytic_bound": pos["all_errors_inside_bound"],
        "negative_no_separation_is_rejected": neg1["ratio_is_constant"] and neg1["error_fails_to_converge_to_zero"],
        "negative_rotating_subspace_is_rejected": neg2["boundary_projector_fails_to_converge"],
    }
    return {
        "status": "near-zero master rigging-limit positive/negative control",
        "passed": bool(all(checks.values())),
        "science_status": "ASYMPTOTIC_PROJECTOR_LIMIT_CONTROL",
        "positive_separated_low_cluster": pos,
        "negative_collapsing_without_separation": neg1,
        "negative_rotating_low_subspace": neg2,
        "checks": checks,
        "theorem_tested": (
            "If lambda_r/lambda_{r+1}->0, tau=1/sqrt(lambda_r lambda_{r+1}) makes tau lambda_r->0 and tau lambda_{r+1}->infinity, so the heat kernel separates the low cluster without a fitted threshold. Eigenvalue separation alone is insufficient if the low projector fails to converge."
        ),
        "claim_boundary": "Mathematical regulator-limit control only; no BQG near-zero cluster is supplied by this gate.",
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
