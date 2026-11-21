#!/usr/bin/env python3
"""Exact positive control for boundary overlap with an enlarged zero sector.

The control proves a key linear-algebra fact used after the certified full-rank
32D logical master result: a boundary-restricted master can be positive definite
while the boundary still has nonzero overlap with a zero mode of the enlarged
master.

It also checks the heat-kernel limit P0 = lim_{tau->infty} exp(-tau M) and the
first block-Krylov closure for the minimal 2x2 example.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def run():
    M = np.array([[1.0, -1.0], [-1.0, 1.0]])
    e1 = np.array([[1.0], [0.0]])

    ev, U = np.linalg.eigh(M)
    zero = np.abs(ev) < 1e-12
    P0 = (U[:, zero]) @ (U[:, zero]).T

    boundary_master = float((e1.T @ M @ e1)[0, 0])
    boundary_projector = float((e1.T @ P0 @ e1)[0, 0])

    heat_rows = []
    for tau in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        H = expm(-tau * M)
        x = float((e1.T @ H @ e1)[0, 0])
        heat_rows.append({"tau": tau, "boundary_heat_kernel": x, "error_to_projector": abs(x - boundary_projector)})

    # Block-Krylov from the one-column boundary block closes the full 2D space
    # after one residual because M e1 is not collinear with e1.
    q0 = e1.copy()
    a0 = float((q0.T @ M @ q0)[0, 0])
    r = M @ q0 - q0 * a0
    b1 = float(np.linalg.norm(r))
    q1 = r / b1
    Q = np.column_stack([q0[:, 0], q1[:, 0]])
    T = Q.T @ M @ Q
    krylov_reconstruction_error = float(np.linalg.norm(Q @ T @ Q.T - M))
    krylov_orth_error = float(np.linalg.norm(Q.T @ Q - np.eye(2)))

    checks = {
        "master_positive_semidefinite": bool(np.min(ev) > -1e-12),
        "full_master_has_one_zero_mode": bool(np.sum(zero) == 1),
        "boundary_restriction_has_no_zero": bool(boundary_master > 1e-12),
        "boundary_overlap_with_full_zero_projector_nonzero": bool(abs(boundary_projector - 0.5) < 1e-12),
        "heat_kernel_converges_to_boundary_projector": bool(heat_rows[-1]["error_to_projector"] < 1e-12),
        "one_residual_closes_full_krylov_space": bool(krylov_reconstruction_error < 1e-12 and krylov_orth_error < 1e-12),
    }

    return {
        "status": "exact boundary-projector / block-Krylov positive control",
        "passed": bool(all(checks.values())),
        "master": M.tolist(),
        "master_eigenvalues": ev.tolist(),
        "boundary_vector": [1.0, 0.0],
        "boundary_restricted_master": boundary_master,
        "full_zero_projector": P0.tolist(),
        "boundary_projected_zero_projector": boundary_projector,
        "heat_kernel_scan": heat_rows,
        "block_krylov": {
            "a0": a0,
            "b1": b1,
            "T": T.tolist(),
            "orthonormality_error": krylov_orth_error,
            "reconstruction_error": krylov_reconstruction_error,
        },
        "checks": checks,
        "interpretation": (
            "An empty zero sector inside a boundary carrier does not imply zero overlap of that boundary with the full physical zero sector. "
            "The correct boundary observable is V0^dag P_phys V0, obtainable from the enlarged master heat kernel / block-Krylov representation."
        ),
        "claim_boundary": "Exact linear-algebra positive control only; it does not establish a BQG physical zero mode."
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
