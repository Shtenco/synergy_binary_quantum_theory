#!/usr/bin/env python3
"""Fail-closed W5 joint-Gram and common-kernel gate.

The input must carry raw, row-aligned complex maps A and B, explicit identities
for both output bases, and explicitly derived (never silently defaulted)
normalization coefficients. Rectangular-map nullities include structural
right-null directions when the output dimension is below five.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

W5_BASIS = [0, 2, 4, 6, 8]


def zscalar(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real or [real, imag] pair")
    if isinstance(value, (int, float)):
        z = complex(value)
    elif isinstance(value, list) and len(value) == 2:
        z = complex(float(value[0]), float(value[1]))
    else:
        raise ValueError(f"{name} must be a number or [real, imag]")
    if not (np.isfinite(z.real) and np.isfinite(z.imag)):
        raise ValueError(f"{name} must be finite")
    return z


def zmatrix(value, name):
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 3 or arr.shape[2] != 2:
        raise ValueError(f"{name} must have shape (output_rows, 5, 2) as [real, imag]")
    if arr.shape[0] == 0:
        raise ValueError(f"{name} must have at least one output row")
    if arr.shape[1] != 5:
        raise ValueError(f"{name} must have exactly five W5 input columns")
    if not np.isfinite(arr).all():
        raise ValueError(f"{name} contains non-finite entries")
    return arr[..., 0] + 1j * arr[..., 1]


def hermitian_eigenvalues(g):
    residual = float(np.linalg.norm(g - g.conj().T, ord=2))
    return np.linalg.eigvalsh((g + g.conj().T) / 2), residual


def encode_matrix(m):
    return [[[float(x.real), float(x.imag)] for x in row] for row in m]


def full_right_singular_values(matrix, input_columns, rtol, atol):
    """Return all input-space singular values, padding structural zeros."""
    sv = np.linalg.svd(matrix, compute_uv=False)
    if sv.size < input_columns:
        sv = np.pad(sv, (0, input_columns - sv.size))
    scale = float(sv[0]) if sv.size else 0.0
    tolerance = atol + rtol * max(1.0, scale)
    rank = int(np.count_nonzero(sv > tolerance))
    return sv, tolerance, rank, input_columns - rank


def compute_joint_gram(data, rtol=1e-10, atol=1e-12):
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0 or atol < 0:
        raise ValueError("rtol and atol must be finite and nonnegative")
    if rtol == 0 and atol == 0:
        raise ValueError("at least one of rtol or atol must be positive")
    if data.get("input_basis_K") != W5_BASIS:
        raise ValueError("input_basis_K must exactly equal [0,2,4,6,8]")

    common_id = data.get("output_basis_id")
    a_id = data.get("A_output_basis_id")
    b_id = data.get("B_output_basis_id")
    if not all(isinstance(x, str) and x.strip() for x in (common_id, a_id, b_id)):
        raise ValueError("output_basis_id, A_output_basis_id and B_output_basis_id are all required")
    if a_id != common_id or b_id != common_id:
        raise ValueError("A and B must explicitly identify the same ordered output basis")

    normalization = data.get("normalization")
    if not isinstance(normalization, dict) or "c_E" not in normalization or "c_L" not in normalization:
        raise ValueError("normalization.c_E and normalization.c_L must be explicit; implicit unit defaults are forbidden")
    if not isinstance(data.get("normalization_source"), str) or not data["normalization_source"].strip():
        raise ValueError("normalization_source must identify the frozen canonical convention or derivation")

    A = zmatrix(data["A"], "A")
    B = zmatrix(data["B"], "B")
    if A.shape != B.shape:
        raise ValueError(f"A and B must have identical shape; got {A.shape} and {B.shape}")

    cE = zscalar(normalization["c_E"], "normalization.c_E")
    cL = zscalar(normalization["c_L"], "normalization.c_L")

    GEE = A.conj().T @ A
    GEL = A.conj().T @ B
    GLE = B.conj().T @ A
    GLL = B.conj().T @ B
    Hsum = cE * A + cL * B
    Gsum = Hsum.conj().T @ Hsum
    Hstack = np.vstack((A, B))
    Gstack = Hstack.conj().T @ Hstack

    eig_EE, herm_EE = hermitian_eigenvalues(GEE)
    eig_EL, herm_EL = hermitian_eigenvalues(GEL + GLE)
    eig_LL, herm_LL = hermitian_eigenvalues(GLL)
    eig_sum, herm_sum = hermitian_eigenvalues(Gsum)
    eig_stack, herm_stack = hermitian_eigenvalues(Gstack)

    sv_sum, tol_sum, rank_sum, kernel_sum = full_right_singular_values(
        Hsum, A.shape[1], rtol, atol
    )
    sv_stack, tol_stack, rank_stack, kernel_stack = full_right_singular_values(
        Hstack, A.shape[1], rtol, atol
    )

    return {
        "status": "COMPUTED",
        "claim_scope": "finite W5 maps only; not a full-theory or continuum result",
        "input_basis_K": W5_BASIS,
        "output_basis_id": common_id,
        "A_output_basis_id": a_id,
        "B_output_basis_id": b_id,
        "output_rows": int(A.shape[0]),
        "input_columns": int(A.shape[1]),
        "normalization": {
            "c_E": [cE.real, cE.imag],
            "c_L": [cL.real, cL.imag],
            "source": data["normalization_source"],
        },
        "gram_matrices": {
            "A_dagger_A": encode_matrix(GEE),
            "A_dagger_B": encode_matrix(GEL),
            "B_dagger_A": encode_matrix(GLE),
            "B_dagger_B": encode_matrix(GLL),
            "Hsum_dagger_Hsum": encode_matrix(Gsum),
            "stacked_dagger_stacked": encode_matrix(Gstack),
        },
        "spectra": {
            "A_dagger_A_eigenvalues": eig_EE.tolist(),
            "A_dagger_B_plus_B_dagger_A_eigenvalues": eig_EL.tolist(),
            "B_dagger_B_eigenvalues": eig_LL.tolist(),
            "Hsum_gram_eigenvalues": eig_sum.tolist(),
            "stacked_gram_eigenvalues": eig_stack.tolist(),
            "Hsum_all_five_right_singular_values": sv_sum.tolist(),
            "stacked_all_five_right_singular_values": sv_stack.tolist(),
        },
        "diagnostics": {
            "hermiticity_residuals": {
                "GEE": herm_EE,
                "GEL_plus_GLE": herm_EL,
                "GLL": herm_LL,
                "Gsum": herm_sum,
                "Gstack": herm_stack,
            },
            "Hsum_rank_at_tolerance": rank_sum,
            "Hsum_kernel_dimension_at_tolerance": kernel_sum,
            "stacked_rank_at_tolerance": rank_stack,
            "stacked_kernel_dimension_at_tolerance": kernel_stack,
            "Hsum_singular_tolerance": tol_sum,
            "stacked_singular_tolerance": tol_stack,
            "Hsum_sigma_min_including_structural_zeros": float(sv_sum[-1]),
            "stacked_sigma_min_including_structural_zeros": float(sv_stack[-1]),
        },
        "interpretation": {
            "sum_operator": "A common zero exists in the represented W5 sector iff Hsum has a nontrivial numerical right kernel under the declared tolerance.",
            "stacked_operator": "A common zero of A and B exists iff the vertically stacked map has a nontrivial numerical right kernel.",
            "warning": "A_dagger_B requires raw A and B with identical output states, row ordering, phases and normalization. Euclidean Gram data alone cannot reconstruct it.",
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("--output", type=Path, default=Path("BQG_LORENTZIAN_W5_JOINT_GRAM.json"))
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--atol", type=float, default=1e-12)
    args = ap.parse_args()

    data = json.loads(args.input.read_text(encoding="utf-8"))
    result = compute_joint_gram(data, rtol=args.rtol, atol=args.atol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    d = result["diagnostics"]
    print(json.dumps({
        "status": result["status"],
        "A_shape": [result["output_rows"], result["input_columns"]],
        "Hsum_rank": d["Hsum_rank_at_tolerance"],
        "sum_kernel_dim": d["Hsum_kernel_dimension_at_tolerance"],
        "stacked_rank": d["stacked_rank_at_tolerance"],
        "stacked_kernel_dim": d["stacked_kernel_dimension_at_tolerance"],
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"GATE_ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
