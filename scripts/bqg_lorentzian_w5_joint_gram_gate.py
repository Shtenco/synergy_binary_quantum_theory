#!/usr/bin/env python3
"""W5 Lorentzian joint-Gram gate.

Input JSON must contain raw, row-aligned complex amplitude maps A and B:
{
  "input_basis_K": [0,2,4,6,8],
  "output_basis_id": "same ordered output basis for A and B",
  "A": [[[re,im], ... five columns ...], ... output rows ...],
  "B": [[[re,im], ... five columns ...], ... output rows ...],
  "normalization": {"c_E": [1,0], "c_L": [1,0]}
}
The raw maps, not just A^*A, are required to determine A^*B.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np


def zscalar(value):
    if isinstance(value, (int, float)):
        return complex(value)
    if isinstance(value, list) and len(value) == 2:
        return complex(float(value[0]), float(value[1]))
    raise ValueError("complex scalar must be a number or [real, imag]")


def zmatrix(value, name):
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 3 or arr.shape[2] != 2:
        raise ValueError(f"{name} must have shape (output_rows, 5, 2) as [real, imag]")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("--output", type=Path, default=Path("BQG_LORENTZIAN_W5_JOINT_GRAM.json"))
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--atol", type=float, default=1e-12)
    args = ap.parse_args()

    data = json.loads(args.input.read_text(encoding="utf-8"))
    basis = data.get("input_basis_K")
    if basis != [0, 2, 4, 6, 8]:
        raise ValueError("input_basis_K must exactly equal [0,2,4,6,8]")
    if not data.get("output_basis_id"):
        raise ValueError("output_basis_id is required to certify row alignment")
    A = zmatrix(data["A"], "A")
    B = zmatrix(data["B"], "B")
    if A.shape != B.shape:
        raise ValueError(f"A and B must have identical shape; got {A.shape} and {B.shape}")
    if data.get("A_output_basis_id", data["output_basis_id"]) != data["output_basis_id"]:
        raise ValueError("A output basis identifier mismatch")
    if data.get("B_output_basis_id", data["output_basis_id"]) != data["output_basis_id"]:
        raise ValueError("B output basis identifier mismatch")

    norm = data.get("normalization", {})
    cE = zscalar(norm.get("c_E", [1, 0]))
    cL = zscalar(norm.get("c_L", [1, 0]))
    GEE = A.conj().T @ A
    GEL = A.conj().T @ B
    GLE = B.conj().T @ A
    GLL = B.conj().T @ B
    Hsum = cE * A + cL * B
    Gsum = Hsum.conj().T @ Hsum
    Hstack = np.vstack((A, B))
    Gstack = Hstack.conj().T @ Hstack

    eig_EE, herm_EE = hermitian_eigenvalues(GEE)
    eig_LL, herm_LL = hermitian_eigenvalues(GLL)
    eig_sum, herm_sum = hermitian_eigenvalues(Gsum)
    eig_stack, herm_stack = hermitian_eigenvalues(Gstack)
    sv_sum = np.linalg.svd(Hsum, compute_uv=False)
    sv_stack = np.linalg.svd(Hstack, compute_uv=False)
    tol_sum = args.atol + args.rtol * max(1.0, float(sv_sum[0]))
    tol_stack = args.atol + args.rtol * max(1.0, float(sv_stack[0]))

    result = {
        "status": "COMPUTED",
        "claim_scope": "finite W5 maps only; not a full-theory or continuum result",
        "input_basis_K": basis,
        "output_basis_id": data["output_basis_id"],
        "output_rows": int(A.shape[0]),
        "input_columns": 5,
        "normalization": {"c_E": [cE.real, cE.imag], "c_L": [cL.real, cL.imag]},
        "gram_matrices": {
            "A_dagger_A": encode_matrix(GEE),
            "A_dagger_B": encode_matrix(GEL),
            "B_dagger_A": encode_matrix(GLE),
            "B_dagger_B": encode_matrix(GLL),
            "Hsum_dagger_Hsum": encode_matrix(Gsum),
            "stacked_dagger_stacked": encode_matrix(Gstack)
        },
        "spectra": {
            "A_dagger_A_eigenvalues": eig_EE.tolist(),
            "B_dagger_B_eigenvalues": eig_LL.tolist(),
            "Hsum_gram_eigenvalues": eig_sum.tolist(),
            "stacked_gram_eigenvalues": eig_stack.tolist(),
            "Hsum_singular_values": sv_sum.tolist(),
            "stacked_singular_values": sv_stack.tolist()
        },
        "diagnostics": {
            "hermiticity_residuals": {"GEE": herm_EE, "GLL": herm_LL, "Gsum": herm_sum, "Gstack": herm_stack},
            "Hsum_kernel_dimension_at_tolerance": int(np.count_nonzero(sv_sum <= tol_sum)),
            "stacked_kernel_dimension_at_tolerance": int(np.count_nonzero(sv_stack <= tol_stack)),
            "Hsum_singular_tolerance": tol_sum,
            "stacked_singular_tolerance": tol_stack,
            "Hsum_sigma_min": float(sv_sum[-1]),
            "stacked_sigma_min": float(sv_stack[-1])
        },
        "interpretation": {
            "sum_operator": "A common zero exists in the represented W5 sector iff Hsum has a nontrivial numerical kernel under the declared tolerance.",
            "stacked_operator": "A common zero of A and B exists iff the vertically stacked map has a nontrivial numerical kernel.",
            "warning": "A^*B requires raw A and B with the same output basis, row ordering, phases, and normalization. A^*A alone does not determine A^*B."
        }
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "A_shape": list(A.shape), "sigma_min_sum": result["diagnostics"]["Hsum_sigma_min"], "sum_kernel_dim": result["diagnostics"]["Hsum_kernel_dimension_at_tolerance"], "sigma_min_stacked": result["diagnostics"]["stacked_sigma_min"], "stacked_kernel_dim": result["diagnostics"]["stacked_kernel_dimension_at_tolerance"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"GATE_ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
