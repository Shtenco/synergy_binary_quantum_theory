#!/usr/bin/env python3
"""Fail-closed finite-dimensional master-constraint/common-kernel certificate.

Infrastructure only: this accepts explicitly supplied matrices. It does NOT
load or infer the full W5 Lorentzian constraints and cannot certify BQG physics.
All arithmetic is exact SymPy arithmetic; floating-point inputs are rejected.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Sequence
import sympy as sp


def _exact_matrix(value) -> sp.Matrix:
    matrix = sp.Matrix(value)
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError("floating-point entries are forbidden; use exact rationals/algebraics")
    return matrix.applyfunc(sp.simplify)


def _positive_definite_hermitian(metric: sp.Matrix) -> bool:
    if metric.rows != metric.cols or any(sp.simplify(x) != 0 for x in metric - metric.H):
        return False
    # Sylvester criterion for Hermitian matrices: leading principal minors > 0.
    for k in range(1, metric.rows + 1):
        det = sp.simplify(sp.expand_complex(metric[:k, :k].det()))
        if det.is_real is not True or det.is_positive is not True:
            return False
    return True


def _zero(matrix: sp.Matrix) -> bool:
    return all(sp.simplify(sp.expand_complex(x)) == 0 for x in matrix)


def build_master_constraint(constraints: Sequence[sp.Matrix], label_metric: sp.Matrix,
                            input_metric: sp.Matrix | None = None) -> dict:
    """Build M=C†(G⊗I_K)C and the exact orthogonal common-kernel projector.

    Every C_A must map the same finite input H to the same finite output fiber K.
    G must be exactly Hermitian positive definite.
    """
    if not constraints:
        raise ValueError("at least one constraint matrix is required")
    Cs = [_exact_matrix(C) for C in constraints]
    G = _exact_matrix(label_metric)
    n, k, m = Cs[0].cols, Cs[0].rows, len(Cs)
    Q = sp.eye(n) if input_metric is None else _exact_matrix(input_metric)
    if Q.shape != (n, n) or not _positive_definite_hermitian(Q):
        raise ValueError("input metric is not provably Hermitian positive definite")
    if n == 0 or k == 0:
        raise ValueError("input and output dimensions must be positive")
    if any(C.cols != n or C.rows != k for C in Cs):
        raise ValueError("all constraints must share the same input H and output fiber K")
    if G.shape != (m, m):
        raise ValueError("label metric dimension must equal number of constraints")
    if not _positive_definite_hermitian(G):
        raise ValueError("label metric is not provably Hermitian positive definite")

    stacked = sp.zeros(m * k, n)
    for a, C in enumerate(Cs):
        stacked[a*k:(a+1)*k, :] = C
    pairing = sp.kronecker_product(G, sp.eye(k))
    quadratic_form = sp.simplify(stacked.H * pairing * stacked)
    # The coordinate matrix of the self-adjoint operator is Q^{-1} A.
    # It is Q-Hermitian, not generally Hermitian in the Euclidean metric.
    master = sp.simplify(Q.inv() * quadratic_form)
    if not _zero(Q * master - master.H * Q):
        raise ArithmeticError("constructed master operator is not self-adjoint in input metric")

    kernel = stacked.nullspace()
    rank_c, rank_m = stacked.rank(), master.rank()
    if rank_c != rank_m:
        raise ArithmeticError("rank(C†(G⊗I)C) != rank(C); common-kernel identity failed")
    if kernel:
        N = sp.Matrix.hstack(*kernel)
        projector = sp.simplify(N * (N.H * Q * N).inv() * N.H * Q)
    else:
        N = sp.zeros(n, 0)
        projector = sp.zeros(n, n)

    checks = {
        "master_metric_self_adjoint": _zero(Q * master - master.H * Q),
        "master_annihilates_kernel": _zero(master * N),
        "projector_metric_self_adjoint": _zero(projector.H * Q - Q * projector),
        "projector_idempotent": _zero(projector * projector - projector),
        "constraints_annihilate_projector": _zero(stacked * projector),
        "master_annihilates_projector": _zero(master * projector),
        "rank_nullity": len(kernel) == n - rank_c,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"finite projector certificate failed: {checks}")

    positive_eigenvalues = []
    for eigenvalue, multiplicity in master.eigenvals().items():
        ev = sp.simplify(sp.expand_complex(eigenvalue))
        if ev.is_zero is True:
            continue
        if ev.is_positive is not True:
            raise ArithmeticError(f"eigenvalue not provably positive: {ev}")
        positive_eigenvalues.extend([ev] * multiplicity)

    r = len(kernel)
    return {
        "input_dimension": n, "output_fiber_dimension": k, "constraint_count": m,
        "input_metric": Q,
        "constraint_rank": rank_c, "master_rank": rank_m,
        "physical_kernel_dimension": r,
        "physical_projective_dimension": r - 1 if r else None,
        "physical_fs_real_dimension": 2 * (r - 1) if r else None,
        "master_matrix": master, "projector_matrix": projector,
        "kernel_basis": N, "positive_eigenvalues": tuple(positive_eigenvalues),
        "checks": checks,
        "claim_boundary": (
            "Exact finite certificate for supplied matrices and declared input metric only. It does not "
            "establish that they are the full BQG/W5 Lorentzian constraints, nor "
            "a continuum limit, HDA, TT graviton sector, or Einstein IR limit."
        ),
    }


def _serialize_matrix(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


def run_reference_cases() -> dict:
    C1, C2 = sp.Matrix([[1, 0, 0]]), sp.Matrix([[0, 1, 0]])
    result = build_master_constraint([C1, C2], sp.eye(2))
    return {
        "status": "exact finite master-kernel reference gate",
        "passed": all(result["checks"].values()),
        "physical_kernel_dimension": result["physical_kernel_dimension"],
        "master_rank": result["master_rank"],
        "master_matrix": _serialize_matrix(result["master_matrix"]),
        "projector_matrix": _serialize_matrix(result["projector_matrix"]),
        "checks": result["checks"], "claim_boundary": result["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = run_reference_cases()
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
