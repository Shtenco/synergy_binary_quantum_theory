#!/usr/bin/env python3
"""Executable adapter from a finite constraint family to the existing BQG history/source stack.

This gate does NOT invent a new physicalization formalism. It factors the already
registered route into one reusable numerical interface:

    {C_A} -> M=sum C_A^dagger C_A -> P_phys
          -> physical observables -> connected source Hessian
          -> metric response -> tangent Gamma^(2).

The selftest uses the existing C8 relational positive control, but crucially the
projector used downstream is reconstructed *only* from the combined constraint
C=I-G. It is then regression-tested against the analytic relational projector.
The same functions are intended to accept the future common-habitat BQG bundle
{H_E,v, H_L,v, D_I} without changing the source/observable machinery.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

TOL = 2e-10


def spectral_projector_from_constraints(constraints, zero_rtol=1e-10):
    if not constraints:
        raise ValueError("empty constraint family")
    mats = [np.asarray(c, dtype=complex) for c in constraints]
    n = mats[0].shape[0]
    if any(c.shape != (n, n) for c in mats):
        raise ValueError("all constraints must be square and act on one common habitat")
    M = np.zeros((n, n), dtype=complex)
    for c in mats:
        M += c.conj().T @ c
    M = 0.5 * (M + M.conj().T)
    evals, evecs = np.linalg.eigh(M)
    scale = max(float(np.max(np.abs(evals))), 1.0)
    zero = np.abs(evals) <= zero_rtol * scale
    Q0 = evecs[:, zero]
    P0 = Q0 @ Q0.conj().T
    positive = evals[~zero]
    gap = float(np.min(positive)) if positive.size else None
    return M, P0, Q0, evals, gap


def physical_compression(Q0, O):
    O = np.asarray(O, dtype=complex)
    return Q0.conj().T @ O @ Q0


def maximally_mixed_connected_hessian(Q0, observables):
    r = Q0.shape[1]
    if r == 0:
        raise ValueError("physical zero sector is empty")
    ops = [physical_compression(Q0, O) for O in observables]
    means = np.array([np.trace(o).real / r for o in ops])
    C = np.zeros((len(ops), len(ops)), dtype=float)
    for a, oa in enumerate(ops):
        for b, ob in enumerate(ops):
            sym = 0.5 * (oa @ ob + ob @ oa)
            C[a, b] = float((np.trace(sym).real / r) - means[a] * means[b])
    return C, means, ops


def c8_control():
    n = 8
    S = np.zeros((n, n), dtype=complex)
    for t in range(n):
        S[(t + 1) % n, t] = 1.0
    J = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=complex)
    G = np.kron(S, J)
    C = np.eye(16, dtype=complex) - G

    blocks = []
    Rt = np.eye(2, dtype=complex)
    for _ in range(n):
        blocks.append(Rt / np.sqrt(n))
        Rt = Rt @ J
    V = np.vstack(blocks)
    P_rel = V @ V.conj().T

    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)

    def relationalize(O):
        out = np.zeros((16, 16), dtype=complex)
        Rt = np.eye(2, dtype=complex)
        for t in range(n):
            block = Rt @ O @ np.linalg.inv(Rt)
            out[2*t:2*t+2, 2*t:2*t+2] = block
            Rt = Rt @ J
        return out

    return C, G, V, P_rel, [relationalize(X), relationalize(Z), relationalize(Y)]


def q2_metric_jacobian():
    s3 = np.sqrt(3.0)
    MX = np.array([[s3/2, 0, s3/2], [0, -s3/2, -s3/2], [s3/2, -s3/2, 0]], dtype=float)
    MZ = np.array([[0.5, 1, -0.5], [1, 0.5, -0.5], [-0.5, -0.5, -1]], dtype=float)
    return np.column_stack([MX.reshape(-1), MZ.reshape(-1)])


def run():
    C, G, V, P_rel, rel_ops = c8_control()
    M, P0, Q0, evals, gap = spectral_projector_from_constraints([C])

    checks = {}
    checks["constraint_reconstructed_projector_rank_is_2"] = Q0.shape[1] == 2
    checks["constraint_reconstructed_projector_matches_analytic_Prel"] = float(np.linalg.norm(P0 - P_rel)) < 2e-9
    checks["constraint_annihilates_physical_projector"] = float(np.linalg.norm(C @ P0)) < 2e-9
    checks["projector_is_idempotent"] = float(np.linalg.norm(P0 @ P0 - P0)) < 2e-9
    checks["projector_is_hermitian"] = float(np.linalg.norm(P0 - P0.conj().T)) < 2e-9

    Csrc, means, compressed = maximally_mixed_connected_hessian(Q0, rel_ops)
    checks["physical_source_means_vanish"] = float(np.max(np.abs(means))) < 2e-9
    checks["XYZ_connected_hessian_is_identity"] = float(np.linalg.norm(Csrc - np.eye(3))) < 2e-9

    B = q2_metric_jacobian()
    shape_cov = Csrc[:2, :2]
    Cmetric = B @ shape_cov @ B.T
    Cplus = np.linalg.pinv(Cmetric, rcond=1e-12)
    eig = np.linalg.eigvalsh(Cmetric)
    nz = eig[eig > 1e-9]
    gamma_nz = np.linalg.eigvalsh(Cplus)
    gamma_nz = gamma_nz[gamma_nz > 1e-9]
    checks["metric_response_rank_is_2"] = np.linalg.matrix_rank(Cmetric, tol=1e-9) == 2
    checks["metric_response_nonzero_eigenvalues_are_9_over_2"] = float(np.max(np.abs(nz - 4.5))) < 2e-9
    checks["metric_Gamma2_nonzero_eigenvalues_are_2_over_9"] = float(np.max(np.abs(gamma_nz - (2.0/9.0)))) < 2e-9

    pauli_sq = [float(np.linalg.norm(o @ o - np.eye(2))) for o in compressed]
    checks["compressed_relational_sources_square_to_identity"] = max(pauli_sq) < 2e-9

    # Normalize numpy scalar booleans before JSON serialization.
    checks = {k: bool(v) for k, v in checks.items()}
    passed = bool(all(checks.values()))
    return {
        "status": "constraint-to-existing-history/source physicalization adapter selftest",
        "passed": passed,
        "input_scope": "C8 relational positive control reconstructed from C=I-G; production input is the common-habitat BQG {H_E,H_L,D} family",
        "hilbert_dimension": 16,
        "constraint_count": 1,
        "physical_dimension": int(Q0.shape[1]),
        "master_eigenvalue_min": float(np.min(evals)),
        "master_first_positive_gap": gap,
        "projector_error_to_existing_relational_projector": float(np.linalg.norm(P0 - P_rel)),
        "connected_source_hessian_XZY": Csrc.tolist(),
        "metric_response_nonzero_eigenvalues": [float(x) for x in nz],
        "metric_Gamma2_nonzero_eigenvalues": [float(x) for x in gamma_nz],
        "checks": checks,
        "production_contract": {
            "required_same_habitat_constraints": ["H_E,v", "H_L,v", "D_I (unless independently reduced/proved redundant)"],
            "required_metadata": ["regulator", "basis/habitat identity", "normalization/order convention", "constraint-family hash"],
            "downstream_reused_without_redefinition": ["physical projector", "relational/gauge-invariant source interface", "Z[J]", "W[J]", "metric response/Gamma2", "scalar/TT/FLRW sector maps", "real-observable dictionary"],
            "forbidden_shortcuts": [
                "boundary-compressed master used as full physical projector",
                "constraint spectral z renamed physical omega",
                "missing H_L or D replaced by phenomenological fit",
                "source/background normalization chosen after seeing DM/DE outcome"
            ]
        },
        "next_actual_input": "Replace C8 C=I-G by the common regulated BQG graph-changing constraint bundle and reuse this downstream interface unchanged."
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=Path("verification_results/BQG_PHYSICAL_HISTORY_ADAPTER.json"))
    args = ap.parse_args()
    out = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
