#!/usr/bin/env python3
"""Extend the certified Peter-Weyl H=H_E0+H_E1 block-Lanczos chain to B3.

Input column artifacts contain exact sparse states a_i=H|i>, b_i=H^2|i>.
Column mode applies the same frozen H once more: c_i=H^3|i>.
Assembly reconstructs K, Lambda, B1, B2 and derives B3^dag B3.

The B3 residual is computed twice:
  1. from the moment identity H6-C1^dag C1;
  2. directly in sparse state space after reconstructing Q1 and subtracting Q1 C1.
A closure claim is legal only if those two routes agree and every provenance,
positivity and regulator guard passes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_higher_shell_lambda_gate as HS

N = HS.NLOGICAL
TOL = 1e-10


def herm(M):
    return (M + M.conj().T) / 2


def matrix_json(M):
    return [[[float(z.real), float(z.imag)] for z in row] for row in np.asarray(M)]


def matrix_from_json(rows):
    return np.array([[complex(float(z[0]), float(z[1])) for z in row] for row in rows], dtype=complex)


def sqrt_psd(M, negative_tol=2e-8):
    M = herm(M)
    ev, U = np.linalg.eigh(M)
    scale = max(float(np.max(np.abs(ev))), 1.0)
    if float(np.min(ev)) < -negative_tol * scale:
        raise RuntimeError(f"matrix is not PSD: min={float(np.min(ev))}, scale={scale}")
    S = (U * np.sqrt(np.clip(ev, 0.0, None))) @ U.conj().T
    return herm(S), ev


def sparse_linear_combo(columns, coeffs):
    out = {}
    for col, c in zip(columns, coeffs):
        if abs(c) > 1e-14:
            HS.AN.sparse_add(out, col, complex(c))
    return out


def sparse_subtract_combo(state, basis, coeffs):
    out = dict(state)
    for col, c in zip(basis, coeffs):
        if abs(c) > 1e-14:
            HS.AN.sparse_add(out, col, -complex(c))
    return out


def extend_column(source: Path):
    raw = source.read_bytes()
    d = json.loads(raw)
    idx = int(d["column"])
    if not 0 <= idx < N:
        raise ValueError(f"column outside [0,{N-1}]: {idx}")
    HS.AN.ZVM.patch_and_clear()
    second = HS.rows_to_state(d["second_state"])
    third = HS.apply_H_state(second)
    return {
        "status": "exact Peter-Weyl spectral-history H3 column",
        "column": idx,
        "label": d["label"],
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "operator": "H=H_E,0+H_E,1 with frozen sine ordering",
        "Jmax_used": HS.JMAX2_SECOND_HIT_SAFE / 2,
        "first_order_projection_norm": float(d["first_order_projection_norm"]),
        "first_support": int(d["first_support"]),
        "second_support": int(d["second_support"]),
        "third_support": len(third),
        "first_max_spin": float(d["first_max_spin"]),
        "second_max_spin": float(d["second_max_spin"]),
        "third_max_spin": float(HS.max_spin(third)),
        "first_state": d["first_state"],
        "second_state": d["second_state"],
        "third_state": HS.state_to_rows(third),
    }


def load_columns(directory: Path):
    files = sorted(directory.glob("column_*_b3.json"))
    cols = {}
    for p in files:
        d = json.loads(p.read_text(encoding="utf-8"))
        i = int(d["column"])
        if i in cols:
            raise RuntimeError(f"duplicate column {i}")
        cols[i] = d
    missing = [i for i in range(N) if i not in cols]
    if missing:
        raise RuntimeError(f"missing B3 columns: {missing}")
    return [cols[i] for i in range(N)]


def assemble(directory: Path, reference: Path | None = None):
    ordered = load_columns(directory)
    first = [HS.rows_to_state(d["first_state"]) for d in ordered]
    second = [HS.rows_to_state(d["second_state"]) for d in ordered]
    third = [HS.rows_to_state(d["third_state"]) for d in ordered]

    K = HS.gram(first)
    H4 = HS.gram(second)
    H6 = HS.gram(third)

    B1, ek = sqrt_psd(K, negative_tol=1e-10)
    kscale = max(float(np.max(np.abs(ek))), 1.0)
    krank = int(np.sum(ek > 1e-10 * kscale))
    if krank != N:
        raise RuntimeError(f"K lost full rank: {krank}")

    B1inv = np.linalg.solve(B1, np.eye(N))
    M2 = herm(H4 - K @ K)
    Lambda = herm(B1inv @ M2 @ B1inv)
    B2, el = sqrt_psd(Lambda, negative_tol=2e-7)
    lscale = max(float(np.max(np.abs(el))), 1.0)
    lrank = int(np.sum(el > 1e-10 * lscale))
    if lrank != N:
        raise RuntimeError(f"Lambda lost full rank: {lrank}")

    # Q1 is reconstructed directly from the exact first-shell sparse states:
    # H Q0 = Q1 B1  =>  Q1 = (H Q0) B1^-1.
    q1 = [sparse_linear_combo(first, B1inv[:, a]) for a in range(N)]
    q1_gram = HS.gram(q1)
    q1_orth_error = float(np.linalg.norm(q1_gram - np.eye(N)))

    # C1 = Q1^dag H^3 Q0.  The moment route and recurrence route must agree.
    C1 = B1inv @ H4
    C1_recurrence = B1 @ K + Lambda @ B1
    c1_identity_error = float(np.linalg.norm(C1 - C1_recurrence))

    # Moment residual.
    R3gram_moment = herm(H6 - C1.conj().T @ C1)

    # Independent direct sparse residual: r_j = H^3|j> - Q1 C1[:,j].
    r3_direct = [sparse_subtract_combo(third[j], q1, C1[:, j]) for j in range(N)]
    R3gram_direct = HS.gram(r3_direct)
    r3_direct_norms = [HS.AN.sparse_norm(r) for r in r3_direct]
    r3_direct_supports = [len(r) for r in r3_direct]
    r3_cross_error = float(
        np.linalg.norm(R3gram_direct - R3gram_moment)
        / max(np.linalg.norm(R3gram_direct), 1.0)
    )
    R3gram = herm((R3gram_direct + R3gram_moment) / 2)
    er = np.linalg.eigvalsh(R3gram).real
    rscale = max(float(np.max(np.abs(er))), 1.0)

    P21 = B2 @ B1
    P21inv = np.linalg.solve(P21, np.eye(N))
    B3gram = herm(P21inv.conj().T @ R3gram @ P21inv)
    e3 = np.linalg.eigvalsh(B3gram).real
    e3scale = max(float(np.max(np.abs(e3))), 1.0)
    rank_tol = TOL * e3scale
    rank3 = int(np.sum(e3 > rank_tol))
    B3, _ = sqrt_psd(B3gram, negative_tol=5e-7)
    b3_reconstruction_error = float(np.linalg.norm(B3.conj().T @ B3 - B3gram))

    ref_errors = {}
    if reference is not None:
        ref = json.loads(reference.read_text(encoding="utf-8"))
        mats = ref["matrices_common_logical_basis"]
        for name, cur in (("K", K), ("Lambda", Lambda), ("B1", B1), ("B2", B2)):
            old = matrix_from_json(mats[name])
            ref_errors[name] = float(np.linalg.norm(cur - old) / max(np.linalg.norm(old), 1e-30))

    first_proj_max = max(float(d["first_order_projection_norm"]) for d in ordered)
    third_max_spin = max(float(d["third_max_spin"]) for d in ordered)
    third_support_min = min(int(d["third_support"]) for d in ordered)
    third_support_max = max(int(d["third_support"]) for d in ordered)

    psd_ok = float(np.min(er)) > -5e-7*rscale and float(np.min(e3)) > -5e-7*e3scale
    reference_ok = (not ref_errors) or max(ref_errors.values()) < 2e-9
    regulator_safe = third_max_spin <= 2.0 + 1e-12
    identities_ok = (
        q1_orth_error < 2e-8
        and c1_identity_error < 2e-8 * max(float(np.linalg.norm(C1)), 1.0)
        and r3_cross_error < 5e-8
    )
    passed = bool(
        first_proj_max < 1e-12
        and psd_ok
        and reference_ok
        and regulator_safe
        and identities_ok
        and b3_reconstruction_error < 2e-8
    )

    closure_certified = bool(passed and rank3 == 0)
    if not passed:
        science_status = "INVALID_B3_EVIDENCE"
        closure_statement = (
            "B3 evidence failed at least one algebra/provenance/regulator guard; "
            "no closure or non-closure claim is certified."
        )
    elif rank3 == 0:
        science_status = "FINITE_EUCLIDEAN_HISTORY_CLOSED_DEPTH_2"
        closure_statement = (
            "B3=0 within the frozen rank resolution and all independent guards pass: "
            "the finite Euclidean parity Krylov history closes at depth 2."
        )
    else:
        science_status = "FINITE_EUCLIDEAN_HISTORY_OPEN_AFTER_B3"
        closure_statement = (
            "B3 has nonzero resolved rank with all guards passing: the finite Euclidean "
            "parity Krylov history is not closed at depth 2; the next exact shell is required."
        )

    return {
        "status":"actual Peter-Weyl parity block-Lanczos spectral-history B3 gate",
        "passed":passed,
        "science_status":science_status,
        "operator":"H=H_E,0+H_E,1",
        "seed_dimension":N,
        "column_count":len(ordered),
        "provenance":{
            "source_column_sha256":{str(d["column"]):d["source_sha256"] for d in ordered},
            "reference_matrix_relative_errors":ref_errors,
        },
        "regulator":{
            "Jmax_used":HS.JMAX2_SECOND_HIT_SAFE/2,
            "first_order_projection_max":first_proj_max,
            "third_max_spin":third_max_spin,
            "third_support_min":third_support_min,
            "third_support_max":third_support_max,
            "third_hit_inside_cutoff":regulator_safe,
        },
        "moments":{
            "K_equals_PH2P_eigen_min":float(np.min(ek)),
            "K_equals_PH2P_eigen_max":float(np.max(ek)),
            "H4_frobenius_norm":float(np.linalg.norm(H4)),
            "H6_frobenius_norm":float(np.linalg.norm(H6)),
            "H6_hermiticity_error":float(np.linalg.norm(H6-H6.conj().T)),
        },
        "independent_sparse_audit":{
            "Q1_orthogonality_error":q1_orth_error,
            "R3_direct_vs_moment_relative_error":r3_cross_error,
            "R3_direct_norm_min":min(r3_direct_norms, default=0.0),
            "R3_direct_norm_max":max(r3_direct_norms, default=0.0),
            "R3_direct_support_min":min(r3_direct_supports, default=0),
            "R3_direct_support_max":max(r3_direct_supports, default=0),
        },
        "B2":{
            "Lambda_eigen_min":float(np.min(el)),
            "Lambda_eigen_max":float(np.max(el)),
            "C1_equals_B1inv_H4_recurrence_error":c1_identity_error,
        },
        "B3":{
            "R3gram_min_eigenvalue":float(np.min(er)),
            "R3gram_max_eigenvalue":float(np.max(er)),
            "B3dagB3_min_eigenvalue":float(np.min(e3)),
            "B3dagB3_max_eigenvalue":float(np.max(e3)),
            "B3dagB3_frobenius_norm":float(np.linalg.norm(B3gram)),
            "rank_tolerance":float(rank_tol),
            "rank":rank3,
            "reconstruction_error":b3_reconstruction_error,
            "closure_certified":closure_certified,
        },
        "closure_statement":closure_statement,
        "matrices_common_logical_basis":{
            "basis_order":"environment K234 major, pair K01 minor",
            "B3dagB3":matrix_json(B3gram),
            "B3":matrix_json(B3),
        },
        "claim_boundary":(
            "This is an exact/finite Euclidean constraint spectral-history result for H_E0+H_E1. "
            "It is not the global Lorentzian physical projector/history, physical omega, connected W[J], "
            "or graviton 1PI kernel."
        ),
    }


def write(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--extend-column",type=Path)
    mode.add_argument("--assemble-b3-dir",type=Path)
    ap.add_argument("--reference",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    if args.extend_column is not None:
        out=extend_column(args.extend_column)
        write(args.output,out)
        print(json.dumps({k:out[k] for k in ("column","second_support","third_support","third_max_spin")},indent=2))
        return 0
    out=assemble(args.assemble_b3_dir,args.reference)
    write(args.output,out)
    print(json.dumps({"passed":out["passed"],"science_status":out["science_status"],"B3":out["B3"],"regulator":out["regulator"]},indent=2))
    return 0 if out["passed"] else 1

if __name__=="__main__":
    raise SystemExit(main())
