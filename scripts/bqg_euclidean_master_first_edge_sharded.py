#!/usr/bin/env python3
"""Sharded first spectral edge of the actual finite K5 Euclidean BQG master.

For the complete 32-state all-j=1/2 logical seed block Q0, define

    M_E = sum_v H_v^dag H_v = sum_v H_v^2

using the repository's Hermitian sine ordering.  Column mode computes one exact
sparse state M_E|i>. Assembly forms A0=Q0^dag M_E Q0, mu2=(M_EQ0)^dag(M_EQ0),
and the first block-Lanczos residual Gram

    R1^dag R1 = mu2 - A0^dag A0 = B1^dag B1.

The calculation is sharded only for execution; the assembled operator and
stopping rule are identical to the monolithic definition.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_logical_anisotropy_gate as AN
import peter_weyl_euclidean_sine_ordering_gate as SINE
import peter_weyl_higher_shell_lambda_gate as HS
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

N = 32
TOL = 1e-10
JMAX2 = 5


def herm(M):
    return (M + M.conj().T) / 2


def sparse_add(dst, src, scale=1.0, tol=1e-11):
    for k, a in src.items():
        z = dst.get(k, 0j) + scale * a
        if abs(z) > tol:
            dst[k] = z
        elif k in dst:
            del dst[k]


def sparse_inner(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum(np.conj(v) * b.get(k, 0j) for k, v in a.items())


def sparse_norm(a):
    return math.sqrt(float(sum(abs(v) ** 2 for v in a.values())))


def gram(cols_a, cols_b=None):
    if cols_b is None:
        cols_b = cols_a
    G = np.zeros((len(cols_a), len(cols_b)), dtype=complex)
    for i, a in enumerate(cols_a):
        for j, b in enumerate(cols_b):
            G[i, j] = sparse_inner(a, b)
    return G


def combine(columns, coeffs, tol=1e-12):
    out = {}
    for col, c in zip(columns, coeffs):
        if abs(c) > tol:
            sparse_add(out, col, complex(c))
    return out


def logical_basis():
    keys = []
    labels = []
    for env in AN.ENV_STATES:
        for pair in AN.PAIR_STATES:
            keys.append(AN.logical_key(pair[0], pair[1], env))
            labels.append({"environment_K234": list(env), "pair_K01": list(pair)})
    if len(keys) != N:
        raise RuntimeError(f"expected {N} logical states, got {len(keys)}")
    return keys, labels


def apply_master(state):
    out = {}
    for v in PW.VERT:
        hv = SINE.safe_H_sine(state, v, JMAX2)
        h2v = SINE.safe_H_sine(hv, v, JMAX2)
        sparse_add(out, h2v)
    return out


def compute_column(index: int):
    if not 0 <= index < N:
        raise ValueError(f"column outside [0,{N-1}]: {index}")
    ZVM.patch_and_clear()
    keys, labels = logical_basis()
    q = {keys[index]: 1.0 + 0j}
    mq = apply_master(q)
    rows = HS.state_to_rows(mq)
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return {
        "status": "exact sharded finite K5 Euclidean master column",
        "column": index,
        "label": labels[index],
        "operator": "M_E=sum_v (H_E,v^sine)^dag H_E,v^sine=sum_v (H_E,v^sine)^2",
        "constraint_vertices": list(PW.VERT),
        "Jmax_used": JMAX2 / 2,
        "source_commit": os.environ.get("GITHUB_SHA", "LOCAL_OR_UNKNOWN"),
        "support": len(mq),
        "norm": sparse_norm(mq),
        "max_spin": HS.max_spin(mq),
        "state_sha256": hashlib.sha256(payload).hexdigest(),
        "master_state": rows,
    }


def load_columns(directory: Path):
    files = sorted(directory.glob("master_column_*.json"))
    cols = {}
    for p in files:
        d = json.loads(p.read_text(encoding="utf-8"))
        i = int(d["column"])
        if i in cols:
            raise RuntimeError(f"duplicate master column {i}")
        cols[i] = d
    missing = [i for i in range(N) if i not in cols]
    if missing:
        raise RuntimeError(f"missing master columns: {missing}")
    return [cols[i] for i in range(N)]


def assemble(directory: Path):
    ordered = load_columns(directory)
    keys, labels = logical_basis()
    q0 = [{key: 1.0 + 0j} for key in keys]
    mq0 = [HS.rows_to_state(d["master_state"]) for d in ordered]

    q0_gram = gram(q0)
    q0_orth_error = float(np.linalg.norm(q0_gram - np.eye(N)))
    A0_raw = gram(q0, mq0)
    A0_herm_error = float(np.linalg.norm(A0_raw - A0_raw.conj().T))
    A0 = herm(A0_raw)
    a0_eigs = np.linalg.eigvalsh(A0).real
    mu2 = herm(gram(mq0))

    residual = []
    for j, mq in enumerate(mq0):
        r = dict(mq)
        for i, qi in enumerate(q0):
            sparse_add(r, qi, -A0[i, j])
        residual.append(r)

    R1_direct = herm(gram(residual))
    R1_moment = herm(mu2 - A0.conj().T @ A0)
    moment_identity_error = float(
        np.linalg.norm(R1_direct - R1_moment)
        / max(np.linalg.norm(R1_direct), 1.0)
    )

    ew, U = np.linalg.eigh(R1_direct)
    scale = max(float(np.max(np.abs(ew))), float(np.linalg.norm(mu2)), 1.0)
    rank_tol = TOL * scale
    keep = np.where(ew > rank_tol)[0]
    rank = int(len(keep))
    negative_min = float(np.min(ew))

    supports = [int(d["support"]) for d in ordered]
    norms = [float(d["norm"]) for d in ordered]
    max_spins = [float(d["max_spin"]) for d in ordered]
    residual_norm = math.sqrt(max(float(np.trace(R1_direct).real), 0.0))
    mq_norm = math.sqrt(sum(x * x for x in norms))
    residual_relative = residual_norm / max(mq_norm, 1e-30)

    recurrence_error = 0.0
    q1_orth_error = 0.0
    q0q1_error = 0.0
    B1_shape = [0, N]
    if rank:
        lam = ew[keep]
        Ur = U[:, keep]
        q1 = [combine(residual, Ur[:, a] / math.sqrt(float(lam[a]))) for a in range(rank)]
        B1 = np.diag(np.sqrt(lam)) @ Ur.conj().T
        B1_shape = list(B1.shape)
        q1_orth_error = float(np.linalg.norm(gram(q1) - np.eye(rank)))
        q0q1_error = float(np.linalg.norm(gram(q0, q1)))
        errors = []
        for j, mq in enumerate(mq0):
            recon = {}
            for i, qi in enumerate(q0):
                sparse_add(recon, qi, A0[i, j])
            for a, qa in enumerate(q1):
                sparse_add(recon, qa, B1[a, j])
            diff = dict(mq)
            sparse_add(diff, recon, -1.0)
            errors.append(sparse_norm(diff) / max(sparse_norm(mq), 1e-30))
        recurrence_error = max(errors, default=0.0)

    regulator_safe = max(max_spins, default=0.0) <= 1.5 + 1e-12
    psd_ok = (
        float(np.min(a0_eigs)) > -2e-8 * max(float(np.max(np.abs(a0_eigs))), 1.0)
        and negative_min > -2e-8 * scale
    )
    algebra_pass = bool(
        q0_orth_error < 1e-12
        and A0_herm_error < 2e-8
        and psd_ok
        and moment_identity_error < 2e-7
        and q1_orth_error < 2e-7
        and q0q1_error < 2e-7
        and recurrence_error < 2e-7
        and regulator_safe
    )

    closure_certified = bool(algebra_pass and rank == 0)
    if not algebra_pass:
        science_status = "INVALID_EUCLIDEAN_MASTER_FIRST_EDGE"
        closure_statement = "At least one algebra/PSD/regulator guard failed; no closure statement is certified."
    elif rank == 0:
        science_status = "FINITE_EUCLIDEAN_MASTER_HISTORY_CLOSED_DEPTH_0"
        closure_statement = "B1 has zero resolved rank with all guards passing: the 32-state seed sector is invariant under M_E."
    else:
        science_status = "FINITE_EUCLIDEAN_MASTER_HISTORY_OPEN_AFTER_B1"
        closure_statement = "B1 has nonzero resolved rank with all guards passing: M_E generates a new shell and the master Lanczos chain must continue."

    return {
        "status": "sharded actual K5 Euclidean master first spectral edge",
        "passed": algebra_pass,
        "science_status": science_status,
        "operator": "M_E=sum_v (H_E,v^sine)^dag H_E,v^sine=sum_v (H_E,v^sine)^2",
        "seed": "complete 32-state all-j=1/2 logical K5 sector",
        "seed_dimension": N,
        "constraint_vertices": len(PW.VERT),
        "labels": labels,
        "Jmax_used": JMAX2 / 2,
        "provenance": {
            "column_source_commits": {str(d["column"]): d.get("source_commit") for d in ordered},
            "column_state_sha256": {str(d["column"]): d["state_sha256"] for d in ordered},
        },
        "regulator": {
            "master_image_max_spin": max(max_spins, default=0.0),
            "master_image_inside_cutoff": regulator_safe,
            "support_min": min(supports, default=0),
            "support_max": max(supports, default=0),
            "norm_min": min(norms, default=0.0),
            "norm_max": max(norms, default=0.0),
        },
        "A0": {
            "hermiticity_error": A0_herm_error,
            "eigenvalue_min": float(np.min(a0_eigs)),
            "eigenvalue_max": float(np.max(a0_eigs)),
        },
        "first_edge": {
            "mu2_minus_A0dagA0_relative_identity_error": moment_identity_error,
            "R1_gram_min_eigenvalue": negative_min,
            "R1_gram_max_eigenvalue": float(np.max(ew)),
            "rank_tolerance": float(rank_tol),
            "B1_rank": rank,
            "B1_shape": B1_shape,
            "Q1_orthogonality_error": q1_orth_error,
            "Q0_Q1_overlap_error": q0q1_error,
            "recurrence_relative_error": recurrence_error,
            "next_residual_relative": residual_relative,
            "closure_certified": closure_certified,
        },
        "closure_statement": closure_statement,
        "claim_boundary": (
            "Finite Euclidean common-kernel master spectral data only. This does not by itself close "
            "the Lorentzian/global refinement-compatible physical history, connected W[J], physical omega, "
            "or the graviton/scalar effective kernel."
        ),
    }


def write(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--column", type=int)
    mode.add_argument("--assemble-dir", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.column is not None:
        out = compute_column(args.column)
    else:
        out = assemble(args.assemble_dir)
    write(args.output, out)
    if args.column is not None:
        print(json.dumps({k: out[k] for k in ("column", "support", "norm", "max_spin", "state_sha256")}, indent=2))
        return 0
    print(json.dumps({"passed": out["passed"], "science_status": out["science_status"], "first_edge": out["first_edge"], "regulator": out["regulator"]}, indent=2))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
