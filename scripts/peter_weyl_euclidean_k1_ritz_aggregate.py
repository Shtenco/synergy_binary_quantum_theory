#!/usr/bin/env python3
"""Aggregate five target-node shards into the Euclidean K1 Ritz master.

Consumes five NPZ files from peter_weyl_euclidean_k1_ritz_shard.py, checks that
the independently recomputed first-layer Grams agree, sums D^(w), whitens on
supp(G), and reports the generated-layer generalized/Ritz spectrum. No target
spectrum is hard-coded.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

TOL = 5.0e-8


def load_shards(input_dir: Path):
    files = sorted(input_dir.glob("k1_ritz_node_*.npz"))
    if len(files) != 5:
        raise ValueError(f"expected 5 node shards, got {len(files)}: {files}")
    rows = []
    seen = set()
    for p in files:
        z = np.load(p)
        node = int(z["target_node"][0])
        if node in seen:
            raise ValueError(f"duplicate target node {node}")
        seen.add(node)
        rows.append((node, z["G"], z["D"], z["labels"]))
    rows.sort(key=lambda x: x[0])
    if seen != set(range(5)):
        raise ValueError(f"missing target nodes: {set(range(5)) - seen}")
    return rows


def run(input_dir: Path) -> dict:
    shards = load_shards(input_dir)
    G0 = shards[0][1]
    labels0 = shards[0][3]

    g_errors = []
    label_equal = []
    for node, G, D, labels in shards:
        rel = float(np.linalg.norm(G - G0) / max(np.linalg.norm(G0), 1e-30))
        g_errors.append(rel)
        label_equal.append(bool(np.array_equal(labels, labels0)))

    G = (G0 + G0.conj().T) / 2.0
    Dparts = [(D + D.conj().T) / 2.0 for _, _, D, _ in shards]
    D = sum(Dparts)
    D = (D + D.conj().T) / 2.0

    gev, U = np.linalg.eigh(G)
    gscale = max(float(np.max(np.abs(gev))), 1.0)
    gtol = 1.0e-10 * gscale
    keep = gev > gtol
    rankG = int(np.sum(keep))
    if rankG == 0:
        raise RuntimeError("empty G support")

    Ur = U[:, keep]
    gr = gev[keep]
    Dr = Ur.conj().T @ D @ Ur
    invs = 1.0 / np.sqrt(gr)
    M_odd = (invs[:, None] * Dr) * invs[None, :]
    M_odd = (M_odd + M_odd.conj().T) / 2.0
    odd_ev = np.linalg.eigvalsh(M_odd)

    # Reconstruct bare five-node boundary master from the five diagonal 32x32
    # blocks of the first-layer Gram.
    M_B = np.zeros((32, 32), complex)
    for v in range(5):
        M_B += G[32*v:32*(v+1), 32*v:32*(v+1)]
    M_B = (M_B + M_B.conj().T) / 2.0
    b_ev = np.linalg.eigvalsh(M_B)

    full_ev = np.sort(np.concatenate([b_ev, odd_ev]))
    bmin = float(np.min(b_ev))
    omin = float(np.min(odd_ev))
    kmin = float(np.min(full_ev))
    dress_ratio = kmin / bmin if bmin > 0 else None

    D_ev = np.linalg.eigvalsh(D)
    dscale = max(float(np.max(np.abs(D_ev))), 1.0)
    oddscale = max(float(np.max(np.abs(odd_ev))), 1.0)

    part_traces = np.array([np.trace(X).real for X in Dparts], float)
    trmean = max(abs(float(np.mean(part_traces))), 1.0)
    target_trace_spread = float((np.max(part_traces) - np.min(part_traces)) / trmean)

    # Numerical nullities are descriptive only.
    btol = 1e-10 * max(float(np.max(np.abs(b_ev))), 1.0)
    otol = 1e-10 * oddscale
    ftol = 1e-10 * max(float(np.max(np.abs(full_ev))), 1.0)

    checks = {
        "five_shards_present": len(shards) == 5,
        "labels_identical_across_shards": bool(all(label_equal)),
        "first_Gram_recomputed_consistently": max(g_errors) < 2e-12,
        "G_positive_semidefinite": float(np.min(gev)) > -TOL * gscale,
        "D_positive_semidefinite": float(np.min(D_ev)) > -TOL * dscale,
        "whitened_odd_Ritz_positive_semidefinite": float(np.min(odd_ev)) > -TOL * oddscale,
        "boundary_master_positive_semidefinite": float(np.min(b_ev)) > -TOL * max(float(np.max(np.abs(b_ev))), 1.0),
    }

    return {
        "status": "aggregated Euclidean K1 Peter-Weyl Ritz master",
        "passed": bool(all(checks.values())),
        "science_status": "FINITE_EUCLIDEAN_K1_RITZ_MASTER",
        "Jmax": 2.5,
        "first_layer": {
            "labelled_columns": 160,
            "rank": rankG,
            "nullity": 160 - rankG,
            "Gram_rank_tolerance": gtol,
            "smallest_positive_Gram_eigenvalue": float(np.min(gr)),
            "largest_Gram_eigenvalue": float(np.max(gr)),
            "support_condition_number": float(np.max(gr) / np.min(gr)),
            "cross_shard_relative_errors": g_errors,
            "total_K1_dimension_by_parity": 32 + rankG,
        },
        "second_layer_master": {
            "summed_D_min_eigenvalue": float(np.min(D_ev)),
            "summed_D_max_eigenvalue": float(np.max(D_ev)),
            "target_node_D_traces": [float(x) for x in part_traces],
            "target_node_trace_relative_spread": target_trace_spread,
        },
        "boundary_master": {
            "dimension": 32,
            "min_eigenvalue": bmin,
            "max_eigenvalue": float(np.max(b_ev)),
            "rank": int(np.sum(b_ev > btol)),
            "nullity": int(np.sum(b_ev <= btol)),
            "eigenvalues": [float(x) for x in b_ev],
        },
        "generated_odd_Ritz_master": {
            "dimension": rankG,
            "min_eigenvalue": omin,
            "max_eigenvalue": float(np.max(odd_ev)),
            "rank": int(np.sum(odd_ev > otol)),
            "nullity": int(np.sum(odd_ev <= otol)),
            "eigenvalues": [float(x) for x in odd_ev],
        },
        "full_K1_Ritz": {
            "dimension": int(len(full_ev)),
            "min_eigenvalue": kmin,
            "max_eigenvalue": float(np.max(full_ev)),
            "rank": int(np.sum(full_ev > ftol)),
            "nullity": int(np.sum(full_ev <= ftol)),
            "eigenvalues": [float(x) for x in full_ev],
            "dressing_minimum_ratio_to_boundary": float(dress_ratio) if dress_ratio is not None else None,
        },
        "checks": checks,
        "interpretation": (
            "The generalized eigenvalues are the Euclidean normal-master Ritz spectrum on the first odd constraint-generated Peter-Weyl layer. "
            "A decrease of the minimum is a dressing diagnostic only; a near-zero candidate requires tighter pruning, tangential constraints, Lorentzian completion and refinement/rigging stability."
        ),
        "claim_boundary": (
            "Normal Euclidean K1 Ritz calculation only. No complete physical projector, dark scalar, physical omega, dark matter or dark energy is inferred."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = run(args.input_dir)
    text = json.dumps(out, indent=2)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
