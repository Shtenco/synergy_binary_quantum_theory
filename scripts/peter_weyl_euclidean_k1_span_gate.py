#!/usr/bin/env python3
"""Preregistered first Euclidean constraint-generated Peter-Weyl span.

Computes the 160x160 Gram matrix of H_v^E |b_i> for five K5 nodes and 32
all-j=1/2 logical boundary states. The Gram rank is a scientific output and is
not used as a pass target.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW

JMAX2 = 5
PRUNE = 1.0e-8
TOL = 3.0e-9


def sparse_inner(a: dict, b: dict) -> complex:
    if len(a) <= len(b):
        return sum(np.conj(x) * b.get(k, 0.0j) for k, x in a.items())
    return np.conj(sparse_inner(b, a))


def is_boundary_key(key) -> bool:
    spins, Ks = key
    return all(s == 1 for s in spins) and all(k in (0, 2) for k in Ks)


def boundary_projection_norm(state: dict) -> float:
    z = sum(abs(a) ** 2 for k, a in state.items() if is_boundary_key(k))
    return float(np.sqrt(z))


def run() -> dict:
    basis = PW.basis_full_jhalf()
    if len(basis) != 32:
        raise RuntimeError(f"expected 32 boundary states, got {len(basis)}")

    labels = []
    images = []
    supports = []
    boundary_overlap_max = 0.0

    for v in range(5):
        for i, key in enumerate(basis):
            out = PW.prune_state(PW.apply_H_cached_state({key: 1.0 + 0.0j}, v, JMAX2), PRUNE)
            labels.append((v, i))
            images.append(out)
            supports.append(len(out))
            boundary_overlap_max = max(boundary_overlap_max, boundary_projection_norm(out))

    n = len(images)
    G = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i, n):
            z = sparse_inner(images[i], images[j])
            G[i, j] = z
            G[j, i] = np.conj(z)

    herm = float(np.linalg.norm(G - G.conj().T))
    Gh = (G + G.conj().T) / 2.0
    ev = np.linalg.eigvalsh(Gh)
    scale = max(float(np.max(np.abs(ev))), 1.0)
    rank_tol = 1.0e-10 * scale
    rank = int(np.sum(ev > rank_tol))
    nullity = n - rank
    pos = ev[ev > rank_tol]

    # 5x5 block diagnostics, each block is 32x32.
    block_fro = np.zeros((5, 5), float)
    block_trace = np.zeros((5, 5), complex)
    node_blocks = []
    for v in range(5):
        for w in range(5):
            block = G[32 * v:32 * (v + 1), 32 * w:32 * (w + 1)]
            block_fro[v, w] = float(np.linalg.norm(block))
            block_trace[v, w] = np.trace(block)
        node_blocks.append(G[32 * v:32 * (v + 1), 32 * v:32 * (v + 1)])

    Mdiag = sum(node_blocks)
    Mdiag = (Mdiag + Mdiag.conj().T) / 2.0
    mev = np.linalg.eigvalsh(Mdiag)

    diagonal_sq = float(sum(np.linalg.norm(node_blocks[v]) ** 2 for v in range(5)))
    total_sq = float(np.linalg.norm(G) ** 2)
    cross_sq = max(total_sq - diagonal_sq, 0.0)

    diagonal_traces = np.array([np.trace(B).real for B in node_blocks])
    tr_mean = max(abs(float(np.mean(diagonal_traces))), 1.0)
    trace_rel_spread = float((np.max(diagonal_traces) - np.min(diagonal_traces)) / tr_mean)

    # Algebraic two-node combined-image reconstruction for three deterministic pairs.
    pair_checks = []
    pair_max_rel = 0.0
    for v, w in ((0, 1), (0, 2), (3, 4)):
        reconstructed = (
            G[32*v:32*(v+1), 32*v:32*(v+1)]
            + G[32*w:32*(w+1), 32*w:32*(w+1)]
            + G[32*v:32*(v+1), 32*w:32*(w+1)]
            + G[32*w:32*(w+1), 32*v:32*(v+1)]
        )
        direct = np.zeros((32, 32), complex)
        summed = []
        for i in range(32):
            acc = dict(images[32*v+i])
            for k, a in images[32*w+i].items():
                acc[k] = acc.get(k, 0.0j) + a
                if abs(acc[k]) < 1.0e-14:
                    del acc[k]
            summed.append(acc)
        for i in range(32):
            for j in range(i, 32):
                z = sparse_inner(summed[i], summed[j])
                direct[i, j] = z
                direct[j, i] = np.conj(z)
        rel = float(np.linalg.norm(direct - reconstructed) / max(np.linalg.norm(direct), 1.0e-30))
        pair_max_rel = max(pair_max_rel, rel)
        pair_checks.append({"nodes": [v, w], "relative_reconstruction_error": rel})

    checks = {
        "all_160_columns_evaluated": n == 160,
        "Gram_hermitian": herm < TOL,
        "Gram_positive_semidefinite": float(np.min(ev)) > -TOL * scale,
        "outgoing_odd_layer_orthogonal_to_boundary": boundary_overlap_max < 1.0e-10,
        "combined_pair_Gram_reconstruction": pair_max_rel < 2.0e-12,
        "node_diagonal_trace_covariance": trace_rel_spread < 3.0e-9,
    }

    return {
        "status": "preregistered Euclidean K1 Peter-Weyl outgoing-span Gram",
        "passed": bool(all(checks.values())),
        "science_status": "FINITE_CONSTRAINT_GENERATED_K1_SPAN",
        "Jmax": JMAX2 / 2.0,
        "prune_threshold": PRUNE,
        "boundary_dimension": 32,
        "labelled_outgoing_columns": n,
        "generated_odd_dimension": rank,
        "total_K1_dimension_by_parity": 32 + rank,
        "Gram": {
            "dimension": n,
            "rank": rank,
            "nullity": nullity,
            "rank_tolerance": rank_tol,
            "hermiticity_error": herm,
            "eigenvalue_min": float(np.min(ev)),
            "eigenvalue_max": float(np.max(ev)),
            "smallest_positive_eigenvalue": float(np.min(pos)) if len(pos) else None,
            "condition_number_on_support": float(np.max(pos) / np.min(pos)) if len(pos) else None,
            "eigenvalues": [float(x) for x in ev],
        },
        "diagonal_master_reconstructed_from_G1": {
            "eigenvalue_min": float(np.min(mev)),
            "eigenvalue_max": float(np.max(mev)),
            "trace": float(np.trace(Mdiag).real),
        },
        "block_frobenius_norms": block_fro.tolist(),
        "block_traces": [[[float(z.real), float(z.imag)] for z in row] for row in block_trace],
        "diagonal_block_trace_relative_spread": trace_rel_spread,
        "diagonal_Gram_weight_fraction": diagonal_sq / max(total_sq, 1.0e-30),
        "cross_node_Gram_weight_fraction": cross_sq / max(total_sq, 1.0e-30),
        "support_statistics": {
            "min": int(min(supports)),
            "max": int(max(supports)),
            "mean": float(np.mean(supports)),
        },
        "boundary_projection_max_norm": boundary_overlap_max,
        "pair_reconstruction_checks": pair_checks,
        "pair_reconstruction_max_relative_error": pair_max_rel,
        "checks": checks,
        "interpretation": (
            "rank(G1) is the actual dimension of the first odd Peter-Weyl layer generated from the frozen 32-state boundary by the five Euclidean constraints. "
            "The next Ritz/master calculation should use this measured span rather than the naive 160-column upper bound."
        ),
        "claim_boundary": (
            "This is a one-action Euclidean habitat compression diagnostic. It is not the enlarged zero projector, not a Lorentzian physical state and not a dark-sector observable."
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
