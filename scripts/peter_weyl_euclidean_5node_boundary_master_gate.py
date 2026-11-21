#!/usr/bin/env python3
"""Preregistered five-node Euclidean Peter-Weyl boundary master.

Builds the canonical identity-metric master compressed to the 32-dimensional
all-j=1/2 q=2 K5 boundary carrier:

    M_ij = sum_v < H_v b_i | H_v b_j >.

The outgoing states are the full sparse Peter-Weyl states. They are not first
projected back to the logical carrier. Rank/nullity is reported, not used as a
pass criterion.
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
TOL = 2.0e-9


def sparse_inner(a: dict, b: dict) -> complex:
    if len(a) > len(b):
        a, b = b, a
        return np.conj(sparse_inner(a, b))
    return sum(np.conj(x) * b.get(k, 0.0j) for k, x in a.items())


def sparse_norm(a: dict) -> float:
    return float(np.sqrt(max(sparse_inner(a, a).real, 0.0)))


def is_boundary_key(key) -> bool:
    spins, Ks = key
    return all(s == 1 for s in spins) and all(k in (0, 2) for k in Ks)


def projected_boundary_norm(state: dict) -> float:
    return sparse_norm({k: a for k, a in state.items() if is_boundary_key(k)})


def gram(images: list[dict]) -> np.ndarray:
    n = len(images)
    G = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(i, n):
            z = sparse_inner(images[i], images[j])
            G[i, j] = z
            G[j, i] = np.conj(z)
    return G


def run() -> dict:
    basis = PW.basis_full_jhalf()
    if len(basis) != 32:
        raise RuntimeError(f"expected 32 boundary states, got {len(basis)}")

    images_by_node: dict[int, list[dict]] = {}
    node_grams: dict[int, np.ndarray] = {}
    node_rows = []
    direct_return_max = 0.0

    for v in range(5):
        imgs = []
        supports = []
        norms = []
        for key in basis:
            out = PW.apply_H_cached_state({key: 1.0 + 0.0j}, v, JMAX2)
            out = PW.prune_state(out, PRUNE)
            imgs.append(out)
            supports.append(len(out))
            norms.append(sparse_norm(out))
            direct_return_max = max(direct_return_max, projected_boundary_norm(out))

        Gv = gram(imgs)
        images_by_node[v] = imgs
        node_grams[v] = Gv
        node_rows.append({
            "node": v,
            "trace": float(np.trace(Gv).real),
            "frobenius_norm": float(np.linalg.norm(Gv)),
            "support_min": int(min(supports)),
            "support_max": int(max(supports)),
            "support_mean": float(np.mean(supports)),
            "image_norm_min": float(min(norms)),
            "image_norm_max": float(max(norms)),
            "image_norm_mean": float(np.mean(norms)),
            "node_Gram_hermiticity_error": float(np.linalg.norm(Gv - Gv.conj().T)),
        })

    M = sum(node_grams.values())
    M = (M + M.conj().T) / 2.0
    herm = float(np.linalg.norm(M - M.conj().T))
    ev = np.linalg.eigvalsh(M)
    scale = max(float(np.max(np.abs(ev))), 1.0)
    rank_tol = 1.0e-10 * scale
    rank = int(np.sum(ev > rank_tol))
    nullity = len(ev) - rank
    positive = ev[ev > rank_tol]

    traces = np.array([r["trace"] for r in node_rows], float)
    trace_mean = max(abs(float(np.mean(traces))), 1.0)
    trace_rel_spread = float((np.max(traces) - np.min(traces)) / trace_mean)

    # An independent PSD witness: random deterministic coefficient vectors must
    # agree between c^dag M c and the sum of outgoing state norms.
    rng = np.random.default_rng(5092026)
    quadratic_checks = []
    quad_max_rel = 0.0
    for _ in range(5):
        c = rng.normal(size=32) + 1j * rng.normal(size=32)
        c /= np.linalg.norm(c)
        lhs = float(np.vdot(c, M @ c).real)
        rhs = 0.0
        for v in range(5):
            acc = {}
            for ci, img in zip(c, images_by_node[v]):
                for k, a in img.items():
                    acc[k] = acc.get(k, 0.0j) + ci * a
            rhs += sparse_inner(acc, acc).real
        rhs = float(rhs)
        rel = abs(lhs - rhs) / max(abs(rhs), 1.0e-30)
        quad_max_rel = max(quad_max_rel, rel)
        quadratic_checks.append({"matrix_form": lhs, "sum_image_norms": rhs, "relative_error": rel})

    checks = {
        "five_nodes_evaluated": len(node_rows) == 5,
        "master_hermitian": herm < TOL,
        "master_positive_semidefinite": float(np.min(ev)) > -TOL * scale,
        "direct_boundary_return_zero_by_parity": direct_return_max < 1.0e-10,
        "node_trace_permutation_covariance": trace_rel_spread < 2.0e-9,
        "quadratic_form_identity": quad_max_rel < 2.0e-10,
    }

    return {
        "status": "preregistered five-node Euclidean Peter-Weyl boundary master",
        "passed": bool(all(checks.values())),
        "science_status": "FINITE_EUCLIDEAN_BOUNDARY_MASTER",
        "Jmax": JMAX2 / 2.0,
        "prune_threshold": PRUNE,
        "boundary_dimension": 32,
        "constraint_nodes": list(range(5)),
        "definition": "M_B[i,j] = sum_v <H_v^E b_i | H_v^E b_j>",
        "master": {
            "dimension": 32,
            "rank": rank,
            "nullity": nullity,
            "rank_tolerance": rank_tol,
            "hermiticity_error": herm,
            "eigenvalue_min": float(np.min(ev)),
            "eigenvalue_max": float(np.max(ev)),
            "smallest_positive_eigenvalue": float(np.min(positive)) if len(positive) else None,
            "condition_number_on_support": float(np.max(positive) / np.min(positive)) if len(positive) else None,
            "eigenvalues": [float(x) for x in ev],
            "trace": float(np.trace(M).real),
            "frobenius_norm": float(np.linalg.norm(M)),
        },
        "per_node": node_rows,
        "node_trace_relative_spread": trace_rel_spread,
        "first_action_boundary_projection_max_norm": direct_return_max,
        "quadratic_form_regression": quadratic_checks,
        "quadratic_form_max_relative_error": quad_max_rel,
        "checks": checks,
        "interpretation": (
            "Nullity is a scientific output, not a pass target. Nullity=0 means only that the complete five-node Euclidean constraint family has no common zero vector inside the bare 32D q=2 boundary carrier. "
            "It does not exclude a zero sector in the enlarged Peter-Weyl/graph-changing habitat, whose relevant boundary observable is B^dag P0 B."
        ),
        "claim_boundary": (
            "Finite Euclidean master compressed to the frozen q=2 boundary only; no Lorentzian projector, continuum rigging map, dark sector, physical frequency or cosmological observable is established."
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
