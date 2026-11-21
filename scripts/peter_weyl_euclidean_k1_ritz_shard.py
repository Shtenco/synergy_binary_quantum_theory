#!/usr/bin/env python3
"""One target-node shard for the preregistered Euclidean K1 Ritz master.

For all 160 first generated columns g_(v,i)=H_v b_i, compute the common first
Gram G and one positive second-action Gram

    D^(w)_[alpha,beta] = < H_w g_alpha | H_w g_beta >.

The five node shards are summed by peter_weyl_euclidean_k1_ritz_aggregate.py.
"""
from __future__ import annotations

import argparse
import gc
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


def sparse_inner(a: dict, b: dict) -> complex:
    if len(a) <= len(b):
        return sum(np.conj(x) * b.get(k, 0.0j) for k, x in a.items())
    return np.conj(sparse_inner(b, a))


def gram(images: list[dict]) -> np.ndarray:
    n = len(images)
    out = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i, n):
            z = sparse_inner(images[i], images[j])
            out[i, j] = z
            out[j, i] = np.conj(z)
    return out


def max_spin(state: dict) -> float:
    if not state:
        return 0.0
    return max(max(key[0]) for key in state) / 2.0


def run(target_node: int, clear_every: int = 8):
    if target_node not in range(5):
        raise ValueError("target_node must be 0..4")

    basis = PW.basis_full_jhalf()
    generated = []
    labels = []
    first_supports = []
    for v in range(5):
        for i, key in enumerate(basis):
            st = PW.prune_state(PW.apply_H_cached_state({key: 1.0 + 0.0j}, v, JMAX2), PRUNE)
            generated.append(st)
            labels.append((v, i))
            first_supports.append(len(st))

    G = gram(generated)

    # Keep the small representation/intertwiner caches, but prevent the large
    # per-input T_cached dictionary from growing without bound during 160 HH
    # columns. Clearing this cache changes performance only, not amplitudes.
    PW.T_cached.cache_clear()
    second = []
    second_supports = []
    second_max_spins = []
    for alpha, st in enumerate(generated):
        out = PW.compose_on_sparse(st, target_node, JMAX2)
        second.append(out)
        second_supports.append(len(out))
        second_max_spins.append(max_spin(out))
        if clear_every > 0 and (alpha + 1) % clear_every == 0:
            PW.T_cached.cache_clear()
            gc.collect()

    D = gram(second)

    Gh = (G + G.conj().T) / 2.0
    Dh = (D + D.conj().T) / 2.0
    gev = np.linalg.eigvalsh(Gh)
    dev = np.linalg.eigvalsh(Dh)
    gscale = max(float(np.max(np.abs(gev))), 1.0)
    dscale = max(float(np.max(np.abs(dev))), 1.0)

    summary = {
        "status": "Euclidean K1 Ritz target-node shard",
        "target_node": target_node,
        "Jmax": JMAX2 / 2.0,
        "columns": len(generated),
        "first_support_min": int(min(first_supports)),
        "first_support_max": int(max(first_supports)),
        "first_support_mean": float(np.mean(first_supports)),
        "second_support_min": int(min(second_supports)),
        "second_support_max": int(max(second_supports)),
        "second_support_mean": float(np.mean(second_supports)),
        "second_max_spin": float(max(second_max_spins)),
        "G_hermiticity_error": float(np.linalg.norm(G - G.conj().T)),
        "D_hermiticity_error": float(np.linalg.norm(D - D.conj().T)),
        "G_min_eigenvalue": float(np.min(gev)),
        "D_min_eigenvalue": float(np.min(dev)),
        "G_PSD": bool(float(np.min(gev)) > -3e-9 * gscale),
        "D_PSD": bool(float(np.min(dev)) > -3e-9 * dscale),
        "D_trace": float(np.trace(D).real),
        "D_frobenius_norm": float(np.linalg.norm(D)),
        "cache_clear_every": clear_every,
        "claim_boundary": "One Euclidean HH target-node Gram only; no physical projector or Lorentzian result.",
    }
    summary["passed"] = bool(
        summary["columns"] == 160
        and summary["G_hermiticity_error"] < 3e-9
        and summary["D_hermiticity_error"] < 3e-9
        and summary["G_PSD"]
        and summary["D_PSD"]
        and summary["second_max_spin"] <= 2.5 + 1e-12
    )
    return G, D, np.asarray(labels, dtype=int), summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-node", type=int, required=True)
    ap.add_argument("--output-npz", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--clear-cache-every", type=int, default=8)
    args = ap.parse_args()

    G, D, labels, summary = run(args.target_node, args.clear_cache_every)
    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        target_node=np.array([args.target_node], dtype=int),
        G=G,
        D=D,
        labels=labels,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
