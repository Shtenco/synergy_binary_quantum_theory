#!/usr/bin/env python3
"""Compare legacy and no-prune Peter-Weyl constraint maps on one carrier.

The map C is the direct sum over constraint nodes:
    C |a> = (H_0|a>, ..., H_{m-1}|a>).
All spectral comparisons are made through small dim-by-dim Gram matrices,
never by materializing the potentially enormous output Hilbert-space matrix.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_master_projector_pilot as PILOT


def stack_columns(columns, nodes, dim):
    """Combine node outputs into a direct-sum output space for each input."""
    stacked = {}
    supports = {}
    for a in range(dim):
        out = {}
        per_node = {}
        for v in nodes:
            col = columns[(v, a)]
            per_node[str(v)] = len(col)
            for key, value in col.items():
                out[(v, key)] = complex(value)
        stacked[a] = out
        supports[str(a)] = per_node
    return stacked, supports


def subtract_stacked(left, right, dim):
    """Return right-left, retaining every exactly nonzero difference."""
    delta = {}
    for a in range(dim):
        keys = set(left[a]) | set(right[a])
        col = {}
        for key in keys:
            value = right[a].get(key, 0j) - left[a].get(key, 0j)
            if value != 0j:
                col[key] = value
        delta[a] = col
    return delta


def gram_from_columns(columns, dim):
    G = np.empty((dim, dim), dtype=complex)
    for a in range(dim):
        for b in range(dim):
            G[a, b] = PILOT.inner_sparse(columns[a], columns[b])
    return G


def spectrum_of_gram(G, rank_rtol):
    raw_hermiticity = float(np.linalg.norm(G - G.conj().T, ord="fro"))
    H = 0.5 * (G + G.conj().T)
    evals = np.linalg.eigvalsh(H)
    scale = max(1.0, float(np.max(np.abs(evals))))
    psd_tol = 1e-10 * scale
    if float(np.min(evals)) < -psd_tol:
        raise ValueError(f"Gram matrix is not PSD: min eigenvalue {float(np.min(evals))}")
    threshold = rank_rtol * scale
    clipped = np.maximum(evals, 0.0)
    singular = np.sqrt(clipped)
    positive = singular[evals > threshold]
    return {
        "gram_hermiticity_defect_frobenius": raw_hermiticity,
        "gram_min_eigenvalue": float(np.min(evals)),
        "gram_max_eigenvalue": float(np.max(evals)),
        "rank_threshold_on_gram_eigenvalues": threshold,
        "estimated_rank": int(np.count_nonzero(evals > threshold)),
        "minimum_singular_value_including_nullity": float(singular[0]),
        "smallest_positive_singular_value": float(np.min(positive)) if len(positive) else None,
        "largest_singular_value": float(singular[-1]),
        "singular_values": [float(x) for x in singular],
        "psd_tolerance": psd_tol,
    }


def column_digest(columns):
    records = []
    for a in sorted(columns):
        for (v, key), value in sorted(columns[a].items(), key=lambda item: repr(item[0])):
            z = complex(value)
            records.append((a, v, repr(key), z.real.hex(), z.imag.hex()))
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def run(domain="low", jmax2=5, rank_rtol=1e-10):
    if not math.isfinite(rank_rtol) or rank_rtol < 0:
        raise ValueError("rank_rtol must be finite and nonnegative")
    nodes = tuple(PILOT.PW.VERT)
    basis = PILOT.carrier_basis(domain)
    dim = len(basis)

    legacy_raw, legacy_supports = PILOT.constraint_columns(
        basis, nodes, jmax2, mode="legacy"
    )
    legacy, legacy_node_supports = stack_columns(legacy_raw, nodes, dim)
    legacy_hash = column_digest(legacy)
    legacy_gram = gram_from_columns(legacy, dim)
    del legacy_raw, legacy_node_supports

    # The producer cache is keyed by pruning mode; clear it between modes so
    # the second calculation cannot accidentally reuse legacy cache entries.
    if hasattr(PILOT.PW.T_cached, "cache_clear"):
        PILOT.PW.T_cached.cache_clear()

    no_prune_raw, no_prune_supports = PILOT.constraint_columns(
        basis, nodes, jmax2, mode="no_prune"
    )
    no_prune, no_prune_node_supports = stack_columns(no_prune_raw, nodes, dim)
    no_prune_hash = column_digest(no_prune)
    no_prune_gram = gram_from_columns(no_prune, dim)
    del no_prune_raw, no_prune_node_supports

    delta = subtract_stacked(legacy, no_prune, dim)
    delta_hash = column_digest(delta)
    delta_gram = gram_from_columns(delta, dim)

    legacy_spec = spectrum_of_gram(legacy_gram, rank_rtol)
    no_prune_spec = spectrum_of_gram(no_prune_gram, rank_rtol)
    delta_spec = spectrum_of_gram(delta_gram, rank_rtol)
    sigma_min = legacy_spec["minimum_singular_value_including_nullity"]
    delta_norm = delta_spec["largest_singular_value"]
    stable = bool(
        legacy_spec["estimated_rank"] == dim
        and delta_norm < sigma_min
    )
    return {
        "status": "legacy versus no-prune finite Peter-Weyl constraint-map comparison",
        "domain": domain,
        "domain_dimension": dim,
        "nodes": list(nodes),
        "Jmax": jmax2 / 2,
        "rank_rtol": rank_rtol,
        "row_space_definition": "direct sum over node labels; row keys are (node, output Peter-Weyl basis key)",
        "legacy": {**legacy_spec, "stacked_column_sha256": legacy_hash,
                   "node_column_supports": legacy_supports},
        "no_prune": {**no_prune_spec, "stacked_column_sha256": no_prune_hash,
                     "node_column_supports": no_prune_supports},
        "difference_no_prune_minus_legacy": {
            **delta_spec, "stacked_column_sha256": delta_hash,
            "interpretation": "operator norm of C_no_prune-C_legacy",
        },
        "weyl_rank_stability_sufficient_condition": {
            "delta_operator_norm": delta_norm,
            "legacy_minimum_singular_value": sigma_min,
            "strict_inequality_holds": bool(delta_norm < sigma_min),
            "legacy_full_column_rank_at_declared_tolerance": bool(legacy_spec["estimated_rank"] == dim),
            "sufficient_rank_stability_certified": stable,
            "note": "Sufficient perturbation bound only; failure to satisfy it does not prove rank loss.",
        },
        "claim_boundary": "Finite carrier and finite Jmax only; float64 numerical result, not exact arithmetic or a continuum theorem.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", choices=("low", "full"), default="low")
    ap.add_argument("--jmax2", type=int, default=5)
    ap.add_argument("--rank-rtol", type=float, default=1e-10)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = run(args.domain, args.jmax2, args.rank_rtol)
    output = args.output or Path(
        f"verification_results/PETER_WEYL_CONSTRAINT_STABILITY_{args.domain}_J{args.jmax2}.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
