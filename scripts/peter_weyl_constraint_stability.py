#!/usr/bin/env python3
"""Compare legacy and no-prune Peter-Weyl constraint maps on one carrier.

The map C is the direct sum over constraint nodes:
    C |a> = (H_0|a>, ..., H_{m-1}|a>).
The legacy Gram diagnostics are retained for comparison, while the original
constraint map is also analyzed by streaming Householder TSQR + direct SVD.
The output Hilbert-space matrix is never materialized densely.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import heapq
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_master_projector_pilot as PILOT
import peter_weyl_streaming_svd as STREAM_SVD


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


def spectrum_of_gram(G, gram_eigenvalue_rtol):
    """Gram-only diagnostic; its tolerance is relative to eigenvalues, not sigmas."""
    if not math.isfinite(gram_eigenvalue_rtol) or gram_eigenvalue_rtol < 0:
        raise ValueError("gram_eigenvalue_rtol must be finite and nonnegative")
    raw_hermiticity = float(np.linalg.norm(G - G.conj().T, ord="fro"))
    H = 0.5 * (G + G.conj().T)
    evals = np.linalg.eigvalsh(H)
    scale = float(np.max(np.abs(evals))) if len(evals) else 0.0
    psd_tol = 100.0 * np.finfo(float).eps * scale
    if len(evals) and float(np.min(evals)) < -psd_tol:
        raise ValueError(f"Gram matrix is not PSD: min eigenvalue {float(np.min(evals))}")
    threshold = gram_eigenvalue_rtol * scale
    clipped = np.maximum(evals, 0.0)
    singular = np.sqrt(clipped)
    positive = singular[evals > threshold]
    return {
        "gram_hermiticity_defect_frobenius": raw_hermiticity,
        "gram_min_eigenvalue": float(np.min(evals)) if len(evals) else 0.0,
        "gram_max_eigenvalue": float(np.max(evals)) if len(evals) else 0.0,
        "gram_eigenvalue_scale": scale,
        "rank_threshold_on_gram_eigenvalues": threshold,
        "rank_threshold_resolves_roundoff": bool(threshold > psd_tol),
        "estimated_rank": int(np.count_nonzero(evals > threshold)),
        "minimum_singular_value_including_nullity": float(singular[0]) if len(singular) else 0.0,
        "smallest_positive_singular_value": float(np.min(positive)) if len(positive) else None,
        "largest_singular_value": float(singular[-1]) if len(singular) else 0.0,
        "singular_values": [float(x) for x in singular],
        "psd_tolerance": psd_tol,
        "rank_tolerance_semantics": "relative Gram eigenvalue threshold; equivalent sigma threshold is sqrt(rtol)",
        "singular_value_method": "sqrt(eigvalsh(C^dagger C)); Gram formation squares condition number",
        "gram_spectral_resolution_floor": float(psd_tol),
    }


def iter_operator_row_blocks(columns, dim, block_rows=4096):
    """Yield deterministic dense row blocks from sparse columns, never full C.

    Sparse keys are (node, output-state). Sorted key references are merged
    deterministically; only one block_rows x dim dense numerical buffer is
    allocated at a time (the sparse input dictionaries remain resident).
    """
    if not isinstance(block_rows, int) or block_rows <= 0:
        raise ValueError("block_rows must be a positive integer")
    key_streams = [iter(sorted(columns[a], key=repr)) for a in range(dim)]
    merged = heapq.merge(*key_streams, key=repr)
    block = np.empty((block_rows, dim), dtype=np.complex128)
    used = 0
    previous = None
    have_previous = False
    for key in merged:
        if have_previous and key == previous:
            continue
        block[used, :] = [columns[a].get(key, 0j) for a in range(dim)]
        used += 1
        previous = key
        have_previous = True
        if used == block_rows:
            yield block.copy()
            used = 0
    if used:
        yield block[:used].copy()
    elif not have_previous:
        # An empty sparse support is the zero operator. Represent it by one
        # explicit zero row so the SVD routine can return its full nullspace.
        yield np.zeros((1, dim), dtype=np.complex128)


def direct_svd_from_columns(columns, dim, singular_rtol=1e-10, block_rows=4096):
    result = STREAM_SVD.streaming_svd(
        iter_operator_row_blocks(columns, dim, block_rows), dim, singular_rtol
    )
    return {
        "rank": result.rank,
        "singular_value_rtol": result.rank_rtol,
        "singular_values": [float(x) for x in result.singular_values],
        "right_singular_vectors": result.right_singular_vectors.tolist(),
        "row_count": result.row_count,
        "block_count": result.block_count,
        "input_scale": result.input_scale,
        "max_abs_component": result.max_abs_component,
        "method": result.method,
        "exact_rank_certified": False,
        "classification": "NUMERICAL_ONLY_NO_UPSTREAM_ERROR_BOUND",
    }


def column_digest(columns):
    records = []
    for a in sorted(columns):
        for (v, key), value in sorted(columns[a].items(), key=lambda item: repr(item[0])):
            z = complex(value)
            records.append((a, v, repr(key), z.real.hex(), z.imag.hex()))
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def run(domain="low", jmax2=5, rank_rtol=1e-10,
        singular_rtol=1e-10, block_rows=4096):
    if not math.isfinite(rank_rtol) or rank_rtol < 0:
        raise ValueError("rank_rtol must be finite and nonnegative")
    if not math.isfinite(singular_rtol) or singular_rtol < 0:
        raise ValueError("singular_rtol must be finite and nonnegative")
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
    gc.collect()

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

    # Direct singular values of C, not square roots of Gram eigenvalues.
    legacy_direct = direct_svd_from_columns(legacy, dim, singular_rtol, block_rows)
    no_prune_direct = direct_svd_from_columns(no_prune, dim, singular_rtol, block_rows)
    delta_direct = direct_svd_from_columns(delta, dim, singular_rtol, block_rows)

    legacy_spec = spectrum_of_gram(legacy_gram, rank_rtol)
    no_prune_spec = spectrum_of_gram(no_prune_gram, rank_rtol)
    delta_spec = spectrum_of_gram(delta_gram, rank_rtol)
    sigma_min = legacy_direct["singular_values"][-1]
    delta_norm = delta_direct["singular_values"][0]
    stable = bool(
        legacy_direct["rank"] == dim
        and delta_norm < sigma_min
    )
    return {
        "status": "legacy versus no-prune finite Peter-Weyl constraint-map comparison",
        "domain": domain,
        "domain_dimension": dim,
        "nodes": list(nodes),
        "Jmax": jmax2 / 2,
        "gram_eigenvalue_rtol": rank_rtol,
        "singular_value_rtol": singular_rtol,
        "streaming_block_rows": block_rows,
        "row_space_definition": "direct sum over node labels; row keys are (node, output Peter-Weyl basis key)",
        "legacy": {**legacy_spec, "direct_svd": legacy_direct,
                   "stacked_column_sha256": legacy_hash,
                   "node_column_supports": legacy_supports},
        "no_prune": {**no_prune_spec, "direct_svd": no_prune_direct,
                     "stacked_column_sha256": no_prune_hash,
                     "node_column_supports": no_prune_supports},
        "difference_no_prune_minus_legacy": {
            **delta_spec, "direct_svd": delta_direct,
            "stacked_column_sha256": delta_hash,
            "interpretation": "operator norm of C_no_prune-C_legacy",
        },
        "weyl_rank_stability_sufficient_condition": {
            "delta_operator_norm": delta_norm,
            "legacy_minimum_singular_value": sigma_min,
            "strict_inequality_holds": bool(delta_norm < sigma_min),
            "legacy_full_column_rank_at_declared_tolerance": bool(legacy_direct["rank"] == dim),
            "rank_source": "direct streaming SVD of C",
            "sufficient_rank_stability_for_computed_float_maps": stable,
            "exact_operator_rank_certified": False,
            "scope": "represented float64 maps only; upstream operator-construction and subtraction error bounds are not included",
            "note": "Weyl sufficient condition on the represented maps only; failure to satisfy it does not prove rank loss.",
        },
        "claim_boundary": "Finite carrier and finite Jmax only; float64 numerical result, not exact arithmetic or a continuum theorem.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", choices=("low", "full"), default="low")
    ap.add_argument("--jmax2", type=int, default=5)
    ap.add_argument("--rank-rtol", type=float, default=1e-10,
                    help="relative Gram-eigenvalue tolerance (not singular-value tolerance)")
    ap.add_argument("--singular-rtol", type=float, default=1e-10,
                    help="relative singular-value tolerance for direct TSQR/SVD")
    ap.add_argument("--block-rows", type=int, default=4096)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = run(args.domain, args.jmax2, args.rank_rtol,
                 args.singular_rtol, args.block_rows)
    output = args.output or Path(
        f"verification_results/PETER_WEYL_CONSTRAINT_STABILITY_{args.domain}_J{args.jmax2}.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
