#!/usr/bin/env python3
"""Pair independently executed W5 cache-on/cache-off shards fail-closed."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from bqg_v10_w5_amp_tol_shadow_compare import compare_cache_runs, closure_matrix_available

NODES = tuple(range(5))
SEED_KS = (0, 2, 4, 6, 8)
DEFAULT_JMAX2 = (10, 12, 14)
MODES = ("on", "off")


def aggregate(directory: Path, jmax2_values=DEFAULT_JMAX2):
    expected = {
        (jmax2, node, seed_k)
        for jmax2 in jmax2_values
        for node in NODES
        for seed_k in SEED_KS
    }
    found = {}
    errors = []
    for path in sorted(directory.glob("W5_CACHE_MATRIX_J*_NODE*_K*_MODE*_RUN*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
            nodes = report.get("source_nodes")
            seeds = report.get("seed_Ks")
            if not isinstance(nodes, list) or len(nodes) != 1:
                raise ValueError("expected one source node")
            if not isinstance(seeds, list) or len(seeds) != 1:
                raise ValueError("expected one seed K")
            key = (int(report["Jmax2"]), int(nodes[0]), int(seeds[0]))
            mode = report.get("cache_mode")
            if key not in expected or mode not in MODES:
                errors.append(f"{path.name}: unexpected shard identity {(key, mode)}")
                continue
            slot = (key, mode)
            if slot in found:
                errors.append(f"{path.name}: duplicate shard identity {slot}")
                continue
            found[slot] = report
            checks = {
                "status": report.get("status") == "PASS_ADJOINT_CONSISTENCY",
                "production_pruning_mode": report.get("pruning_mode") == "production",
                "unbounded_state_closure": report.get("max_states") is None,
                "unbounded_depth_closure": report.get("max_depth") is None,
                "complete_matrix": closure_matrix_available(report),
                "basis_keys": isinstance(report.get("basis_keys_repr"), list),
                "forward_columns": isinstance(report.get("forward_B_columns_sparse"), list),
                "direct_adjoint_columns": isinstance(report.get("direct_Bdag_columns_sparse"), list),
                "lower_level_ledger": isinstance(report.get("lower_level_pruning_ledger"), dict),
                "coefficient_loss_ledger": isinstance(report.get("coefficient_loss_ledger"), dict),
                "coherent_node_sum_ledger": isinstance(report.get("coherent_node_sum_ledger"), list),
                "spectral_sensitivity_ledger": isinstance(report.get("spectral_sensitivity"), dict),
                "diagnostics": isinstance(report.get("diagnostics"), dict),
            }
            failed = [name for name, ok in checks.items() if not ok]
            if failed:
                errors.append(f"{path.name}: failed checks: {', '.join(failed)}")
        except Exception as exc:
            errors.append(f"{path.name}: malformed report: {type(exc).__name__}: {exc}")

    comparisons = []
    for key in sorted(expected):
        on = found.get((key, "on"))
        off = found.get((key, "off"))
        if on is None or off is None:
            continue
        comparison = compare_cache_runs(on, off)
        passed = comparison.get("status") == "PASS_CACHE_EQUIVALENCE"
        if not passed:
            errors.append(f"J{key[0]} node{key[1]} seedK{key[2]}: {comparison.get('status')}")
        comparisons.append({
            "Jmax2": key[0],
            "node": key[1],
            "seed_K": key[2],
            "cache_on_status": on.get("status"),
            "cache_off_status": off.get("status"),
            "comparison": comparison,
            "passed": passed,
        })

    expected_count = len(expected) * len(MODES)
    missing = [
        [key[0], key[1], key[2], mode]
        for key in sorted(expected)
        for mode in MODES
        if (key, mode) not in found
    ]
    if missing:
        errors.append(f"missing {len(missing)} cache-mode shards")
    passed = not errors and len(found) == expected_count and len(comparisons) == len(expected)
    return {
        "schema": "bqg-w5-cache-equivalence-paired-shard-aggregate-v2",
        "status": "PASS_ALL_CACHE_PAIRS" if passed else "INCOMPLETE_OR_FAILED_CACHE_PAIRS",
        "passed": passed,
        "expected_pair_count": len(expected),
        "completed_pair_count": len(comparisons),
        "expected_mode_shard_count": expected_count,
        "valid_mode_shard_count": len(found),
        "Jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "seed_Ks": list(SEED_KS),
        "missing_mode_shards": missing,
        "errors": errors,
        "comparisons": comparisons,
        "scope": (
            "Every cache-on and cache-off run is a separate fresh process and CI job. "
            "Each pair is compared using exact full complex forward/direct-adjoint matrices "
            "and all required ledgers. PASS applies only to the declared finite W5 cutoffs, "
            "nodes and seed-K set."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--Jmax2", type=int, nargs="+", default=list(DEFAULT_JMAX2))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.directory, args.Jmax2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "completed_pair_count": result["completed_pair_count"],
        "expected_pair_count": result["expected_pair_count"],
        "missing_mode_shards": len(result["missing_mode_shards"]),
        "errors": result["errors"][:12],
    }, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
