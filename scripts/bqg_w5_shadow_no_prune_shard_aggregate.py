#!/usr/bin/env python3
"""Aggregate finite W5 production-vs-shadow comparisons by node and cutoff.

The aggregate deliberately remains incomplete as a global pruning proof:
instrumented thresholds are bypassed, but uninstrumented sites and other
regulator boundaries are still outside the claim.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

NODES = tuple(range(5))
JMAX2_VALUES = (10, 12, 14)


def aggregate(directory: Path, jmax2_values=JMAX2_VALUES):
    expected = {(j, n) for j in jmax2_values for n in NODES}
    found = {}
    errors = []
    for path in sorted(directory.glob("W5_FULL_SHADOW_J*_NODE*_RUN*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
            params = report.get("parameters") or {}
            nodes = report.get("source_nodes")
            jmax2 = int(report["Jmax2"])
            if not isinstance(nodes, list) or len(nodes) != 1:
                raise ValueError("expected one source node")
            node = int(nodes[0])
            key = (jmax2, node)
        except Exception as exc:
            errors.append(f"{path.name}: malformed report: {type(exc).__name__}: {exc}")
            continue
        if key not in expected:
            errors.append(f"{path.name}: unexpected shard {key}")
            continue
        if key in found:
            errors.append(f"{path.name}: duplicate shard {key}")
            continue
        found[key] = path.name
        checks = {
            "expected_finite_status": report.get("status") == "INCOMPLETE_NO_PRUNING_REFERENCE",
            "comparison_performed": report.get("comparison_performed") is True,
            "unbounded_closure": params.get("unbounded_closure") is True,
            "instrumented_thresholds_bypassed": report.get("lower_level_pruning_bypassed_at_instrumented_sites") is True,
            "global_bypass_not_claimed": report.get("all_lower_level_pruning_bypassed") is False,
            "unverified_sites_disclosed": report.get("unverified_lower_level_sites_may_remain") is True,
            "production_loss_ledger": isinstance(report.get("production_coefficient_loss_ledger"), dict),
            "shadow_loss_ledger": isinstance(report.get("shadow_coefficient_loss_ledger"), dict),
            "production_coherent_ledger": isinstance(report.get("production_coherent_node_sum_ledger"), (dict, list)),
            "shadow_coherent_ledger": isinstance(report.get("shadow_coherent_node_sum_ledger"), (dict, list)),
        }
        for field in ("forward_matrix_comparison", "direct_adjoint_matrix_comparison"):
            comparison = report.get(field) or {}
            residual = comparison.get("relative_l2_residual")
            checks[field + "_finite"] = (
                isinstance(residual, (int, float)) and math.isfinite(float(residual))
            )
        failed = [name for name, ok in checks.items() if not ok]
        if failed:
            errors.append(f"{path.name}: failed checks: {', '.join(failed)}")

    missing = sorted(expected - set(found))
    if missing:
        errors.append(f"missing {len(missing)} node/cutoff shards")
    complete = not errors and len(found) == len(expected)
    return {
        "schema": "bqg-w5-shadow-no-prune-shard-aggregate-v1",
        "status": (
            "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED"
            if complete else "INCOMPLETE_OR_FAILED_SHADOW_SHARDS"
        ),
        "passed": False,
        "finite_matrix_comparison_complete": complete,
        "comparison_performed": bool(found),
        "expected_shard_count": len(expected),
        "valid_shard_count": len(found),
        "Jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "missing_shards": [list(key) for key in missing],
        "errors": errors,
        "shards": [
            {"Jmax2": key[0], "node": key[1], "artifact": value}
            for key, value in sorted(found.items())
        ],
        "claim_boundary": (
            "Even complete finite matrices do not close the global pruning proof. "
            "Uninstrumented threshold sites, spin-wall/state-depth caps, spectral "
            "conventions and continuum convergence remain open."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--Jmax2", type=int, nargs="+", default=list(JMAX2_VALUES))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.directory, args.Jmax2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "finite_matrix_comparison_complete": result["finite_matrix_comparison_complete"],
        "valid_shard_count": result["valid_shard_count"],
        "expected_shard_count": result["expected_shard_count"],
        "errors": result["errors"][:10],
    }, indent=2))
    return 0 if result["finite_matrix_comparison_complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
