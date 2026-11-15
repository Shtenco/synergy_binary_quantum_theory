#!/usr/bin/env python3
"""Pair W5 production and shadow-no-prune matrices without overstating closure."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available,
    compare_vectors,
    unpack_columns,
)

NODES = tuple(range(5))
DEFAULT_JMAX2 = (10, 12, 14)
MODES = ("production", "shadow_no_prune")


def aggregate(directory: Path, jmax2_values=DEFAULT_JMAX2):
    expected = {(jmax2, node) for jmax2 in jmax2_values for node in NODES}
    found = {}
    errors = []
    for path in sorted(directory.glob("W5_SHADOW_MATRIX_J*_NODE*_MODE*_RUN*.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
            nodes = report.get("source_nodes")
            if not isinstance(nodes, list) or len(nodes) != 1:
                raise ValueError("expected exactly one source node")
            key = (int(report["Jmax2"]), int(nodes[0]))
            mode = report.get("pruning_mode")
            if key not in expected or mode not in MODES:
                errors.append(f"{path.name}: unexpected shard identity {(key, mode)}")
                continue
            slot = (key, mode)
            if slot in found:
                errors.append(f"{path.name}: duplicate shard identity {slot}")
                continue
            found[slot] = report
            matrix_complete = closure_matrix_available(
                report, shadow=(mode == "shadow_no_prune")
            )
            checks = {
                "complete_finite_matrix": matrix_complete,
                "unbounded_state_closure": report.get("max_states") is None,
                "unbounded_depth_closure": report.get("max_depth") is None,
                "all_five_seed_Ks": report.get("seed_Ks") == [0, 2, 4, 6, 8],
                "cache_on_for_both_modes": report.get("cache_mode") == "on",
                "forward_columns": isinstance(report.get("forward_B_columns_sparse"), list),
                "direct_adjoint_columns": isinstance(report.get("direct_Bdag_columns_sparse"), list),
                "coefficient_loss_ledger": isinstance(report.get("coefficient_loss_ledger"), dict),
                "lower_level_pruning_ledger": isinstance(report.get("lower_level_pruning_ledger"), dict),
                "coherent_node_sum_ledger": isinstance(report.get("coherent_node_sum_ledger"), list),
            }
            if mode == "production":
                checks["production_operator_gate"] = report.get("status") == "PASS_ADJOINT_CONSISTENCY"
            else:
                checks["shadow_operator_gate"] = (
                    report.get("status") == "INCOMPLETE_NO_PRUNING_COVERAGE"
                    and report.get("algebraic_status") == "PASS_ADJOINT_CONSISTENCY"
                )
            failed = [name for name, ok in checks.items() if not ok]
            if failed:
                errors.append(f"{path.name}: failed checks: {', '.join(failed)}")
        except Exception as exc:
            errors.append(f"{path.name}: malformed report: {type(exc).__name__}: {exc}")

    comparisons = []
    for key in sorted(expected):
        production = found.get((key, "production"))
        shadow = found.get((key, "shadow_no_prune"))
        if production is None or shadow is None:
            continue
        forward = compare_vectors(unpack_columns(production), unpack_columns(shadow))
        adjoint = compare_vectors(
            unpack_columns(production, "direct_Bdag_columns_sparse"),
            unpack_columns(shadow, "direct_Bdag_columns_sparse"),
        )
        finite = all(
            isinstance(item.get("relative_l2_residual"), (int, float))
            and math.isfinite(float(item["relative_l2_residual"]))
            for item in (forward, adjoint)
        )
        if not finite:
            errors.append(f"J{key[0]} node{key[1]}: non-finite matrix comparison")
        comparisons.append({
            "Jmax2": key[0],
            "node": key[1],
            "production_status": production.get("status"),
            "shadow_status": shadow.get("status"),
            "forward_matrix_comparison": forward,
            "direct_adjoint_matrix_comparison": adjoint,
            "production_pruning_event_count": (production.get("lower_level_pruning_ledger") or {}).get("event_count"),
            "shadow_pruning_event_count": (shadow.get("lower_level_pruning_ledger") or {}).get("event_count"),
            "production_coefficient_loss_ledger": production.get("coefficient_loss_ledger"),
            "shadow_coefficient_loss_ledger": shadow.get("coefficient_loss_ledger"),
            "production_coherent_node_sum_ledger": production.get("coherent_node_sum_ledger"),
            "shadow_coherent_node_sum_ledger": shadow.get("coherent_node_sum_ledger"),
        })

    expected_modes = len(expected) * len(MODES)
    missing = [
        [key[0], key[1], mode]
        for key in sorted(expected)
        for mode in MODES
        if (key, mode) not in found
    ]
    if missing:
        errors.append(f"missing {len(missing)} production/shadow mode shards")
    complete = not errors and len(found) == expected_modes and len(comparisons) == len(expected)
    return {
        "schema": "bqg-w5-shadow-no-prune-paired-shard-aggregate-v2",
        "status": (
            "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED"
            if complete else "INCOMPLETE_OR_FAILED_SHADOW_SHARDS"
        ),
        "passed": False,
        "finite_matrix_comparison_complete": complete,
        "expected_pair_count": len(expected),
        "completed_pair_count": len(comparisons),
        "expected_mode_shard_count": expected_modes,
        "valid_mode_shard_count": len(found),
        "Jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "missing_mode_shards": missing,
        "errors": errors,
        "comparisons": comparisons,
        "claim_boundary": (
            "This paired finite matrix diagnostic never marks the global pruning proof "
            "complete. Uninstrumented thresholds, spin-wall/state-depth regulator closure, "
            "spectral conventions and continuum convergence remain open."
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
        "finite_matrix_comparison_complete": result["finite_matrix_comparison_complete"],
        "completed_pair_count": result["completed_pair_count"],
        "expected_pair_count": result["expected_pair_count"],
        "missing_mode_shards": len(result["missing_mode_shards"]),
        "errors": result["errors"][:12],
    }, indent=2))
    return 0 if result["finite_matrix_comparison_complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
