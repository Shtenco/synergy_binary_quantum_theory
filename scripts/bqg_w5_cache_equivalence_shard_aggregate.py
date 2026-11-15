#!/usr/bin/env python3
"""Fail-closed aggregation of sharded W5 cache-on/cache-off matrix controls."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

NODES = tuple(range(5))
SEED_KS = (0, 2, 4, 6, 8)
DEFAULT_JMAX2 = (10, 12, 14)


def _load(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"{path.name}: invalid JSON: {type(exc).__name__}: {exc}") from exc


def aggregate(directory: Path, jmax2_values=DEFAULT_JMAX2):
    expected = {
        (jmax2, node, seed_k)
        for jmax2 in jmax2_values
        for node in NODES
        for seed_k in SEED_KS
    }
    found = {}
    errors = []
    paths = sorted(directory.glob("W5_CACHE_EQUIVALENCE_J*_NODE*_K*_RUN*.json"))
    for path in paths:
        report = _load(path)
        params = report.get("parameters") or {}
        try:
            nodes = params.get("nodes")
            seeds = params.get("seed_Ks")
            if not isinstance(nodes, list) or len(nodes) != 1:
                raise ValueError("expected exactly one source node")
            if not isinstance(seeds, list) or len(seeds) != 1:
                raise ValueError("expected exactly one seed K")
            key = (int(params["Jmax2"]), int(nodes[0]), int(seeds[0]))
        except Exception as exc:
            errors.append(f"{path.name}: malformed shard identity: {exc}")
            continue
        if key not in expected:
            errors.append(f"{path.name}: unexpected shard identity {key}")
            continue
        if key in found:
            errors.append(f"{path.name}: duplicate shard identity {key}")
            continue
        found[key] = path.name
        comparison = report.get("comparison") or {}
        checks = {
            "status": report.get("status") == "PASS_CACHE_EQUIVALENCE",
            "comparison_performed": report.get("comparison_performed") is True,
            "unbounded_closure": params.get("unbounded_closure") is True,
            "cache_on_matrix_complete": report.get("cache_on", {}).get("status") in {"PASS_ADJOINT_CONSISTENCY", "FAIL_ADJOINT_OR_LEAKAGE_GATE"},
            "cache_off_matrix_complete": report.get("cache_off", {}).get("status") in {"PASS_ADJOINT_CONSISTENCY", "FAIL_ADJOINT_OR_LEAKAGE_GATE"},
            "comparison_status": comparison.get("status") == "PASS_CACHE_EQUIVALENCE",
            "comparison_completed": comparison.get("comparison_completed") is True,
            "matrix_exactly_equal": comparison.get("matrix_exactly_equal") is True,
            "equivalence_passed": comparison.get("equivalence_passed") is True,
            "observed_ledgers_exactly_equal": comparison.get("observed_ledger_exactly_equal") is True,
            "static_coverage_complete": comparison.get("static_coverage_complete") is True,
        }
        failed = [name for name, ok in checks.items() if not ok]
        if failed:
            errors.append(f"{path.name}: failed checks: {', '.join(failed)}")

    missing = sorted(expected - set(found))
    if missing:
        errors.append(f"missing {len(missing)} shards")
    passed = not errors and len(found) == len(expected)
    return {
        "schema": "bqg-w5-cache-equivalence-shard-aggregate-v1",
        "status": "PASS_ALL_CACHE_SHARDS" if passed else "INCOMPLETE_OR_FAILED_CACHE_SHARDS",
        "passed": passed,
        "comparison_performed": bool(found),
        "expected_shard_count": len(expected),
        "valid_shard_count": len(found),
        "jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "seed_Ks": list(SEED_KS),
        "missing_shards": [list(key) for key in missing],
        "errors": errors,
        "shards": [
            {"Jmax2": key[0], "node": key[1], "seed_K": key[2], "artifact": name}
            for key, name in sorted(found.items())
        ],
        "scope": (
            "Every (Jmax2, source node, seed K) shard must independently produce "
            "complete unbounded cache-on and cache-off matrices, exact forward and "
            "direct-adjoint equality, and matching runtime/static ledgers. This is "
            "a finite W5 cache-equivalence control, not a continuum or self-adjointness proof."
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
        "valid_shard_count": result["valid_shard_count"],
        "expected_shard_count": result["expected_shard_count"],
        "missing_shards": len(result["missing_shards"]),
        "errors": result["errors"][:10],
    }, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
