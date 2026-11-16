#!/usr/bin/env python3
"""Aggregate channel-sharded W5 production vs shadow_no_prune sparse outputs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

MODES = ("production", "shadow_no_prune")
NODES = tuple(range(5))
CHANNELS = ((0, 2), (2, 0))
M2_VALUES = (-2, 0, 2)


def _norm(values):
    result = 0.0
    for value in values:
        result = math.hypot(result, float(value))
        if not math.isfinite(result):
            raise ValueError("non-finite sparse-vector norm")
    return result


def _decode_entries(column, *, context):
    entries = column.get("entries")
    if not isinstance(entries, list):
        raise ValueError(f"{context}: missing sparse entries")
    result = {}
    for row in entries:
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError(f"{context}: malformed sparse entry")
        key, real, imag = row
        if not isinstance(key, str):
            raise ValueError(f"{context}: output key must be a string")
        real, imag = float(real), float(imag)
        if not math.isfinite(real) or not math.isfinite(imag):
            raise ValueError(f"{context}: non-finite output amplitude")
        if key in result:
            raise ValueError(f"{context}: duplicate output key {key!r}")
        if real != 0.0 or imag != 0.0:
            result[key] = complex(real, imag)
    return result


def _load_report(directory: Path, mode: str, node: int, channel):
    k12, k34 = channel
    path = directory / f"direct_all_M_{mode}_node{node}_channel{k12}_{k34}.json"
    if not path.is_file():
        raise FileNotFoundError(f"missing required artifact: {path.name}")
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("schema") != "bqg-w5-direct-all-M-su2-commutator-v1":
        raise ValueError(f"{path.name}: unexpected schema")
    if report.get("status") != "PASS_DIRECT_ALL_M_SU2_COMMUTATOR" or report.get("passed") is not True:
        raise ValueError(f"{path.name}: underlying SU(2) gate did not pass")
    if report.get("pruning_mode") != mode:
        raise ValueError(f"{path.name}: pruning_mode mismatch")
    if report.get("source_nodes") != [node]:
        raise ValueError(f"{path.name}: source-node mismatch")
    if report.get("Jmax2") != 10:
        raise ValueError(f"{path.name}: unexpected Jmax2")
    declared = report.get("input_charged_sector", {}).get("channels_K12_K34")
    if declared != [[k12, k34]]:
        raise ValueError(f"{path.name}: channel provenance mismatch")
    columns = report.get("operator_columns_sparse")
    if not isinstance(columns, list) or len(columns) != 3:
        raise ValueError(f"{path.name}: expected three complete M columns for this channel")
    decoded = {}
    for column in columns:
        identity = (
            int(column["source_node"]), int(column["K12"]),
            int(column["K34"]), int(column["input_M2"]),
        )
        if identity in decoded:
            raise ValueError(f"{path.name}: duplicate input column {identity}")
        if identity[:3] != (node, k12, k34) or identity[3] not in M2_VALUES:
            raise ValueError(f"{path.name}: unexpected input column {identity}")
        decoded[identity] = _decode_entries(column, context=f"{path.name} {identity}")
    expected = {(node, k12, k34, m2) for m2 in M2_VALUES}
    if set(decoded) != expected:
        raise ValueError(f"{path.name}: incomplete or unexpected input-column set")
    return report, decoded


def compare_column(reference, candidate):
    keys = set(reference) | set(candidate)
    differences, reference_values = [], []
    max_absolute = 0.0
    for key in keys:
        left, right = reference.get(key, 0j), candidate.get(key, 0j)
        magnitude = abs(left - right)
        differences.append(magnitude)
        reference_values.append(abs(left))
        max_absolute = max(max_absolute, magnitude)
    difference_norm = _norm(differences)
    reference_norm = _norm(reference_values)
    return {
        "reference_support": len(reference),
        "candidate_support": len(candidate),
        "support_symmetric_difference": len(set(reference) ^ set(candidate)),
        "difference_l2_norm": difference_norm,
        "reference_l2_norm": reference_norm,
        "max_absolute_difference": max_absolute,
        "relative_l2_error": difference_norm / max(reference_norm, 1e-300),
    }


def compare_directory(directory: Path, *, tolerance=1e-10):
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be finite and nonnegative")
    reports, decoded = {}, {}
    for node in NODES:
        for channel in CHANNELS:
            for mode in MODES:
                key = (mode, node, channel)
                reports[key], decoded[key] = _load_report(directory, mode, node, channel)

    pairs, passed = [], True
    for node in NODES:
        for channel in CHANNELS:
            production = decoded[("production", node, channel)]
            shadow = decoded[("shadow_no_prune", node, channel)]
            for identity in sorted(production):
                metrics = compare_column(production[identity], shadow[identity])
                row_passed = metrics["relative_l2_error"] <= tolerance
                passed = passed and row_passed
                pairs.append({
                    "source_node": node, "K12": channel[0], "K34": channel[1],
                    "input_M2": identity[3], "passed": row_passed, **metrics,
                })
    return {
        "schema": "bqg-w5-direct-all-m-shadow-equivalence-v2",
        "status": "PASS_PRODUCTION_SHADOW_EQUIVALENCE" if passed else "FAIL_PRODUCTION_SHADOW_EQUIVALENCE",
        "passed": passed,
        "comparison_performed": True,
        "nodes": list(NODES),
        "channels": [list(x) for x in CHANNELS],
        "pruning_modes": list(MODES),
        "Jmax2": 10,
        "input_columns_compared": len(pairs),
        "expected_input_columns": len(NODES) * len(CHANNELS) * len(M2_VALUES),
        "relative_l2_tolerance": tolerance,
        "max_relative_l2_error": max((row["relative_l2_error"] for row in pairs), default=0.0),
        "summed_production_case_pruning_event_observations": sum(
            int(row.get("pruning_events", 0))
            for (mode, _, _), report in reports.items() if mode == "production"
            for row in report.get("cases", [])
        ),
        "summed_shadow_case_pruning_event_observations": sum(
            int(row.get("pruning_events", 0))
            for (mode, _, _), report in reports.items() if mode == "shadow_no_prune"
            for row in report.get("cases", [])
        ),
        "pairs": pairs,
        "scope": (
            "Compares complete sparse H_E output maps for the declared J=1 "
            "recoupling channels and all M at all five K5 source nodes. "
            "This is a finite seed/cutoff control, not a global pruning theorem, "
            "all-spin SU(2) proof, or continuum self-adjointness result."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare_directory(args.directory, tolerance=args.tolerance)
    except Exception as exc:
        result = {
            "schema": "bqg-w5-direct-all-m-shadow-equivalence-v2",
            "status": "INCOMPLETE_PRODUCTION_SHADOW_EQUIVALENCE",
            "passed": False,
            "comparison_performed": False,
            "error": f"{type(exc).__name__}: {exc}",
            "scope": "No equivalence claim is possible without all twenty valid channel/node/mode artifacts.",
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "comparison_performed": result["comparison_performed"],
        "input_columns_compared": result.get("input_columns_compared", 0),
        "max_relative_l2_error": result.get("max_relative_l2_error"),
        "output": str(args.output),
    }, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
