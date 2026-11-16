#!/usr/bin/env python3
"""Aggregate adjacent-M-sharded W5 production vs shadow_no_prune sparse outputs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

MODES = ("production", "shadow_no_prune")
NODES = tuple(range(5))
CHANNELS = ((0, 2), (2, 0))
M_PAIRS = ((-2, 0), (0, 2))
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


def _pair_tag(pair):
    return "-2_0" if pair == (-2, 0) else "0_2"


def _load_pair(directory: Path, mode: str, node: int, channel, pair):
    k12, k34 = channel
    tag = _pair_tag(pair)
    path = directory / (
        f"direct_all_M_{mode}_node{node}_channel{k12}_{k34}_mpair{tag}.json"
    )
    if not path.is_file():
        raise FileNotFoundError(f"missing required artifact: {path.name}")
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("schema") != "bqg-w5-direct-all-M-su2-commutator-v1":
        raise ValueError(f"{path.name}: unexpected schema")
    if report.get("status") != "PASS_DIRECT_ALL_M_SU2_COMMUTATOR" or report.get("passed") is not True:
        raise ValueError(f"{path.name}: underlying adjacent-M SU(2) gate did not pass")
    if report.get("pruning_mode") != mode:
        raise ValueError(f"{path.name}: pruning_mode mismatch")
    if report.get("source_nodes") != [node]:
        raise ValueError(f"{path.name}: source-node mismatch")
    if report.get("Jmax2") != 10:
        raise ValueError(f"{path.name}: unexpected Jmax2")
    declared = report.get("input_charged_sector", {}).get("channels_K12_K34")
    if declared != [[k12, k34]]:
        raise ValueError(f"{path.name}: channel provenance mismatch")
    if report.get("magnetic_components_tested_M2") != list(pair):
        raise ValueError(f"{path.name}: adjacent-M pair provenance mismatch")
    expected_cases = {("J_plus", pair[0]), ("J_minus", pair[1])}
    cases = report.get("cases")
    if not isinstance(cases, list) or {
        (row.get("direction"), int(row.get("input_M2", 999))) for row in cases
    } != expected_cases:
        raise ValueError(f"{path.name}: both adjacent ladder directions are required")
    if report.get("all_full_output_matrices_nonzero") is not True:
        raise ValueError(f"{path.name}: full sparse operator output is not certified nonzero")
    numeric_gates = (
        ("max_commutator_relative_residual", 1e-8),
        ("max_jz_commutator_relative_residual", 1e-12),
        ("max_wrong_JM_weight_fraction", 1e-18),
        ("max_spectator_non_singlet_weight_fraction", 1e-18),
        ("max_actual_complex_superposition_relative_residual", 1e-10),
    )
    for field, limit in numeric_gates:
        value = float(report.get(field, float("inf")))
        if not math.isfinite(value) or value >= limit:
            raise ValueError(f"{path.name}: failed or missing {field}")
    columns = report.get("operator_columns_sparse")
    if not isinstance(columns, list) or len(columns) != 2:
        raise ValueError(f"{path.name}: expected exactly two sparse M columns")
    decoded = {}
    for column in columns:
        identity = (
            int(column["source_node"]), int(column["K12"]),
            int(column["K34"]), int(column["input_M2"]),
        )
        if identity in decoded:
            raise ValueError(f"{path.name}: duplicate input column {identity}")
        if identity[:3] != (node, k12, k34) or identity[3] not in pair:
            raise ValueError(f"{path.name}: unexpected input column {identity}")
        decoded[identity] = _decode_entries(column, context=f"{path.name} {identity}")
        if not decoded[identity]:
            raise ValueError(f"{path.name}: zero/empty sparse output column {identity}")
    expected = {(node, k12, k34, m2) for m2 in pair}
    if set(decoded) != expected:
        raise ValueError(f"{path.name}: incomplete adjacent-M column set")
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


def compare_directory(directory: Path, *, tolerance=1e-10, duplicate_tolerance=1e-12):
    if (not math.isfinite(tolerance) or tolerance < 0.0
            or not math.isfinite(duplicate_tolerance) or duplicate_tolerance < 0.0):
        raise ValueError("tolerances must be finite and nonnegative")
    reports, decoded = {}, {}
    for node in NODES:
        for channel in CHANNELS:
            for mode in MODES:
                merged = {}
                for pair in M_PAIRS:
                    key = (mode, node, channel, pair)
                    report, columns = _load_pair(directory, mode, node, channel, pair)
                    reports[key] = report
                    for identity, values in columns.items():
                        m2 = identity[3]
                        if m2 in merged:
                            duplicate = compare_column(merged[m2], values)
                            if duplicate["relative_l2_error"] > duplicate_tolerance:
                                raise ValueError(
                                    f"overlapping M2={m2} columns disagree between adjacent-M shards "
                                    f"for mode={mode}, node={node}, channel={channel}"
                                )
                        else:
                            merged[m2] = values
                if set(merged) != set(M2_VALUES):
                    raise ValueError(f"incomplete merged all-M columns for {mode}, node={node}, channel={channel}")
                decoded[(mode, node, channel)] = merged

    pairs, passed = [], True
    for node in NODES:
        for channel in CHANNELS:
            production = decoded[("production", node, channel)]
            shadow = decoded[("shadow_no_prune", node, channel)]
            for m2 in M2_VALUES:
                metrics = compare_column(production[m2], shadow[m2])
                row_passed = metrics["relative_l2_error"] <= tolerance
                passed = passed and row_passed
                pairs.append({
                    "source_node": node, "K12": channel[0], "K34": channel[1],
                    "input_M2": m2, "passed": row_passed, **metrics,
                })
    production_reports = [r for (mode, _, _, _), r in reports.items() if mode == "production"]
    shadow_reports = [r for (mode, _, _, _), r in reports.items() if mode == "shadow_no_prune"]
    return {
        "schema": "bqg-w5-direct-all-m-shadow-equivalence-v3-adjacent-m-shards",
        "status": "PASS_PRODUCTION_SHADOW_EQUIVALENCE" if passed else "FAIL_PRODUCTION_SHADOW_EQUIVALENCE",
        "passed": passed,
        "comparison_performed": True,
        "nodes": list(NODES),
        "channels": [list(x) for x in CHANNELS],
        "adjacent_M_pairs": [list(x) for x in M_PAIRS],
        "pruning_modes": list(MODES),
        "Jmax2": 10,
        "input_columns_compared": len(pairs),
        "expected_input_columns": len(NODES) * len(CHANNELS) * len(M2_VALUES),
        "adjacent_M_overlap_consistency_tolerance": duplicate_tolerance,
        "relative_l2_tolerance": tolerance,
        "max_relative_l2_error": max((row["relative_l2_error"] for row in pairs), default=0.0),
        "summed_production_case_pruning_event_observations": sum(
            int(row.get("pruning_events", 0)) for r in production_reports for row in r.get("cases", [])
        ),
        "summed_shadow_case_pruning_event_observations": sum(
            int(row.get("pruning_events", 0)) for r in shadow_reports for row in r.get("cases", [])
        ),
        "pairs": pairs,
        "scope": (
            "Compares complete sparse H_E output maps for the declared J=1 "
            "recoupling channels and all M at all five K5 source nodes. Each "
            "adjacent-M shard checks both ladder directions for its pair; the "
            "shared M=0 column is cross-checked between shards. This is a finite "
            "seed/cutoff control, not a global pruning theorem, all-spin SU(2) "
            "proof, or continuum self-adjointness result."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--duplicate-tolerance", type=float, default=1e-12)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare_directory(
            args.directory, tolerance=args.tolerance,
            duplicate_tolerance=args.duplicate_tolerance,
        )
    except Exception as exc:
        result = {
            "schema": "bqg-w5-direct-all-m-shadow-equivalence-v3-adjacent-m-shards",
            "status": "INCOMPLETE_PRODUCTION_SHADOW_EQUIVALENCE",
            "passed": False,
            "comparison_performed": False,
            "error": f"{type(exc).__name__}: {exc}",
            "scope": "No equivalence claim is possible without all forty valid adjacent-M artifacts.",
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
