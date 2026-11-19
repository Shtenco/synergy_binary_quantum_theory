#!/usr/bin/env python3
"""Compare the finite G3 Galerkin heat trace with the audited G0-G2 levels.

This is a quantitative finite-level comparison only. It deliberately does not
claim a continuum limit or a stable spectral-dimension plateau.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def compare_reports(refinement_path: Path, g3_path: Path) -> dict:
    refinement = json.loads(refinement_path.read_text(encoding="utf-8"))
    g3 = json.loads(g3_path.read_text(encoding="utf-8"))
    if not refinement.get("passed") or not g3.get("passed"):
        raise ValueError("all finite-level source reports must have passed")
    levels = refinement.get("levels")
    if not isinstance(levels, list):
        raise ValueError("G0-G2 refinement levels are missing")
    by_generation = {int(level["generation"]): level for level in levels}
    if set(by_generation) != {0, 1, 2}:
        raise ValueError("comparison requires exactly the audited G0-G2 levels")
    if int(g3.get("generation", -1)) != 3:
        raise ValueError("G3 report has the wrong generation")
    g3_spectral = g3.get("sparse_galerkin")
    if not isinstance(g3_spectral, dict):
        raise ValueError("G3 Galerkin spectral report is missing")

    g3_heat = {float(row["t"]): float(row["value"])
               for row in g3_spectral["normalized_heat_trace"]}
    g3_dimension = {float(row["t"]): float(row["spectral_dimension"])
                    for row in g3_spectral["heat_trace_normalized"]}
    level_summaries = []
    for generation in range(3):
        level = by_generation[generation]
        gap = level["operator_spectral_gap"]
        level_summaries.append({
            "generation": generation,
            "dimension": int(level["operator_dimension"]),
            "sparse_nnz": int(level["operator_sparse_nnz"]),
            "spectral_gap": gap,
            "spectral_gap_ratio_to_g3": (
                float(gap) / float(g3_spectral["spectral_gap"])
                if gap is not None and g3_spectral["spectral_gap"] is not None
                else None
            ),
            "minimum_eigenvalue": float(level["operator_minimum_eigenvalue"]),
            "maximum_eigenvalue": float(level["operator_largest_eigenvalue"]),
            "trace": float(level["operator_trace"]),
            "zero_mode_count": int(level["operator_zero_mode_count"]),
            "measure_total": float(level["pushforward_measure_total"]),
        })
    level_summaries.append({
        "generation": 3,
        "dimension": int(g3_spectral["dimension"]),
        "sparse_nnz": int(g3_spectral["nonzero_count"]),
        "spectral_gap": g3_spectral["spectral_gap"],
        "spectral_gap_ratio_to_g3": 1.0 if g3_spectral["spectral_gap"] is not None else None,
        "minimum_eigenvalue": float(g3_spectral["minimum_eigenvalue"]),
        "maximum_eigenvalue": float(g3_spectral["maximum_eigenvalue"]),
        "trace": float(g3_spectral["trace"]),
        "zero_mode_count": int(g3_spectral["zero_mode_count"]),
        "measure_total": float(g3_spectral["measure_total"]),
    })

    comparisons = []
    for generation in range(3):
        level = by_generation[generation]
        heat = {float(row["time"]): float(row["normalized_heat_trace"])
                for row in level["heat_trace"]}
        dimensions = {
            float(row["time"]): float(row["effective_spectral_dimension"])
            for row in level["effective_spectral_dimension"]
        }
        shared_heat_times = sorted(set(heat) & set(g3_heat))
        shared_dimension_times = sorted(set(dimensions) & set(g3_dimension))
        if len(shared_heat_times) < 5:
            raise ValueError(f"insufficient shared heat-trace times for G{generation}")
        heat_deltas = [
            {"t": t, "absolute_delta": abs(g3_heat[t] - heat[t]),
             "g3": g3_heat[t], f"g{generation}": heat[t]}
            for t in shared_heat_times
        ]
        dimension_deltas = [
            {"t": t, "absolute_delta": abs(g3_dimension[t] - dimensions[t]),
             "g3": g3_dimension[t], f"g{generation}": dimensions[t]}
            for t in shared_dimension_times
        ]
        comparisons.append({
            "generation": generation,
            "shared_heat_trace_times": shared_heat_times,
            "max_absolute_heat_trace_delta": max(x["absolute_delta"] for x in heat_deltas),
            "heat_trace_deltas": heat_deltas,
            "shared_effective_dimension_times": shared_dimension_times,
            "max_absolute_effective_dimension_delta": (
                max((x["absolute_delta"] for x in dimension_deltas), default=None)
            ),
            "effective_dimension_deltas": dimension_deltas,
        })

    result = {
        "schema": "q2-g0-g3-galerkin-refinement-comparison-v1",
        "status": "PASS_FINITE_G0_G3_REFINEMENT_COMPARISON",
        "passed": True,
        "levels_compared": [0, 1, 2, 3],\n        "level_summaries": level_summaries,
        "operator_semantics": (
            "Mass-normalized finite Galerkin stiffness from PL dual-edge "
            "multiplicities and the pushforward fibre measure; geometric "
            "operator, not physical QFI conductance."
        ),
        "g3_dimension": int(g3_spectral["dimension"]),
        "g3_spectral_gap": g3_spectral["spectral_gap"],
        "g3_zero_mode_count": int(g3_spectral["zero_mode_count"]),
        "g3_measure_total": float(g3_spectral["measure_total"]),
        "comparisons_to_g3": comparisons,
        "stability_conclusion": (
            "Finite-level heat-trace and effective-dimension differences are "
            "quantified. A refinement-stable plateau and controlled continuum "
            "limit are not established by G0-G3."
        ),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refinement", type=Path, required=True)
    parser.add_argument("--g3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare_reports(args.refinement, args.g3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
