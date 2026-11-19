#!/usr/bin/env python3
"""Fail-closed aggregate for seed-sharded paired W5 V09/V10 reports."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

CUTOFFS = (10, 12, 14)
MODES = ("shadow_no_prune", "production")
SEEDS = (0, 2, 4, 6, 8)
TOL = 1e-8


def _decode(rows):
    return {str(k): complex(float(re), float(im)) for k, re, im in rows}


def _norm(state):
    return math.sqrt(sum(abs(v) ** 2 for v in state.values()))


def _relative(a, b):
    err = _norm({k: a.get(k, 0j) - b.get(k, 0j) for k in set(a) | set(b)})
    return err / max(_norm(b), 1e-300)


def aggregate(reports):
    expected = {(j, mode, seed) for j in CUTOFFS for mode in MODES for seed in SEEDS}
    indexed, malformed, duplicates = {}, [], []
    for report in reports:
        try:
            seed_values = tuple(int(x) for x in report["seed_Ks"])
            if len(seed_values) != 1:
                raise ValueError("each report must contain exactly one seed_K")
            key = (int(report["Jmax2"]), str(report["mode"]), seed_values[0])
        except Exception as exc:
            malformed.append(str(exc))
            continue
        if key in indexed:
            duplicates.append(key)
        indexed[key] = report

    missing = sorted(expected - set(indexed))
    unexpected = sorted(set(indexed) - expected)
    per_report_failures = []
    for key, report in indexed.items():
        failed = [name for name, value in report.get("acceptance", {}).items() if value is not True]
        if failed:
            per_report_failures.append({"cell": list(key), "failed_checks": failed})
        rows = report.get("V09_vs_V10_forward_columns", [])
        if len(rows) != 1 or int(rows[0].get("seed_K", -1)) != key[2]:
            per_report_failures.append({"cell": list(key), "failed_checks": ["exactly_one_matching_seed_column"]})

    pruning = []
    pruning_failures = []
    for cutoff in CUTOFFS:
        for seed in SEEDS:
            shadow = indexed.get((cutoff, "shadow_no_prune", seed))
            prod = indexed.get((cutoff, "production", seed))
            if shadow is None or prod is None:
                continue
            sv, pv = shadow.get("V10", {}), prod.get("V10", {})
            sc = shadow.get("V09_vs_V10_forward_columns", [{}])[0]
            pc = prod.get("V09_vs_V10_forward_columns", [{}])[0]
            row = {
                "Jmax2": cutoff, "seed_K": seed,
                "basis_identical": sv.get("basis_id_sha256") == pv.get("basis_id_sha256"),
                "closure_size_identical": sv.get("closure_size") == pv.get("closure_size"),
                "forward_matrix_digest_identical": sv.get("forward_matrix_digest") == pv.get("forward_matrix_digest"),
                "direct_adjoint_matrix_digest_identical": sv.get("direct_adjoint_matrix_digest") == pv.get("direct_adjoint_matrix_digest"),
                "V09_raw_shadow_vs_production": _relative(_decode(sc.get("V09_raw_sparse", [])), _decode(pc.get("V09_raw_sparse", []))),
                "V10_shadow_vs_production": _relative(_decode(sc.get("V10_seed_sparse", [])), _decode(pc.get("V10_seed_sparse", []))),
            }
            row["passed"] = all(row[k] for k in ("basis_identical", "closure_size_identical", "forward_matrix_digest_identical", "direct_adjoint_matrix_digest_identical")) and row["V09_raw_shadow_vs_production"] < TOL and row["V10_shadow_vs_production"] < TOL
            pruning.append(row)
            if not row["passed"]:
                pruning_failures.append(row)

    convergence, convergence_failures = [], []
    for mode in MODES:
        for seed in SEEDS:
            base = indexed.get((CUTOFFS[0], mode, seed))
            if base is None:
                continue
            for cutoff in CUTOFFS[1:]:
                candidate = indexed.get((cutoff, mode, seed))
                if candidate is None:
                    continue
                bv, cv = base.get("V10", {}), candidate.get("V10", {})
                bc = base.get("V09_vs_V10_forward_columns", [{}])[0]
                cc = candidate.get("V09_vs_V10_forward_columns", [{}])[0]
                row = {
                    "mode": mode, "seed_K": seed,
                    "reference_Jmax2": CUTOFFS[0], "candidate_Jmax2": cutoff,
                    "basis_identical": bv.get("basis_id_sha256") == cv.get("basis_id_sha256"),
                    "closure_size_identical": bv.get("closure_size") == cv.get("closure_size"),
                    "forward_matrix_digest_identical": bv.get("forward_matrix_digest") == cv.get("forward_matrix_digest"),
                    "direct_adjoint_matrix_digest_identical": bv.get("direct_adjoint_matrix_digest") == cv.get("direct_adjoint_matrix_digest"),
                    "V09_raw_relative_difference": _relative(_decode(bc.get("V09_raw_sparse", [])), _decode(cc.get("V09_raw_sparse", []))),
                    "V10_seed_relative_difference": _relative(_decode(bc.get("V10_seed_sparse", [])), _decode(cc.get("V10_seed_sparse", []))),
                }
                row["passed"] = all(row[k] for k in ("basis_identical", "closure_size_identical", "forward_matrix_digest_identical", "direct_adjoint_matrix_digest_identical")) and row["V09_raw_relative_difference"] < TOL and row["V10_seed_relative_difference"] < TOL
                convergence.append(row)
                if not row["passed"]:
                    convergence_failures.append(row)

    passed = not missing and not unexpected and not malformed and not duplicates and not per_report_failures and len(pruning) == 15 and not pruning_failures and len(convergence) == 20 and not convergence_failures
    return {
        "schema": "bqg-w5-v09-v10-seed-sharded-aggregate-v1",
        "status": "PASS_PAIRED_V09_V10_SEED_CUTOFF_AUDIT" if passed else "INCOMPLETE_OR_FAILED_PAIRED_V09_V10_AUDIT",
        "passed": passed,
        "expected_report_count": len(expected),
        "received_report_count": len(indexed),
        "missing_cells": [list(x) for x in missing],
        "unexpected_cells": [list(x) for x in unexpected],
        "malformed_reports": malformed,
        "duplicate_cells": [list(x) for x in duplicates],
        "per_report_failures": per_report_failures,
        "production_vs_no_pruning": pruning,
        "pruning_failures": pruning_failures,
        "cutoff_convergence": convergence,
        "cutoff_convergence_failures": convergence_failures,
        "acceptance": {
            "all_30_seed_cutoff_mode_reports_present": not missing and not unexpected,
            "all_per_report_checks_pass": not per_report_failures,
            "all_15_production_shadow_pairs_match": len(pruning) == 15 and not pruning_failures,
            "all_20_per_seed_cutoff_comparisons_match": len(convergence) == 20 and not convergence_failures,
        },
        "scope": "Finite node-0 W5 paired V09/V10 raw KKV structural audit, seed-sharded over K=(0,2,4,6,8), Jmax2=(10,12,14), and production/shadow modes. Not a physical H_L or continuum PASS.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reports = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(args.input_dir.glob("*.json"))]
    result = aggregate(reports)
    result["input_report_count"] = len(reports)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "passed", "expected_report_count", "received_report_count", "missing_cells", "per_report_failures", "pruning_failures", "cutoff_convergence_failures")}, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
