#!/usr/bin/env python3
"""Aggregate paired V09/V10 W5 runs across pruning modes and spin cutoffs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

CUTOFFS = (10, 12, 14)
MODES = ("shadow_no_prune", "production")
SEEDS = (0, 2, 4, 6, 8)


def _decode(rows):
    out = {}
    for key, real, imag in rows:
        out[str(key)] = complex(float(real), float(imag))
    return out


def _norm(state):
    result = 0.0
    for value in state.values():
        result = math.hypot(result, abs(value))
    return result


def _relative(a, b):
    error = 0.0
    for key in set(a) | set(b):
        error = math.hypot(error, abs(a.get(key, 0j) - b.get(key, 0j)))
    return error / max(_norm(b), 1e-300)


def aggregate(reports):
    expected = {(cutoff, mode) for cutoff in CUTOFFS for mode in MODES}
    indexed = {}
    malformed = []
    for report in reports:
        try:
            key = (int(report["Jmax2"]), str(report["mode"]))
        except Exception as exc:
            malformed.append(str(exc))
            continue
        if key in indexed:
            malformed.append(f"duplicate report for {key}")
        indexed[key] = report

    missing = sorted(expected - set(indexed))
    unexpected = sorted(set(indexed) - expected)
    report_acceptance_failures = []
    for key, report in indexed.items():
        failed = [name for name, value in report.get("acceptance", {}).items() if value is not True]
        if failed:
            report_acceptance_failures.append({"cell": key, "failed_checks": failed})

    pruning_comparisons = []
    for cutoff in CUTOFFS:
        shadow = indexed.get((cutoff, "shadow_no_prune"))
        production = indexed.get((cutoff, "production"))
        if shadow is None or production is None:
            continue
        sv10, pv10 = shadow.get("V10", {}), production.get("V10", {})
        pruning_comparisons.append({
            "Jmax2": cutoff,
            "V10_basis_identical": sv10.get("basis_id_sha256") == pv10.get("basis_id_sha256"),
            "V10_closure_size_identical": sv10.get("closure_size") == pv10.get("closure_size"),
            "V10_forward_matrix_digest_identical": sv10.get("forward_matrix_digest") == pv10.get("forward_matrix_digest"),
            "V10_direct_adjoint_matrix_digest_identical": sv10.get("direct_adjoint_matrix_digest") == pv10.get("direct_adjoint_matrix_digest"),
            "seed_column_differences": [],
        })
        shadow_cols = {int(x["seed_K"]): x for x in shadow.get("V09_vs_V10_forward_columns", [])}
        prod_cols = {int(x["seed_K"]): x for x in production.get("V09_vs_V10_forward_columns", [])}
        for seed in SEEDS:
            if seed not in shadow_cols or seed not in prod_cols:
                pruning_comparisons[-1]["seed_column_differences"].append({
                    "seed_K": seed, "status": "MISSING_COLUMN"
                })
                continue
            s, p = shadow_cols[seed], prod_cols[seed]
            pruning_comparisons[-1]["seed_column_differences"].append({
                "seed_K": seed,
                "V09_raw_shadow_vs_production": _relative(_decode(s["V09_raw_sparse"]), _decode(p["V09_raw_sparse"])),
                "V10_shadow_vs_production": _relative(_decode(s["V10_seed_sparse"]), _decode(p["V10_seed_sparse"])),
            })

    cutoff_convergence = []
    for mode in MODES:
        base = indexed.get((CUTOFFS[0], mode))
        if base is None:
            continue
        base_v10 = base.get("V10", {})
        base_cols = {int(x["seed_K"]): x for x in base.get("V09_vs_V10_forward_columns", [])}
        for cutoff in CUTOFFS[1:]:
            candidate = indexed.get((cutoff, mode))
            if candidate is None:
                continue
            cv10 = candidate.get("V10", {})
            row = {
                "mode": mode,
                "reference_Jmax2": CUTOFFS[0],
                "candidate_Jmax2": cutoff,
                "basis_identical": base_v10.get("basis_id_sha256") == cv10.get("basis_id_sha256"),
                "closure_size_identical": base_v10.get("closure_size") == cv10.get("closure_size"),
                "forward_matrix_digest_identical": base_v10.get("forward_matrix_digest") == cv10.get("forward_matrix_digest"),
                "direct_adjoint_matrix_digest_identical": base_v10.get("direct_adjoint_matrix_digest") == cv10.get("direct_adjoint_matrix_digest"),
                "seed_column_differences": [],
            }
            cand_cols = {int(x["seed_K"]): x for x in candidate.get("V09_vs_V10_forward_columns", [])}
            for seed in SEEDS:
                if seed not in base_cols or seed not in cand_cols:
                    row["seed_column_differences"].append({"seed_K": seed, "status": "MISSING_COLUMN"})
                    continue
                a, b = base_cols[seed], cand_cols[seed]
                row["seed_column_differences"].append({
                    "seed_K": seed,
                    "V09_raw_cutoff_relative_difference": _relative(_decode(a["V09_raw_sparse"]), _decode(b["V09_raw_sparse"])),
                    "V10_seed_cutoff_relative_difference": _relative(_decode(a["V10_seed_sparse"]), _decode(b["V10_seed_sparse"])),
                })
            cutoff_convergence.append(row)

    serializer_losses = []
    for key, report in sorted(indexed.items()):
        for column in report.get("V09_columns", []):
            serializer_losses.append({
                "Jmax2": key[0], "mode": key[1], "seed_K": column.get("seed_K"),
                "A": column.get("A_serializer_loss"),
                "B_gauss": column.get("B_gauss_serializer_loss"),
                "B_covariant": column.get("B_covariant_serializer_loss"),
            })

    pruning_failures = []
    for row in pruning_comparisons:
        if (not row["V10_basis_identical"] or not row["V10_closure_size_identical"]
                or not row["V10_forward_matrix_digest_identical"]
                or not row["V10_direct_adjoint_matrix_digest_identical"]):
            pruning_failures.append({"Jmax2": row["Jmax2"], "reason": "production/shadow basis differs"})
        for col in row["seed_column_differences"]:
            if col.get("status") == "MISSING_COLUMN" or col.get("V09_raw_shadow_vs_production", 1.0) > 1e-10 or col.get("V10_shadow_vs_production", 1.0) > 1e-10:
                pruning_failures.append({"Jmax2": row["Jmax2"], "seed_K": col.get("seed_K"), "reason": "production/shadow seed column differs or is missing"})

    convergence_failures = []
    for row in cutoff_convergence:
        if (not row["basis_identical"] or not row["closure_size_identical"]
                or not row["forward_matrix_digest_identical"]
                or not row["direct_adjoint_matrix_digest_identical"]):
            convergence_failures.append({"mode": row["mode"], "Jmax2": row["candidate_Jmax2"], "reason": "reachable basis changed with cutoff"})
        for col in row["seed_column_differences"]:
            if col.get("status") == "MISSING_COLUMN" or col.get("V09_raw_cutoff_relative_difference", 1.0) > 1e-8 or col.get("V10_seed_cutoff_relative_difference", 1.0) > 1e-8:
                convergence_failures.append({"mode": row["mode"], "Jmax2": row["candidate_Jmax2"], "seed_K": col.get("seed_K"), "reason": "seed column not cutoff-stable"})

    passed = not (
        missing or unexpected or malformed or report_acceptance_failures
        or pruning_failures or convergence_failures
    )
    return {
        "schema": "bqg-w5-v09-v10-paired-aggregate-v1",
        "status": "PASS_PAIRED_V09_V10_PRUNING_AND_CUTOFF_AUDIT" if passed else "INCOMPLETE_OR_FAILED_PAIRED_V09_V10_AUDIT",
        "passed": passed,
        "expected_report_count": len(expected),
        "received_report_count": len(indexed),
        "missing_cells": [list(x) for x in missing],
        "unexpected_cells": [list(x) for x in unexpected],
        "malformed_reports": malformed,
        "report_acceptance_failures": report_acceptance_failures,
        "production_vs_no_pruning": pruning_comparisons,
        "pruning_failures": pruning_failures,
        "cutoff_convergence": cutoff_convergence,
        "cutoff_convergence_failures": convergence_failures,
        "V09_serializer_loss_by_seed_and_run": serializer_losses,
        "acceptance": {
            "all_six_node0_mode_cutoff_reports_present": not missing and not unexpected,
            "all_per_report_checks_pass": not report_acceptance_failures,
            "production_matches_no_pruning_seed_columns": not pruning_failures,
            "five_seed_columns_stable_across_Jmax2_10_12_14": not convergence_failures,
        },
        "scope": "Node 0, five W5 seed states, finite raw KKV structural audit only; not a physical H_L or continuum PASS.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reports = [json.loads(path.read_text(encoding="utf-8")) for path in sorted(args.input_dir.glob("*.json"))]
    result = aggregate(reports)
    result["input_report_count"] = len(reports)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
