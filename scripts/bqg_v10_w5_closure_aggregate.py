#!/usr/bin/env python3
"""Aggregate W5 closure cells and test finite-cutoff matrix convergence."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

EXPECTED_NODES = (0, 1, 2, 3, 4)
EXPECTED_SEED_KS = (0, 2, 4, 6, 8)
EXPECTED_JMAX2 = (10, 12, 14)


def _finite_positive(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value)) and float(value) > 0.0


def aggregate_records(records, nodes=EXPECTED_NODES, seed_ks=EXPECTED_SEED_KS, jmax2_values=EXPECTED_JMAX2):
    expected = {(j, n, k) for j in jmax2_values for n in nodes for k in seed_ks}
    indexed = {}
    duplicates = []
    malformed = []
    for record in records:
        try:
            ns = record["source_nodes"]
            ks = record["seed_Ks"]
            if len(ns) != 1 or len(ks) != 1:
                raise ValueError("each record must contain exactly one node and one seed K")
            key = (int(record["Jmax2"]), int(ns[0]), int(ks[0]))
        except Exception as exc:
            malformed.append(str(exc))
            continue
        if key in indexed:
            duplicates.append(key)
        indexed[key] = record

    present = set(indexed)
    missing = sorted(expected - present)
    unexpected = sorted(present - expected)
    cell_failures = []
    for key in sorted(expected & present):
        r = indexed[key]
        d = r.get("direct_adjoint_relative_residual")
        leakage = r.get("max_component_leakage")
        max_spin = r.get("max_spin_reached")
        max_intermediate = r.get("max_intermediate_spin_reached")
        wall_count = r.get("spin_wall_omitted_branch_count")
        if r.get("status") != "PASS_ADJOINT_CONSISTENCY" or r.get("algebraic_status") != "PASS_ADJOINT_CONSISTENCY":
            cell_failures.append({"cell": key, "reason": "closure/algebraic status is not PASS"})
        if not r.get("matrices_present") or not r.get("basis_id_sha256"):
            cell_failures.append({"cell": key, "reason": "complete matrices or basis fingerprint missing"})
        if not _finite_positive(r.get("forward_frobenius_norm")) or not _finite_positive(r.get("direct_adjoint_frobenius_norm")):
            cell_failures.append({"cell": key, "reason": "forward/adjoint operator norm is zero or invalid"})
        if not isinstance(d, (int, float)) or not math.isfinite(float(d)) or float(d) >= 1e-8:
            cell_failures.append({"cell": key, "reason": "direct-adjoint residual is not below 1e-8"})
        if not isinstance(leakage, (int, float)) or not math.isfinite(float(leakage)) or float(leakage) >= 1e-9:
            cell_failures.append({"cell": key, "reason": "component leakage is not below 1e-9"})
        if r.get("scalar_projection_rejected_terms") != 0:
            cell_failures.append({"cell": key, "reason": "scalar projection rejected terms are nonzero"})
        if not isinstance(max_spin, (int, float)) or float(max_spin) > 3.5 + 1e-12:
            cell_failures.append({"cell": key, "reason": "maximum spin exceeds 3.5"})
        if not isinstance(max_intermediate, (int, float)) or float(max_intermediate) > 3.5 + 1e-12:
            cell_failures.append({"cell": key, "reason": "intermediate spin exceeds 3.5"})
        if wall_count != 0:
            cell_failures.append({"cell": key, "reason": "spin-wall omitted branches are nonzero"})
        if r.get("sine_engine_gate_passed") is not True:
            cell_failures.append({"cell": key, "reason": "underlying sine engine gate did not pass"})

    convergence_failures = []
    for node in nodes:
        for seed_k in seed_ks:
            cells = [indexed.get((j, node, seed_k)) for j in jmax2_values]
            if any(c is None for c in cells):
                continue
            reference = cells[0]
            for j, candidate in zip(jmax2_values[1:], cells[1:]):
                if candidate.get("basis_id_sha256") != reference.get("basis_id_sha256"):
                    convergence_failures.append({"node": node, "seed_K": seed_k, "Jmax2": j, "reason": "reachable basis fingerprint changed with cutoff"})
                if candidate.get("closure_size") != reference.get("closure_size"):
                    convergence_failures.append({"node": node, "seed_K": seed_k, "Jmax2": j, "reason": "closure size changed with cutoff"})
                for field in ("forward_matrix_digest_10sig", "direct_adjoint_matrix_digest_10sig"):
                    if candidate.get(field) != reference.get(field):
                        convergence_failures.append({"node": node, "seed_K": seed_k, "Jmax2": j, "reason": f"{field} changed with cutoff"})
                for field in ("forward_frobenius_norm", "direct_adjoint_frobenius_norm"):
                    a, b = reference.get(field), candidate.get(field)
                    if not _finite_positive(a) or not _finite_positive(b):
                        convergence_failures.append({"node": node, "seed_K": seed_k, "Jmax2": j, "reason": f"{field} missing or invalid"})
                    elif abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), 1e-300) > 1e-8:
                        convergence_failures.append({"node": node, "seed_K": seed_k, "Jmax2": j, "reason": f"{field} not cutoff-stable"})

    passed = not (missing or unexpected or duplicates or malformed or cell_failures or convergence_failures)
    return {
        "schema": "bqg-w5-closure-matrix-aggregate-v1",
        "status": "PASS_W5_CLOSURE_MATRIX_AND_CUTOFF_CONVERGENCE" if passed else "INCOMPLETE_OR_FAILED_W5_CLOSURE_MATRIX",
        "passed": passed,
        "expected_cell_count": len(expected),
        "received_unique_cell_count": len(present),
        "missing_cells": [list(x) for x in missing],
        "unexpected_cells": [list(x) for x in unexpected],
        "duplicate_cells": [list(x) for x in duplicates],
        "malformed_records": malformed,
        "cell_failures": cell_failures,
        "cutoff_convergence_failures": convergence_failures,
        "cutoffs_Jmax2": list(jmax2_values),
        "nodes": list(nodes),
        "seed_Ks": list(seed_ks),
        "acceptance": {
            "all_75_production_cells_pass": not cell_failures and not missing and not unexpected,
            "basis_and_matrix_fingerprints_stable_across_cutoffs": not convergence_failures,
            "matrix_fingerprint_precision": "10 significant digits",
            "direct_adjoint_relative_residual_max": 1e-8,
            "component_leakage_max": 1e-9,
            "max_spin": 3.5,
            "spin_wall_omitted_branches": 0,
        },
        "scope": "Finite W5-seeded reachable closure only; does not establish the physical Lorentzian ordering, continuum limit, or full BQG constraint.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = sorted(args.input_dir.glob("*.summary.json"))
    records = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    report = aggregate_records(records)
    report["summary_files_read"] = len(paths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
