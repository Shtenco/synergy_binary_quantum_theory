#!/usr/bin/env python3
"""Create a compact, cutoff-comparable evidence record for one W5 closure cell."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def matrix_digest(columns):
    h = hashlib.sha256()
    h.update(f"columns:{len(columns)}\n".encode())
    for column_index, column in enumerate(columns):
        for row_index, real, imag in sorted(column, key=lambda item: int(item[0])):
            h.update(
                f"{column_index}:{int(row_index)}:{float(real):.10g}:{float(imag):.10g}\n".encode()
            )
    return h.hexdigest()


def summarize(result):
    diagnostics = result.get("diagnostics") or result.get("partial_diagnostics") or {}
    forward = result.get("forward_B_columns_sparse")
    adjoint = result.get("direct_Bdag_columns_sparse")
    complete_matrices = isinstance(forward, list) and isinstance(adjoint, list)
    return {
        "schema": "bqg-w5-closure-cell-summary-v1",
        "status": result.get("status"),
        "algebraic_status": result.get("algebraic_status"),
        "source_nodes": result.get("source_nodes", []),
        "seed_Ks": result.get("seed_Ks", []),
        "Jmax2": result.get("Jmax2"),
        "pruning_mode": result.get("pruning_mode"),
        "cache_mode": result.get("cache_mode"),
        "closure_size": result.get("closure_size"),
        "closure_size_discovered": result.get("closure_size_discovered"),
        "closure_depth_max": result.get("closure_depth_max"),
        "basis_id_sha256": result.get("basis_id_sha256"),
        "forward_frobenius_norm": diagnostics.get("forward_frobenius_norm"),
        "direct_adjoint_frobenius_norm": diagnostics.get("direct_adjoint_frobenius_norm"),
        "direct_adjoint_relative_residual": diagnostics.get("direct_adjoint_relative_residual"),
        "max_component_leakage": diagnostics.get("max_component_leakage"),
        "scalar_projection_rejected_terms": diagnostics.get("scalar_projection_rejected_terms"),
        "max_spin_reached": diagnostics.get("max_spin_reached"),
        "max_intermediate_spin_reached": diagnostics.get("max_intermediate_spin_reached"),
        "spin_wall_omitted_branch_count": diagnostics.get("spin_wall_omitted_branch_count"),
        "sine_engine_gate_passed": diagnostics.get("sine_engine_gate_passed"),
        "lower_level_pruning_coverage_complete": result.get("lower_level_pruning_coverage_complete"),
        "matrices_present": complete_matrices,
        "forward_matrix_digest_10sig": matrix_digest(forward) if complete_matrices else None,
        "direct_adjoint_matrix_digest_10sig": matrix_digest(adjoint) if complete_matrices else None,
        "claim_scope": result.get("claim_scope"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = json.loads(args.input.read_text(encoding="utf-8"))
    record = summarize(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
