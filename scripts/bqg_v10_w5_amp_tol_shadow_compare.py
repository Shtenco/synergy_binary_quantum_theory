#!/usr/bin/env python3
"""Isolated V10 output-threshold sensitivity comparison.

This is deliberately NOT called a no-pruning reference evaluator: lower-level
operator pruning remains active in the pinned production modules.  It runs the
same V10 source twice in separate Python processes (fresh import namespaces and
caches), once at the requested amp_tol and once at amp_tol=0, then compares the
complete encoded complex sparse columns.  A full no-pruning verdict is refused
until every lower-level pruning site has an instrumented bypass/ledger.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def run_one(script: Path, output: Path, amp_tol: float, args) -> dict:
    command = [
        sys.executable, str(script),
        "--Jmax2", str(args.Jmax2),
        "--max-states", str(args.max_states),
        "--max-depth", str(args.max_depth),
        "--amp-tol", repr(amp_tol),
        "--output", str(output),
    ]
    command.extend(["--nodes", *[str(node) for node in args.nodes]])
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if not output.exists():
        raise RuntimeError(
            f"V10 subprocess produced no JSON (exit={completed.returncode}): "
            f"{completed.stderr[-4000:]}"
        )
    result = json.loads(output.read_text(encoding="utf-8"))
    result["_process"] = {
        "returncode": completed.returncode,
        "command": command,
        "stderr_tail": completed.stderr[-2000:],
    }
    return result


def unpack_columns(result: dict) -> dict:
    basis = result["basis_keys_repr"]
    columns = result["forward_B_columns_sparse"]
    decoded = {}
    for input_index, encoded in enumerate(columns):
        input_key = basis[input_index]
        vector = {}
        for output_index, real, imag in encoded:
            output_key = basis[output_index]
            value = complex(real, imag)
            if value != 0j:
                vector[output_key] = value
        decoded[input_key] = vector
    return decoded


def compare_vectors(reference: dict, candidate: dict) -> dict:
    input_keys = set(reference) | set(candidate)
    residual2 = 0.0
    reference2 = 0.0
    max_abs = 0.0
    support_symmetric_difference = 0
    column_rows = []
    for input_key in sorted(input_keys):
        a = reference.get(input_key, {})
        b = candidate.get(input_key, {})
        keys = set(a) | set(b)
        diff2 = sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys)
        ref2 = sum(abs(v) ** 2 for v in a.values())
        residual2 += diff2
        reference2 += ref2
        max_abs = max(
            max_abs,
            max((abs(a.get(k, 0j) - b.get(k, 0j)) for k in keys), default=0.0),
        )
        support_symmetric_difference += len(set(a) ^ set(b))
        column_rows.append({
            "input_key": input_key,
            "reference_support": len(a),
            "candidate_support": len(b),
            "support_symmetric_difference": len(set(a) ^ set(b)),
            "absolute_l2_residual": math.sqrt(diff2),
        })
    residual = math.sqrt(residual2)
    scale = math.sqrt(reference2)
    return {
        "input_column_count_union": len(input_keys),
        "support_symmetric_difference_total": support_symmetric_difference,
        "absolute_l2_residual": residual,
        "relative_l2_residual": residual / max(scale, 1e-300),
        "max_coefficient_absolute_residual": max_abs,
        "per_input_column": column_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--max-states", type=int, default=500)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--amp-tol", type=float, default=1e-12)
    parser.add_argument("--nodes", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.amp_tol <= 0:
        parser.error("--amp-tol must be positive; zero is reserved for the shadow run")
    script = Path(__file__).with_name("bqg_v10_w5_lorentzian_adjoint_closure.py")
    with tempfile.TemporaryDirectory(prefix="bqg-v10-shadow-") as temp:
        root = Path(temp)
        production = run_one(script, root / "production.json", args.amp_tol, args)
        zero_output = run_one(script, root / "zero-output-threshold.json", 0.0, args)

    complete_statuses = {
        "PASS_ADJOINT_CONSISTENCY",
        "FAIL_ADJOINT_OR_LEAKAGE_GATE",
    }
    if production.get("status") not in complete_statuses or zero_output.get("status") not in complete_statuses:
        result = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "production_status": production.get("status"),
            "zero_output_threshold_status": zero_output.get("status"),
            "comparison_performed": False,
        }
        exit_code = 2
    else:
        comparison = compare_vectors(
            unpack_columns(zero_output), unpack_columns(production)
        )
        result = {
            "status": "INCOMPLETE_NO_PRUNING_REFERENCE",
            "comparison_performed": True,
            "production_amp_tol": args.amp_tol,
            "shadow_amp_tol": 0.0,
            "production_status": production.get("status"),
            "zero_output_threshold_status": zero_output.get("status"),
            "comparison": comparison,
            "scope": (
                "Isolated-process comparison of V10's final coherent output "
                "threshold only. Lower-level amplitude pruning in imported "
                "operator modules remains active; this is not a no-pruning "
                "reference and cannot support ZERO_LOSS."
            ),
            "lower_level_pruning_active": True,
            "fresh_processes_and_caches": True,
            "same_source_file": str(script),
            "parameters": {
                "nodes": args.nodes,
                "Jmax2": args.Jmax2,
                "max_states": args.max_states,
                "max_depth": args.max_depth,
            },
        }
        exit_code = 2

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "comparison_performed": result.get("comparison_performed", False),
        "output": str(args.output),
    }, indent=2))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
