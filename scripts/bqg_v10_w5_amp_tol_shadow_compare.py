#!/usr/bin/env python3
"""Isolated production-versus-bounded-no-pruning V10 comparison.

Runs separate fresh processes. The shadow process sets amp_tol=0 and enables
the context-local no-pruning bypass at every explicitly wired lower-level
amplitude threshold. The comparison remains explicitly incomplete until the
full site inventory, non-amplitude regulator boundaries, and closure limits
are independently certified.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


COMPLETE_CLOSURE_STATUSES = {
    "PASS_ADJOINT_CONSISTENCY",
    "FAIL_ADJOINT_OR_LEAKAGE_GATE",
}


def closure_matrix_available(result: dict, *, shadow: bool = False) -> bool:
    status = result.get("status")
    if status in COMPLETE_CLOSURE_STATUSES:
        return True
    # Full shadow matrices remain labelled INCOMPLETE until the pruning-site
    # coverage proof is complete, but are valid for bounded diagnostics.
    return bool(
        shadow
        and status == "INCOMPLETE_NO_PRUNING_COVERAGE"
        and result.get("pruning_mode") == "shadow_no_prune"
        and result.get("algebraic_status") in COMPLETE_CLOSURE_STATUSES
        and isinstance(result.get("forward_B_columns_sparse"), list)
    )


def run_one(script: Path, output: Path, amp_tol: float, args,
             pruning_mode: str = "production") -> dict:
    command = [
        sys.executable, str(script),
        "--Jmax2", str(args.Jmax2),
        "--max-states", str(args.max_states),
        "--max-depth", str(args.max_depth),
        "--amp-tol", repr(amp_tol),
        "--pruning-mode", pruning_mode,
        "--output", str(output),
    ]
    command.extend(["--nodes", *[str(node) for node in args.nodes]])
    output.unlink(missing_ok=True)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, check=False,
            timeout=args.timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        # A timed-out process may have left a partial JSON file. Never treat
        # that file as a completed sparse matrix.
        output.unlink(missing_ok=True)
        stderr = exc.stderr or ""
        stdout = exc.stdout or ""
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        return {
            "status": "INCOMPLETE_TIMEOUT",
            "_process": {
                "command": command,
                "timeout_seconds": args.timeout_seconds,
                "stdout_tail": stdout[-2000:],
                "stderr_tail": stderr[-4000:],
                "partial_output_discarded": True,
            },
        }
    if not output.exists():
        return {
            "status": "INCOMPLETE_NO_RESULT_JSON",
            "_process": {
                "returncode": completed.returncode,
                "command": command,
                "stderr_tail": completed.stderr[-4000:],
                "stdout_tail": completed.stdout[-2000:],
            },
        }
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
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    if args.amp_tol <= 0:
        parser.error("--amp-tol must be positive; zero is reserved for the shadow run")
    script = Path(__file__).with_name("bqg_v10_w5_lorentzian_adjoint_closure.py")
    with tempfile.TemporaryDirectory(prefix="bqg-v10-shadow-") as temp:
        root = Path(temp)
        production = run_one(script, root / "production.json", args.amp_tol, args)
        zero_output = run_one(
            script, root / "shadow-no-prune.json", 0.0, args,
            pruning_mode="shadow_no_prune",
        )

    if not closure_matrix_available(production) or not closure_matrix_available(zero_output, shadow=True):
        timed_out = any(
            row.get("status") == "INCOMPLETE_TIMEOUT"
            for row in (production, zero_output)
        )
        result = {
            "status": (
                "INCOMPLETE_TIMEOUT" if timed_out
                else "INCOMPLETE_CLOSURE_LIMIT_REACHED"
            ),
            "production_status": production.get("status"),
            "shadow_no_prune_status": zero_output.get("status"),
            "production_process": production.get("_process"),
            "shadow_process": zero_output.get("_process"),
            "timeout_seconds_per_process": args.timeout_seconds,
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
            "shadow_pruning_mode": zero_output.get("pruning_mode"),
            "production_lower_level_pruning": production.get("lower_level_pruning_ledger"),
            "shadow_lower_level_pruning": zero_output.get("lower_level_pruning_ledger"),
            "production_status": production.get("status"),
            "shadow_no_prune_status": zero_output.get("status"),
            "shadow_algebraic_status": zero_output.get("algebraic_status"),
            "comparison": comparison,
            "scope": (
                "Isolated-process comparison of production thresholds against "
                "the bounded shadow_no_prune mode. Explicitly instrumented "
                "lower-level amplitude thresholds are bypassed in the shadow "
                "process and record would-prune events. This remains incomplete "
                "until every reachable threshold site and non-amplitude "
                "regulator boundary is independently certified."
            ),
            "lower_level_pruning_bypassed_at_instrumented_sites": True,
            "all_lower_level_pruning_bypassed": False,
            "unverified_lower_level_sites_may_remain": True,
            "unclosed_scope": [
                "Any threshold site not yet wired to PRUNE.keep",
                "The zero-aware spectral eigenvalue convention (operator-definition axis, not amplitude pruning)",
                "Spin-wall and state/depth closure caps",
                "Continuum and regulator convergence beyond this finite closure",
            ],
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
