#!/usr/bin/env python3
"""Run fresh-process W5 cache-on/cache-off equivalence on a bounded closure."""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available,
    compare_cache_runs,
    run_one,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=int, nargs="+", default=[0])
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--max-states", type=int, default=500)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--unbounded-closure", action="store_true")
    parser.add_argument("--amp-tol", type=float, default=1e-12)
    parser.add_argument("--spectral-zero-factor", type=float, default=1000.0)
    parser.add_argument("--timeout-seconds", type=float, default=900.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    if args.amp_tol <= 0:
        parser.error("--amp-tol must be positive")
    if args.Jmax2 < 0 or args.max_states < 1 or args.max_depth < 0:
        parser.error("invalid closure bounds")

    script = Path(__file__).with_name("bqg_v10_w5_lorentzian_adjoint_closure.py")
    with tempfile.TemporaryDirectory(prefix="bqg-w5-cache-equivalence-") as temp:
        root = Path(temp)
        common = argparse.Namespace(
            Jmax2=args.Jmax2,
            max_states=args.max_states,
            max_depth=args.max_depth,
            nodes=args.nodes,
            timeout_seconds=args.timeout_seconds,
            unbounded_closure=args.unbounded_closure,
            spectral_zero_factor=args.spectral_zero_factor,
            cache_mode="on",
        )
        cache_on = run_one(
            script, root / "cache-on.json", args.amp_tol, common
        )
        cache_off_args = argparse.Namespace(**vars(common))
        cache_off_args.cache_mode = "off"
        cache_off = run_one(
            script, root / "cache-off.json", args.amp_tol, cache_off_args
        )

    if (
        closure_matrix_available(cache_on)
        and closure_matrix_available(cache_off)
    ):
        comparison = compare_cache_runs(cache_on, cache_off)
        status = comparison["status"]
        exit_code = 0 if status == "PASS_CACHE_EQUIVALENCE" else 1
        comparison_performed = True
    else:
        comparison = {
            "status": "INCOMPLETE_CACHE_EQUIVALENCE",
            "comparison_performed": False,
            "reason": (
                "Both fresh-process runs must produce complete finite closure "
                "matrices before exact matrix-and-ledger equivalence can be claimed."
            ),
        }
        status = comparison["status"]
        exit_code = 2
        comparison_performed = False

    result = {
        "status": status,
        "comparison_performed": comparison_performed,
        "cache_on": {
            "status": cache_on.get("status"),
            "closure_size": cache_on.get("closure_size"),
            "basis_id_sha256": cache_on.get("basis_id_sha256"),
            "cache_policy": cache_on.get("lower_level_pruning_ledger", {}).get("cache_policy"),
            "event_count": cache_on.get("lower_level_pruning_ledger", {}).get("event_count"),
            "process": cache_on.get("_process"),
        },
        "cache_off": {
            "status": cache_off.get("status"),
            "closure_size": cache_off.get("closure_size"),
            "basis_id_sha256": cache_off.get("basis_id_sha256"),
            "cache_policy": cache_off.get("lower_level_pruning_ledger", {}).get("cache_policy"),
            "event_count": cache_off.get("lower_level_pruning_ledger", {}).get("event_count"),
            "process": cache_off.get("_process"),
        },
        "comparison": comparison,
        "parameters": {
            "nodes": args.nodes,
            "Jmax2": args.Jmax2,
            "max_states": None if args.unbounded_closure else args.max_states,
            "max_depth": None if args.unbounded_closure else args.max_depth,
            "unbounded_closure": args.unbounded_closure,
            "amp_tol": args.amp_tol,
            "spectral_zero_factor": args.spectral_zero_factor,
            "timeout_seconds_per_process": args.timeout_seconds,
        },
        "scope": (
            "Fresh-process cache-on/cache-off comparison of forward and direct-adjoint "
            "complex matrices plus lower-level pruning, coefficient-loss, coherent-sum, "
            "spin-wall/closure and spectral-sensitivity ledgers. A capped or timed-out "
            "matrix is incomplete and cannot pass."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "comparison_performed": comparison_performed,
        "cache_on_status": cache_on.get("status"),
        "cache_off_status": cache_off.get("status"),
        "output": str(args.output),
    }, indent=2))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
