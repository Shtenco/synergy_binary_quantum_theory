#!/usr/bin/env python3
"""Run the full finite-W5 Jmax2 x amp_tol sparse-matrix convergence grid.

Each grid point is executed in a fresh process. A shadow_no_prune matrix is
computed once per Jmax2 and compared with every production amp_tol at that
wall. Production matrices are also compared across amp_tol and adjacent
Jmax2 values using full complex sparse support, not norms alone.

The shadow shares the production operator implementation and remains a
partial reference until the independent derivation and pruning-site coverage
proof are complete. This runner never upgrades that shadow into ZERO_LOSS.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
from pathlib import Path

from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available,
    compare_vectors,
    run_one,
    unpack_columns,
)


def matrix_digest(matrix: dict) -> str:
    rows = []
    for input_key in sorted(matrix):
        for output_key, value in sorted(matrix[input_key].items()):
            z = complex(value)
            rows.append([input_key, output_key, float(z.real), float(z.imag)])
    payload = json.dumps(rows, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def run_grid(args) -> dict:
    walls = sorted(set(int(value) for value in args.Jmax2_grid))
    tolerances = sorted(set(float(value) for value in args.amp_tol_grid))
    if not walls or any(value < 10 for value in walls):
        raise ValueError("Jmax2 grid must be nonempty and every wall must be >=10")
    if not tolerances or any(not math.isfinite(value) or value <= 0 for value in tolerances):
        raise ValueError("amp_tol grid must be nonempty and finite-positive")
    if not math.isfinite(args.spectral_zero_factor) or args.spectral_zero_factor <= 0:
        raise ValueError("spectral_zero_factor must be finite and positive")

    grid_rows = []
    production_matrices = {}
    production_adjoint_matrices = {}
    shadow_matrices = {}
    shadow_adjoint_matrices = {}
    with tempfile.TemporaryDirectory(prefix="bqg-v10-convergence-grid-") as temp:
        root = Path(temp)
        script = Path(__file__).with_name("bqg_v10_w5_lorentzian_adjoint_closure.py")
        for wall in walls:
            wall_args = argparse.Namespace(**vars(args))
            wall_args.Jmax2 = wall
            shadow = run_one(
                script, root / f"shadow-J{wall}.json", 0.0, wall_args,
                pruning_mode="shadow_no_prune",
            )
            shadow_ok = closure_matrix_available(shadow, shadow=True)
            if shadow_ok:
                shadow_matrices[wall] = unpack_columns(shadow)
                shadow_adjoint_matrices[wall] = unpack_columns(
                    shadow, "direct_Bdag_columns_sparse"
                )

            for tol in tolerances:
                production = run_one(
                    script, root / f"production-J{wall}-tol{tol:.17g}.json",
                    tol, wall_args, pruning_mode="production",
                )
                production_ok = closure_matrix_available(production)
                row = {
                    "Jmax2": wall,
                    "amp_tol": tol,
                    "production_status": production.get("status"),
                    "shadow_status": shadow.get("status"),
                    "production_matrix_complete": production_ok,
                    "shadow_matrix_complete": shadow_ok,
                    "production_process": production.get("_process"),
                    "shadow_process": shadow.get("_process"),
                    "production_loss_summary": production.get("coefficient_loss_ledger", {}).get("summary"),
                    "shadow_loss_summary": shadow.get("coefficient_loss_ledger", {}).get("summary"),
                    "production_diagnostics": production.get("diagnostics", production.get("partial_diagnostics")),
                    "shadow_diagnostics": shadow.get("diagnostics", shadow.get("partial_diagnostics")),
                    "production_spectral_sensitivity": production.get("spectral_sensitivity"),
                    "shadow_spectral_sensitivity": shadow.get("spectral_sensitivity"),
                }
                if production_ok:
                    matrix = unpack_columns(production)
                    adjoint_matrix = unpack_columns(production, "direct_Bdag_columns_sparse")
                    production_matrices[(wall, tol)] = matrix
                    production_adjoint_matrices[(wall, tol)] = adjoint_matrix
                    row["production_matrix_sha256"] = matrix_digest(matrix)
                    row["production_adjoint_matrix_sha256"] = matrix_digest(adjoint_matrix)
                    row["production_closure_size"] = production.get("closure_size")
                    row["production_basis_id_sha256"] = production.get("basis_id_sha256")
                if shadow_ok:
                    row["shadow_matrix_sha256"] = matrix_digest(shadow_matrices[wall])
                    row["shadow_adjoint_matrix_sha256"] = matrix_digest(shadow_adjoint_matrices[wall])
                    row["shadow_closure_size"] = shadow.get("closure_size")
                    row["shadow_basis_id_sha256"] = shadow.get("basis_id_sha256")
                if production_ok and shadow_ok:
                    row["production_vs_shadow_forward"] = compare_vectors(
                        shadow_matrices[wall], production_matrices[(wall, tol)]
                    )
                    row["production_vs_shadow_direct_adjoint"] = compare_vectors(
                        shadow_adjoint_matrices[wall], production_adjoint_matrices[(wall, tol)]
                    )
                else:
                    row["production_vs_shadow_forward"] = None
                    row["production_vs_shadow_direct_adjoint"] = None
                grid_rows.append(row)

    amp_convergence = []
    for wall in walls:
        available = [(tol, production_matrices[(wall, tol)])
                     for tol in tolerances if (wall, tol) in production_matrices]
        if len(available) >= 2:
            baseline_tol, baseline = available[0]
            for tol, matrix in available[1:]:
                amp_convergence.append({
                    "Jmax2": wall,
                    "baseline_amp_tol": baseline_tol,
                    "candidate_amp_tol": tol,
                    "forward_comparison": compare_vectors(baseline, matrix),
                    "direct_adjoint_comparison": compare_vectors(
                        production_adjoint_matrices[(wall, baseline_tol)],
                        production_adjoint_matrices[(wall, tol)],
                    ),
                })

    wall_convergence = []
    for tol in tolerances:
        for wall_a, wall_b in zip(walls, walls[1:]):
            if ((wall_a, tol) not in production_matrices
                    or (wall_b, tol) not in production_matrices):
                continue
            wall_convergence.append({
                "amp_tol": tol,
                "lower_Jmax2": wall_a,
                "higher_Jmax2": wall_b,
                "forward_comparison": compare_vectors(
                    production_matrices[(wall_a, tol)],
                    production_matrices[(wall_b, tol)],
                ),
                "direct_adjoint_comparison": compare_vectors(
                    production_adjoint_matrices[(wall_a, tol)],
                    production_adjoint_matrices[(wall_b, tol)],
                ),
            })

    expected_cells = len(walls) * len(tolerances)
    completed_cells = sum(
        row["production_matrix_complete"] and row["shadow_matrix_complete"]
        for row in grid_rows
    )
    all_production_complete = len(production_matrices) == expected_cells
    all_shadow_complete = len(shadow_matrices) == len(walls)
    return {
        "status": (
            "COMPLETED_FINITE_GRID_SHADOW_REFERENCE_INCOMPLETE"
            if completed_cells == expected_cells
            else "INCOMPLETE_CONVERGENCE_GRID"
        ),
        "Jmax2_grid": walls,
        "amp_tol_grid": tolerances,
        "source_nodes": list(args.nodes),
        "closure": {
            "unbounded": bool(args.unbounded_closure),
            "max_states": None if args.unbounded_closure else args.max_states,
            "max_depth": None if args.unbounded_closure else args.max_depth,
        },
        "cache_mode": args.cache_mode,
        "spectral_zero_factor": args.spectral_zero_factor,
        "timeout_seconds_per_process": args.timeout_seconds,
        "expected_cells": expected_cells,
        "completed_production_and_shadow_cells": int(completed_cells),
        "all_production_matrices_complete": all_production_complete,
        "all_shadow_matrices_complete": all_shadow_complete,
        "grid": grid_rows,
        "amp_tol_convergence": amp_convergence,
        "Jmax2_convergence": wall_convergence,
        "reference_scope": (
            "shadow_no_prune uses the same operator implementation and bypasses "
            "only explicitly instrumented amplitude thresholds; it is not an "
            "independent derivation and is not a ZERO_LOSS certificate."
        ),
        "zero_loss_proved": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Jmax2-grid", type=int, nargs="+", required=True)
    parser.add_argument("--amp-tol-grid", type=float, nargs="+", required=True)
    parser.add_argument("--max-states", type=int, default=500)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--unbounded-closure", action="store_true")
    parser.add_argument("--nodes", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--cache-mode", choices=("on", "off"), default="on")
    parser.add_argument("--spectral-zero-factor", type=float, default=1000.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    try:
        result = run_grid(args)
    except ValueError as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "expected_cells": result["expected_cells"],
        "completed_cells": result["completed_production_and_shadow_cells"],
        "output": str(args.output),
    }, indent=2))
    return 0 if result["status"].startswith("COMPLETED_") else 2


if __name__ == "__main__":
    raise SystemExit(main())
