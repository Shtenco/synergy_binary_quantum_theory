#!/usr/bin/env python3
"""Fail-closed cache-on/cache-off comparison for complete finite W5 matrices."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

NODES = tuple(range(5))
SEED_KS = (0, 2, 4, 6, 8)
DEFAULT_JMAX2 = (10, 12, 14)
MATRIX_FIELDS = ("forward_B_columns_sparse", "direct_Bdag_columns_sparse")
TOLERANCE = 1e-12


def _load(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"{path.name}: invalid JSON: {type(exc).__name__}: {exc}") from exc


def _sparse_residual(left, right, basis_dimension):
    """Compare encoded column-sparse matrices without dense materialization."""
    if not isinstance(left, list) or not isinstance(right, list) or len(left) != len(right):
        raise ValueError("sparse matrix column dimensions differ")
    error_sq = 0.0
    scale_sq = 0.0
    exact = True
    for col_index, (left_col, right_col) in enumerate(zip(left, right)):
        if not isinstance(left_col, list) or not isinstance(right_col, list):
            raise ValueError(f"column {col_index} is not a list")
        a, b = {}, {}
        for entry in left_col:
            if not isinstance(entry, list) or len(entry) != 3:
                raise ValueError(f"malformed left sparse entry in column {col_index}")
            row = int(entry[0])
            value = complex(float(entry[1]), float(entry[2]))
            if not (math.isfinite(value.real) and math.isfinite(value.imag)) or row in a:
                raise ValueError(f"invalid or duplicate left entry in column {col_index}")
            a[row] = value
        for entry in right_col:
            if not isinstance(entry, list) or len(entry) != 3:
                raise ValueError(f"malformed right sparse entry in column {col_index}")
            row = int(entry[0])
            value = complex(float(entry[1]), float(entry[2]))
            if not (math.isfinite(value.real) and math.isfinite(value.imag)) or row in b:
                raise ValueError(f"invalid or duplicate right entry in column {col_index}")
            b[row] = value
        if any(row < 0 or row >= basis_dimension for row in set(a) | set(b)):
            raise ValueError("sparse row index is outside the declared closure basis")
        for row in set(a) | set(b):
            av, bv = a.get(row, 0j), b.get(row, 0j)
            delta = av - bv
            error_sq += delta.real * delta.real + delta.imag * delta.imag
            scale_sq += max(
                av.real * av.real + av.imag * av.imag,
                bv.real * bv.real + bv.imag * bv.imag,
            )
            if av != bv:
                exact = False
    absolute = math.sqrt(error_sq)
    relative = absolute / max(math.sqrt(scale_sq), 1e-30)
    return {"exactly_equal": exact, "absolute_l2_residual": absolute, "relative_l2_residual": relative}


def _ledger_projection(report):
    ledger = report.get("lower_level_pruning_ledger") or {}
    return {
        "event_count": ledger.get("event_count"),
        "events": ledger.get("events"),
        "by_site": ledger.get("by_site"),
        "by_site_and_kind": ledger.get("by_site_and_kind"),
        "coverage": ledger.get("coverage"),
    }


def _validate_mode(report, expected_mode, key):
    params = report
    checks = {
        "status": params.get("status") == "PASS_ADJOINT_CONSISTENCY",
        "algebraic_status": params.get("algebraic_status") == "PASS_ADJOINT_CONSISTENCY",
        "source_node": params.get("source_nodes") == [key[1]],
        "seed_Ks": params.get("seed_Ks") == [key[2]],
        "Jmax2": params.get("Jmax2") == key[0],
        "cache_mode": params.get("cache_mode") == expected_mode,
        "pruning_mode": params.get("pruning_mode") == "production",
        "unbounded_closure": params.get("unbounded_closure") is True and params.get("max_states") is None and params.get("max_depth") is None,
        "no_omitted_closure_support": params.get("omitted_new_support_keys") == 0,
        "lower_level_coverage_disclosed": params.get("lower_level_pruning_coverage_complete") is False,
        "basis_hash": isinstance(params.get("basis_id_sha256"), str) and len(params["basis_id_sha256"]) == 64,
        "basis_keys": isinstance(params.get("basis_keys_repr"), list) and bool(params["basis_keys_repr"]),
        "coefficient_loss_ledger": isinstance(params.get("coefficient_loss_ledger"), dict),
        "coherent_node_sum_ledger": isinstance(params.get("coherent_node_sum_ledger"), list),
        "pruning_ledger": isinstance(params.get("lower_level_pruning_ledger"), dict),
    }
    for field in MATRIX_FIELDS:
        checks[field] = isinstance(params.get(field), list) and bool(params[field])
    physical = params.get("physical_status") or {}
    checks["physical_ordering_not_selected"] = physical.get("physical_hermitian_ordering_selected") is False
    checks["full_physical_HL_not_claimed"] = physical.get("full_physical_HL") is False
    diagnostics = params.get("diagnostics") or {}
    residual = diagnostics.get("direct_adjoint_relative_residual")
    checks["finite_direct_adjoint_residual"] = isinstance(residual, (int, float)) and math.isfinite(float(residual))
    return [name for name, passed in checks.items() if not passed]


def aggregate(directory: Path, jmax2_values=DEFAULT_JMAX2, tolerance=TOLERANCE):
    expected = {
        (jmax2, node, seed_k)
        for jmax2 in jmax2_values
        for node in NODES
        for seed_k in SEED_KS
    }
    found = {}
    errors = []
    paths = sorted(directory.glob("W5_CACHE_MATRIX_J*_NODE*_K*_MODE*_RUN*.json"))
    for path in paths:
        report = _load(path)
        try:
            key = (int(report["Jmax2"]), int(report["source_nodes"][0]), int(report["seed_Ks"][0]))
            if len(report["source_nodes"]) != 1 or len(report["seed_Ks"]) != 1:
                raise ValueError("expected exactly one node and one seed K")
            mode = str(report["cache_mode"])
        except Exception as exc:
            errors.append(f"{path.name}: malformed shard identity: {exc}")
            continue
        if key not in expected or mode not in {"on", "off"}:
            errors.append(f"{path.name}: unexpected shard identity {(key, mode)}")
            continue
        slot = (key, mode)
        if slot in found:
            errors.append(f"{path.name}: duplicate shard identity {slot}")
            continue
        found[slot] = (path.name, report)
        failed = _validate_mode(report, mode, key)
        if failed:
            errors.append(f"{path.name}: failed checks: {', '.join(failed)}")

    pair_results = []
    paired_keys = set()
    for key in sorted(expected):
        on_item = found.get((key, "on"))
        off_item = found.get((key, "off"))
        if on_item is None or off_item is None:
            continue
        paired_keys.add(key)
        on_name, on = on_item
        off_name, off = off_item
        comparisons = {}
        same_basis = False
        try:
            same_basis = (
                on.get("basis_id_sha256") == off.get("basis_id_sha256")
                and on.get("basis_keys_repr") == off.get("basis_keys_repr")
            )
            if not same_basis:
                raise ValueError("cache-on/cache-off closure basis identity differs")
            for field in MATRIX_FIELDS:
                comparisons[field] = _sparse_residual(on[field], off[field], len(on["basis_keys_repr"]))
            matrix_ok = all(
                item["relative_l2_residual"] <= tolerance
                for item in comparisons.values()
            )
        except Exception as exc:
            matrix_ok = False
            errors.append(f"{key}: matrix comparison failed: {type(exc).__name__}: {exc}")
        coefficient_ledger_equal = on.get("coefficient_loss_ledger") == off.get("coefficient_loss_ledger")
        coherent_ledger_equal = on.get("coherent_node_sum_ledger") == off.get("coherent_node_sum_ledger")
        pruning_ledger_equal = _ledger_projection(on) == _ledger_projection(off)
        if not matrix_ok:
            errors.append(f"{key}: cache-on/cache-off sparse matrix residual exceeds {tolerance}")
        if not coefficient_ledger_equal:
            errors.append(f"{key}: coefficient loss ledgers differ")
        if not coherent_ledger_equal:
            errors.append(f"{key}: coherent node-sum ledgers differ")
        if not pruning_ledger_equal:
            errors.append(f"{key}: ordered runtime pruning traces/coverage differ")
        pair_results.append({
            "Jmax2": key[0], "node": key[1], "seed_K": key[2],
            "cache_on_artifact": on_name, "cache_off_artifact": off_name,
            "same_basis": same_basis,
            "matrix_comparisons": comparisons,
            "coefficient_loss_ledger_equal": coefficient_ledger_equal,
            "coherent_node_sum_ledger_equal": coherent_ledger_equal,
            "ordered_pruning_trace_and_coverage_equal": pruning_ledger_equal,
            "passed": matrix_ok and coefficient_ledger_equal and coherent_ledger_equal and pruning_ledger_equal,
        })

    missing = []
    for key in sorted(expected):
        for mode in ("on", "off"):
            if (key, mode) not in found:
                missing.append([key[0], key[1], key[2], mode])
    if missing:
        errors.append(f"missing {len(missing)} cache-mode shards")
    passed = not errors and len(pair_results) == len(expected)
    return {
        "schema": "bqg-w5-cache-equivalence-pair-aggregate-v2",
        "status": "PASS_ALL_CACHE_PAIRS" if passed else "INCOMPLETE_OR_FAILED_CACHE_PAIRS",
        "passed": passed,
        "comparison_performed": bool(pair_results),
        "expected_pair_count": len(expected),
        "completed_pair_count": len(pair_results),
        "expected_mode_shard_count": 2 * len(expected),
        "valid_mode_shard_count": len(found),
        "Jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "seed_Ks": list(SEED_KS),
        "tolerance": tolerance,
        "missing_shards": missing,
        "errors": errors,
        "pairs": pair_results,
        "scope": (
            "Finite unbounded-within-fixed-Jmax2 W5 cache-on/cache-off comparison. "
            "Compares complete sparse forward and direct-adjoint matrices, basis identity, "
            "coefficient-loss and coherent node-sum ledgers, and ordered instrumented "
            "pruning traces. Static coverage is limited to the declared instrumented "
            "closure; this is not a global no-pruning proof, all-spin theorem, continuum "
            "limit, or physical self-adjointness result."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--Jmax2", type=int, nargs="+", default=list(DEFAULT_JMAX2))
    parser.add_argument("--tolerance", type=float, default=TOLERANCE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        raise SystemExit("--tolerance must be finite and nonnegative")
    result = aggregate(args.directory, args.Jmax2, args.tolerance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "completed_pair_count": result["completed_pair_count"],
        "expected_pair_count": result["expected_pair_count"],
        "valid_mode_shard_count": result["valid_mode_shard_count"],
        "expected_mode_shard_count": result["expected_mode_shard_count"],
        "missing_shards": len(result["missing_shards"]),
        "errors": result["errors"][:10],
    }, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
