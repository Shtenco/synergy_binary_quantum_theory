#!/usr/bin/env python3
"""Aggregate paired finite W5 production/shadow_no_prune sparse matrices.

A complete aggregate is explicitly NOT a global pruning proof: only the
instrumented threshold sites are bypassed by the shadow mode.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

NODES = tuple(range(5))
JMAX2_VALUES = (10, 12, 14)
SEED_KS = [0, 2, 4, 6, 8]
MATRIX_FIELDS = ("forward_B_columns_sparse", "direct_Bdag_columns_sparse")


def _load(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"{path.name}: invalid JSON: {type(exc).__name__}: {exc}") from exc


def _decode_columns(report, field):
    basis = report.get("basis_keys_repr")
    columns = report.get(field)
    if not isinstance(basis, list) or not isinstance(columns, list) or len(columns) != len(basis):
        raise ValueError(f"{field}: missing or inconsistent basis/column dimensions")
    decoded = {}
    for column_index, entries in enumerate(columns):
        if not isinstance(entries, list):
            raise ValueError(f"{field}: column {column_index} is not a list")
        input_key = basis[column_index]
        values = {}
        for entry in entries:
            if not isinstance(entry, list) or len(entry) != 3:
                raise ValueError(f"{field}: malformed sparse entry in column {column_index}")
            row, real, imag = int(entry[0]), float(entry[1]), float(entry[2])
            if row < 0 or row >= len(basis) or row in values:
                raise ValueError(f"{field}: invalid/duplicate row index in column {column_index}")
            value = complex(real, imag)
            if not (math.isfinite(value.real) and math.isfinite(value.imag)):
                raise ValueError(f"{field}: non-finite coefficient")
            if value != 0j:
                values[basis[row]] = value
        decoded[input_key] = values
    return decoded


def _compare_columns(prod, shadow):
    common_inputs = sorted(set(prod) & set(shadow))
    prod_only_inputs = sorted(set(prod) - set(shadow))
    shadow_only_inputs = sorted(set(shadow) - set(prod))
    error_sq = 0.0
    scale_sq = 0.0
    exact = True
    differing_columns = 0
    support_mismatch_entries = 0
    for input_key in common_inputs:
        a, b = prod[input_key], shadow[input_key]
        if set(a) != set(b):
            differing_columns += 1
            support_mismatch_entries += len(set(a) ^ set(b))
            exact = False
        for output_key in set(a) | set(b):
            av, bv = a.get(output_key, 0j), b.get(output_key, 0j)
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
    return {
        "common_input_columns": len(common_inputs),
        "production_only_input_columns": len(prod_only_inputs),
        "shadow_only_input_columns": len(shadow_only_inputs),
        "differing_common_columns": differing_columns,
        "support_mismatch_entries": support_mismatch_entries,
        "exactly_equal_on_common_domain": exact and not prod_only_inputs and not shadow_only_inputs,
        "absolute_l2_residual_on_common_domain": absolute,
        "relative_l2_residual_on_common_domain": relative,
    }


def _validate(report, mode, key):
    status = "PASS_ADJOINT_CONSISTENCY" if mode == "production" else "INCOMPLETE_NO_PRUNING_COVERAGE"
    basis_keys = report.get("basis_keys_repr")
    basis_digest = (
        hashlib.sha256("\n".join(basis_keys).encode("utf-8")).hexdigest()
        if isinstance(basis_keys, list) and all(isinstance(item, str) for item in basis_keys)
        else None
    )
    checks = {
        "status": report.get("status") == status,
        "algebraic_status": report.get("algebraic_status") == "PASS_ADJOINT_CONSISTENCY",
        "source_node": report.get("source_nodes") == [key[1]],
        "all_seed_Ks": report.get("seed_Ks") == SEED_KS,
        "Jmax2": report.get("Jmax2") == key[0],
        "pruning_mode": report.get("pruning_mode") == mode,
        "unbounded_closure": report.get("unbounded_closure") is True and report.get("max_states") is None and report.get("max_depth") is None,
        "no_omitted_closure_support": report.get("omitted_new_support_keys") == 0,
        "lower_level_coverage_disclosed": report.get("lower_level_pruning_coverage_complete") is False,
        "basis_hash": isinstance(report.get("basis_id_sha256"), str) and report.get("basis_id_sha256") == basis_digest,
        "basis_keys": isinstance(report.get("basis_keys_repr"), list) and bool(report["basis_keys_repr"]),
        "coefficient_loss_ledger": isinstance(report.get("coefficient_loss_ledger"), dict),
        "coherent_node_sum_ledger": isinstance(report.get("coherent_node_sum_ledger"), list),
        "pruning_ledger": isinstance(report.get("lower_level_pruning_ledger"), dict),
    }
    for field in MATRIX_FIELDS:
        checks[field] = isinstance(report.get(field), list) and len(report[field]) == len(report.get("basis_keys_repr", []))
    diag = report.get("diagnostics") or {}
    residual = diag.get("direct_adjoint_relative_residual")
    checks["finite_direct_adjoint_residual"] = isinstance(residual, (int, float)) and math.isfinite(float(residual))
    return [name for name, ok in checks.items() if not ok]


def aggregate(directory: Path, jmax2_values=JMAX2_VALUES):
    expected = {(j, n) for j in jmax2_values for n in NODES}
    found, errors = {}, []
    for path in sorted(directory.glob("W5_SHADOW_MATRIX_J*_NODE*_MODE*_RUN*.json")):
        report = _load(path)
        try:
            nodes = report.get("source_nodes")
            if not isinstance(nodes, list) or len(nodes) != 1:
                raise ValueError("expected exactly one source node")
            key = (int(report["Jmax2"]), int(nodes[0]))
            mode = str(report["pruning_mode"])
        except Exception as exc:
            errors.append(f"{path.name}: malformed shard identity: {exc}")
            continue
        if key not in expected or mode not in {"production", "shadow_no_prune"}:
            errors.append(f"{path.name}: unexpected shard identity {(key, mode)}")
            continue
        slot = (key, mode)
        if slot in found:
            errors.append(f"{path.name}: duplicate shard identity {slot}")
            continue
        found[slot] = (path.name, report)
        failed = _validate(report, mode, key)
        if failed:
            errors.append(f"{path.name}: failed checks: {', '.join(failed)}")

    pairs = []
    for key in sorted(expected):
        production = found.get((key, "production"))
        shadow = found.get((key, "shadow_no_prune"))
        if production is None or shadow is None:
            continue
        prod_name, prod = production
        shadow_name, shadow = shadow
        same_basis = (
            prod.get("basis_id_sha256") == shadow.get("basis_id_sha256")
            and prod.get("basis_keys_repr") == shadow.get("basis_keys_repr")
        )
        comparisons = {}
        try:
            for field in MATRIX_FIELDS:
                comparisons[field] = _compare_columns(
                    _decode_columns(prod, field), _decode_columns(shadow, field)
                )
        except Exception as exc:
            errors.append(f"{key}: matrix comparison failed: {type(exc).__name__}: {exc}")
        pairs.append({
            "Jmax2": key[0], "node": key[1],
            "production_artifact": prod_name, "shadow_artifact": shadow_name,
            "same_closure_basis": same_basis,
            "production_closure_size": len(prod.get("basis_keys_repr", [])),
            "shadow_closure_size": len(shadow.get("basis_keys_repr", [])),
            "matrix_comparisons": comparisons,
            "coefficient_loss_ledgers_equal": prod.get("coefficient_loss_ledger") == shadow.get("coefficient_loss_ledger"),
            "coherent_node_sum_ledgers_equal": prod.get("coherent_node_sum_ledger") == shadow.get("coherent_node_sum_ledger"),
            "production_pruning_events": (prod.get("lower_level_pruning_ledger") or {}).get("event_count"),
            "shadow_pruning_events": (shadow.get("lower_level_pruning_ledger") or {}).get("event_count"),
            "physical_status": prod.get("physical_status"),
            "claim_boundary": "A finite production-vs-instrumented-shadow comparison only; uninstrumented pruning sites remain open.",
        })

    missing = []
    for key in sorted(expected):
        for mode in ("production", "shadow_no_prune"):
            if (key, mode) not in found:
                missing.append([key[0], key[1], mode])
    if missing:
        errors.append(f"missing {len(missing)} node/mode shards")
    complete = not errors and len(pairs) == len(expected)
    return {
        "schema": "bqg-w5-shadow-no-prune-pair-aggregate-v2",
        "status": "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED" if complete else "INCOMPLETE_OR_FAILED_SHADOW_SHARDS",
        "passed": False,
        "finite_matrix_comparison_complete": complete,
        "comparison_performed": bool(pairs),
        "expected_pair_count": len(expected),
        "completed_pair_count": len(pairs),
        "expected_mode_shard_count": 2 * len(expected),
        "valid_mode_shard_count": len(found),
        "Jmax2_values": list(jmax2_values),
        "nodes": list(NODES),
        "missing_shards": missing,
        "errors": errors,
        "pairs": pairs,
        "claim_boundary": (
            "Even a complete finite production/shadow matrix comparison does not close "
            "global pruning coverage. Uninstrumented threshold sites, spin-wall/state-depth "
            "caps outside this run, regulator convergence, all-spin domains and continuum "
            "physicality remain open."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--Jmax2", type=int, nargs="+", default=list(JMAX2_VALUES))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.directory, args.Jmax2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "finite_matrix_comparison_complete": result["finite_matrix_comparison_complete"],
        "completed_pair_count": result["completed_pair_count"],
        "expected_pair_count": result["expected_pair_count"],
        "valid_mode_shard_count": result["valid_mode_shard_count"],
        "expected_mode_shard_count": result["expected_mode_shard_count"],
        "missing_shards": len(result["missing_shards"]),
        "errors": result["errors"][:10],
    }, indent=2))
    return 0 if result["finite_matrix_comparison_complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
