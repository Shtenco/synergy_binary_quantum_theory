#!/usr/bin/env python3
"""Static transitive inventory of amplitude-threshold decisions reachable from V10.

This gate is intentionally separate from runtime loss accounting. It proves that
all threshold-looking branches in the declared V10 local import closure are
either routed through the shared PRUNE instrumentation or explicitly classified
as exact representation rules, spectral operator conventions, acceptance tests,
or independent support caps.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent

EXPECTED_MODULES = {
    "bqg_pruning_audit",
    "bqg_pruning_coverage_inventory",
    "bqg_sparse_linear_accumulator",
    "bqg_v10_w5_lorentzian_adjoint_closure",
    "charged_intertwiner_recoupling_gate",
    "k5_peter_weyl_safe_hda_column",
    "peter_weyl_covariant_K_composition_cached_gate",
    "peter_weyl_covariant_K_composition_gate",
    "peter_weyl_covariant_K_leg_gate",
    "peter_weyl_covariant_K_projection_audit_gate",
    "peter_weyl_covariant_K_sine_composition_gate",
    "peter_weyl_covariant_composition_gate",
    "peter_weyl_covariant_volume_leg_gate",
    "peter_weyl_euclidean_sine_ordering_gate",
    "peter_weyl_lorentzian_K_block_gate",
    "peter_weyl_lorentzian_ordered_triple_gate",
    "peter_weyl_lorentzian_sine_ordered_triple_gate",
    "peter_weyl_q_block_gauge_audit",
    "peter_weyl_volume_cross_rep_audit",
    "peter_weyl_zeroaware_volume_migration_experiment",
}

# Function-level exclusions are deliberate and documented. They are not
# amplitude-pruning sites: each belongs to a distinct mathematical contract.
# This is the private functools cache inside audit_lru_cache(). Its cache
# values include both the numerical result and the complete ordered PRUNE trace;
# the public wrapper replays that trace on cache hits.
AUDITED_CACHE_CLASSIFICATIONS = {
    ("bqg_pruning_audit", "cached"): "audit_lru_cache_internal_value_and_trace_cache",
}

NON_PRUNING_COMPARISON_CLASSIFICATIONS = {
    ("bqg_sparse_linear_accumulator", "sparsify_matrix_columns", "tol < 0"):
        "input_tolerance_parameter_validation",
    ("bqg_sparse_linear_accumulator", "sparsify_matrix_columns", "abs(z) > tol"):
        "post_assembly_fixed_operator_mask_requires_error_bound",
}

NON_PRUNING_CLASSIFICATIONS = {
    ("bqg_v10_w5_lorentzian_adjoint_closure", "_make_spin_wall_probe"):
        "spin_wall_nonzero_branch_accounting",
    ("bqg_v10_w5_lorentzian_adjoint_closure", "probed"):
        "spin_wall_nonzero_branch_accounting",
    ("bqg_v10_w5_lorentzian_adjoint_closure", "_closure_cap_reasons"):
        "state_and_depth_support_caps",
    ("bqg_v10_w5_lorentzian_adjoint_closure", "closure_acceptance"):
        "operator_acceptance_residual_thresholds",
    ("charged_intertwiner_recoupling_gate", "allowed_charged_labels"):
        "exact_SU2_triangle_and_parity_selection",
    ("peter_weyl_covariant_K_composition_cached_gate", "restore_scalar_charge_M"):
        "norm_bounded_SU2_selection_rule_acceptance",
    ("peter_weyl_covariant_K_composition_cached_gate", "physical_acceptance_from_existing_output"):
        "physical_operator_acceptance_residuals",
    ("peter_weyl_lorentzian_sine_ordered_triple_gate", "run"):
        "operator_leakage_acceptance_residuals",
    ("peter_weyl_volume_cross_rep_audit", "one_quartet"):
        "spectral_gap_diagnostic_zero_division_guard",
    ("peter_weyl_zeroaware_volume_migration_experiment", "patch_and_clear"):
        "spectral_zero_factor_configuration_validation",
    ("peter_weyl_covariant_volume_leg_gate", "weight_by_J"):
        "exact_irrep_label_weight_classification",
    ("peter_weyl_lorentzian_ordered_triple_gate", "run"):
        "exact_scalar_sector_label_classification",
    ("peter_weyl_q_block_gauge_audit", "sqrt_abs"):
        "spectral_operator_zero_mode_convention",
    ("peter_weyl_volume_cross_rep_audit", "stable_zeroaware_volume"):
        "spectral_operator_zero_mode_convention",
    ("peter_weyl_zeroaware_volume_migration_experiment", "zeroaware_sqrt_abs"):
        "spectral_operator_zero_mode_convention",
    ("bqg_pruning_audit", "keep"):
        "shared_pruning_policy_primitive",
    ("bqg_pruning_audit", "record_threshold"):
        "shared_pruning_audit_primitive",
    ("bqg_pruning_audit", "record_merge"):
        "shared_pruning_audit_primitive",
    ("bqg_pruning_audit", "record_dropped"):
        "coherent_output_loss_ledger",
}


def local_import_closure(root="bqg_v10_w5_lorentzian_adjoint_closure"):
    pending = [root]
    modules = set()
    while pending:
        name = pending.pop()
        if name in modules:
            continue
        path = SCRIPTS / (name + ".py")
        if not path.is_file():
            continue
        modules.add(name)
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            imported = []
            if isinstance(node, ast.Import):
                imported = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported = [node.module.split(".")[0]]
            for local_name in imported:
                if (SCRIPTS / (local_name + ".py")).is_file():
                    pending.append(local_name)
    return modules


def _plain_lru_cache_decorators(function: ast.AST) -> list[str]:
    """Return non-audited lru_cache decorators regardless of function body."""
    names = [ast.unparse(item) for item in getattr(function, "decorator_list", [])]
    return [
        name for name in names
        if "lru_cache" in name and "audit_lru_cache" not in name
    ]


def _threshold_like(compare: ast.Compare) -> bool:
    source = ast.unparse(compare).lower()
    needles = (
        "abs(", "np.abs(", "amp_tol", "threshold", "tol", "cutoff",
        "epsilon", "tau", "max_states", "max_depth", "wall",
        "leakage", "zero_factor", "zero_tol",
    )
    return any(token in source for token in needles)


def _function_has_prune_instrumentation(function: ast.AST) -> bool:
    text = ast.get_source_segment(_CURRENT_SOURCE, function) or ""
    markers = (
        "PRUNE.record_threshold", "PRUNE.keep", "PRUNE.record_merge",
        "PRUNE.filter_state", "record_dropped(", "prune_state(",
        "add_state(", "add_dict(",
    )
    return any(marker in text for marker in markers)


_CURRENT_SOURCE = ""


def run_inventory() -> dict:
    modules = local_import_closure()
    module_mismatch = sorted(modules.symmetric_difference(EXPECTED_MODULES))
    sites = []
    unclassified = []
    plain_lru_caches = []
    audited_lru_caches = []

    for module in sorted(modules):
        path = SCRIPTS / (module + ".py")
        source = path.read_text(encoding="utf-8")
        global _CURRENT_SOURCE
        _CURRENT_SOURCE = source
        tree = ast.parse(source, filename=str(path))

        functions = [
            node for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]
        for function in functions:
            key = (module, function.name)
            # Inspect every decorator, even when the function contains no
            # threshold-like comparison. Otherwise an uninstrumented cache
            # can disappear from the inventory simply because its body has
            # no local amplitude comparison.

            compares = [
                node for node in ast.walk(function)
                if isinstance(node, ast.Compare) and _threshold_like(node)
            ]
            # Cache decorators are independent of threshold comparisons.
            # Inspect them before the early continue so a cached function with
            # no local comparison cannot evade the transitive cache audit.
            for name in _plain_lru_cache_decorators(function):
                cache_record = {
                    "module": module,
                    "function": function.name,
                    "decorator": name,
                }
                audited_class = AUDITED_CACHE_CLASSIFICATIONS.get(key)
                if audited_class:
                    audited_lru_caches.append({
                        **cache_record,
                        "classification": audited_class,
                    })
                else:
                    plain_lru_caches.append(cache_record)

            if not compares:
                continue

            explicit_class = NON_PRUNING_CLASSIFICATIONS.get(key)
            instrumented = _function_has_prune_instrumentation(function)

            for compare in compares:
                expression = ast.unparse(compare)
                comparison_class = NON_PRUNING_COMPARISON_CLASSIFICATIONS.get(
                    (module, function.name, expression)
                )
                if comparison_class:
                    classification = comparison_class
                elif explicit_class:
                    classification = explicit_class
                elif instrumented:
                    classification = "instrumented_amplitude_threshold"
                else:
                    classification = "UNCLASSIFIED_THRESHOLD_SITE"
                    unclassified.append({
                        "module": module,
                        "function": function.name,
                        "line": function.lineno,
                        "comparisons": [expression],
                    })

                sites.append({
                    "module": module,
                    "function": function.name,
                    "line": compare.lineno,
                    "expression": expression,
                    "classification": classification,
                })

    complete = not module_mismatch and not unclassified and not plain_lru_caches
    return {
        "schema": "bqg-v10-transitive-pruning-coverage-v1",
        "status": "PASS_TRANSITIVE_THRESHOLD_INVENTORY" if complete else "FAIL_TRANSITIVE_THRESHOLD_INVENTORY",
        "complete": bool(complete),
        "root": "bqg_v10_w5_lorentzian_adjoint_closure",
        "module_count": len(modules),
        "modules": sorted(modules),
        "expected_module_set_matches": not module_mismatch,
        "module_set_mismatch": module_mismatch,
        "threshold_site_count": len(sites),
        "threshold_sites": sites,
        "unclassified_sites": unclassified,
        "plain_lru_caches": plain_lru_caches,
        "audited_lru_caches": audited_lru_caches,
        "semantics": (
            "Completeness means every threshold-like comparison in the transitive "
            "V10 local Python import closure is either instrumented by the shared "
            "PRUNE/loss-ledger path or explicitly classified as a non-amplitude "
            "selection, spectral, acceptance, spin-wall, or support-cap rule. "
            "It does not prove zero output loss or convergence."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_inventory()
    text = json.dumps(result, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if result["complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
