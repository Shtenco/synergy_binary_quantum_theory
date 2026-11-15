#!/usr/bin/env python3
"""Transitive amplitude-pruning inventory for the direct W5 SU(2) operator path."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_pruning_coverage_inventory as BASE

ROOT_MODULE = "peter_weyl_w5_direct_all_M_su2_gate"

def _w5_closure(root=ROOT_MODULE):
    return BASE.local_import_closure(root)

# Reuse the audited AST classifier, but scope it to the direct all-M W5 path.
BASE.local_import_closure = _w5_closure
BASE.EXPECTED_MODULES = _w5_closure()
BASE.NON_PRUNING_CLASSIFICATIONS.update({
    (ROOT_MODULE, "run"): "finite_operator_acceptance_residual_thresholds",
    (ROOT_MODULE, "relative_error"): "finite_operator_residual_normalization",
    (ROOT_MODULE, "wrong_sector_fraction"): "finite_SU2_sector_weight_acceptance",
    ("peter_weyl_euclidean_sine_ordering_gate", "run"): "finite_operator_acceptance_residual_thresholds",
    ("peter_weyl_lorentzian_sine_ordered_triple_gate", "run"): "finite_operator_acceptance_residual_thresholds",
    ("peter_weyl_lorentzian_ordered_triple_gate", "run"): "finite_operator_acceptance_residual_thresholds",
    ("peter_weyl_covariant_K_composition_cached_gate", "run"): "finite_cache_equivalence_acceptance",
})

def run_inventory():
    result = BASE.run_inventory()
    result["schema"] = "bqg-w5-su2-transitive-pruning-coverage-v1"
    result["root"] = ROOT_MODULE
    result["semantics"] = (
        "Completeness means every threshold-like comparison and PRUNE decision "
        "in the transitive Python import closure of the direct W5 all-M SU(2) "
        "operator is either routed through the shared PRUNE audit/bypass path "
        "or explicitly classified as a non-pruning representation, operator, "
        "acceptance, or domain rule. This proves instrumentation coverage only; "
        "the separate production-vs-shadow matrix comparison bounds the finite "
        "operator difference."
    )
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_inventory()
    text = json.dumps(result, indent=2)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if result["complete"] else 2

if __name__ == "__main__":
    raise SystemExit(main())
