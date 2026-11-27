#!/usr/bin/env python3
"""Counterexamples showing why hard amplitude pruning cannot define an exact linear W5 map."""
from __future__ import annotations
import json

def prune(z: complex, tol: float) -> complex:
    return z if abs(z) > tol else 0j

def prune_after_each(terms, tol):
    acc = 0j
    for z in terms:
        acc = prune(acc + z, tol)
    return acc

def main():
    tol = 1.0
    x = 1.5 + 0j
    lam = 0.5
    px = prune(x, tol)
    plx = prune(lam * x, tol)
    lam_px = lam * px
    homogeneity_residual = abs(plx - lam_px)

    terms = (0.75 + 0j, 0.75 + 0j, -0.40 + 0j)
    sequential = prune_after_each(terms, tol)
    final_only = prune(sum(terms), tol)

    prod_tol = 1e-11
    prod_x = 2.0 * prod_tol
    prod_lam = 0.25
    prod_residual = abs(
        prune(prod_lam * prod_x, prod_tol)
        - prod_lam * prune(prod_x, prod_tol)
    )
    result = {
        "status": "PASS_COUNTEREXAMPLES",
        "scope": "mathematical pruning contract only; no W5 amplitudes computed",
        "homogeneity": {
            "threshold": tol, "x": x.real, "lambda": lam,
            "prune_x": px.real, "prune_lambda_x": plx.real,
            "lambda_prune_x": lam_px.real,
            "absolute_residual": homogeneity_residual,
            "linear_map_contract_violated": homogeneity_residual > 0,
        },
        "accumulation_order": {
            "threshold": tol, "terms": [z.real for z in terms],
            "threshold_after_each_addition": sequential.real,
            "threshold_only_after_full_sum": final_only.real,
            "absolute_residual": abs(sequential-final_only),
            "addition_order_or_grouping_can_change_result": sequential != final_only,
        },
        "production_scale_example": {
            "threshold": prod_tol, "x": prod_x, "lambda": prod_lam,
            "absolute_residual": prod_residual,
            "linear_map_contract_violated": prod_residual > 0,
        },
    }
    assert result["homogeneity"]["linear_map_contract_violated"]
    assert result["accumulation_order"]["addition_order_or_grouping_can_change_result"]
    assert result["production_scale_example"]["linear_map_contract_violated"]
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
