#!/usr/bin/env python3
"""Offline consistency gate for the LeMaitre–Perche / reduced-BQG cone comparison.

PASS means only that the two published finite reference artefacts are numerically
consistent after a declared dimensionless normalization. It is not a physical
geometrogenesis or full-BQG validation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

# Pinned external reference: Data/planner_neighbor_metric.json in
# AntiPhoton47/Learning-Spacetime-from-Quantum-Info.
LEMAITRE = {
    "run_id": "20261017_141731",
    "master_seed": 20260727,
    "agent_seeds": [399, 410, 982, 979, 839],
    "environment_seeds": [944, 303, 238, 477, 151],
    "evaluation_seed": 12345,
    "evaluation_blob_shas": {
        "agent_0": "64b6fd913c0f6161defd87656c37c8eec2f61b5c",
        "agent_1": "e6fddfa6f6b0bc5988162df8d84cdbddd75fc94f",
        "agent_2": "26a0e0f0ea81bed9d0d8cf0fe6b8d360aae98056",
        "agent_3": "577dab5284b3ce411150966043ab1ed2b0ca6f77",
        "agent_4": "2030e0942839b904771f6c230003f8305f2a2633",
    },
    "g_tt": -0.9998780957794325,
    "g_tx": 7.786653464611928e-08,
    "g_xx": 1.0,
    "fit_point_count": 505,
    "fit_mode": "legacy_unweighted_quadratic",
    "fit_rmse": 0.11668354900549102,
    "source_grid": "planner_softmax_average",
    "metric_percentage": 0.7,
    "radius": 1,
    "notebook_blob_sha": "aa72e00fac1ae24fb33502630e1528e9382091aa",
    "metric_blob_sha": "8d22049c794b0775e867061a60a9d710709c76f1",
    "manifest_blob_sha": "dbd6e969140d8976467f349b2ea89509a01a1a4a",
}

# Reduced BQG free-fermion reference documented in
# BQG_CAUSAL_CONE_AND_W5_SPECTRAL_CERTIFICATE.md and
# scripts/bqg_sharp_analytic_causal_cone_gate.py.
def bqg_front_velocity(mass: float) -> float:
    m = abs(float(mass))
    return math.sqrt(4.0 * m * m + 9.0) - 2.0 * m


def bqg_critical_strip(mass: float) -> float:
    return math.asinh(2.0 * abs(float(mass)) / 3.0)


def lematre_null_slopes() -> tuple[float, float]:
    """Solve g_xx v^2 + 2 g_tx v + g_tt = 0 for v=dx/dt."""
    a = LEMAITRE["g_xx"]
    b = 2.0 * LEMAITRE["g_tx"]
    c = LEMAITRE["g_tt"]
    disc = b * b - 4.0 * a * c
    if not (a > 0.0 and disc > 0.0):
        raise ValueError("reference tensor does not have two real null slopes")
    roots = ((-b - math.sqrt(disc)) / (2.0 * a),
             (-b + math.sqrt(disc)) / (2.0 * a))
    return tuple(sorted(roots))


def run() -> dict:
    checks: list[dict] = []
    def check(name: str, ok: bool, value=None, detail: str = "") -> None:
        checks.append({"name": name, "passed": bool(ok), "value": value, "detail": detail})

    gtt, gtx, gxx = LEMAITRE["g_tt"], LEMAITRE["g_tx"], LEMAITRE["g_xx"]
    det = gtt * gxx - gtx * gtx
    check("Lorentzian_signature", gxx > 0.0 and det < 0.0,
          {"determinant": det, "g_tt": gtt, "g_tx": gtx, "g_xx": gxx},
          "One negative and one positive eigenvalue in 1+1 dimensions.")

    slopes = lematre_null_slopes()
    check("two_real_null_branches", slopes[0] < 0.0 < slopes[1],
          list(slopes), "Roots recomputed from the stored metric tensor.")

    rounded_paper = (-0.99994, 0.99994)
    slope_errors = [abs(slopes[i] - rounded_paper[i]) for i in range(2)]
    check("agreement_with_rounded_paper_slopes", max(slope_errors) < 1.0e-4,
          {"computed": list(slopes), "paper_rounded": list(rounded_paper),
           "absolute_errors": slope_errors},
          "Tolerance accounts for the paper's rounded reported slopes.")

    v0 = bqg_front_velocity(0.0)
    c_normalized = v0 / 3.0
    check("BQG_reduced_massless_front", abs(v0 - 3.0) < 1.0e-13,
          {"v_front_m0": v0, "normalization_speed": 3.0},
          "Exact result for the reduced free-fermion model only.")
    check("dimensionless_cone_match", abs(c_normalized - 1.0) < 1.0e-13
          and max(abs(abs(v) - c_normalized) for v in slopes) < 1.0e-4,
          {"BQG_reduced_c_over_v0": c_normalized,
           "LeMaitre_absolute_slopes": [abs(v) for v in slopes]},
          "Comparison uses each model's dimensionless massless speed; no physical unit calibration is implied.")

    massive_rows = []
    for mass in (0.2, 0.5, 1.0, 2.0, 5.0):
        mu = bqg_critical_strip(mass)
        v = bqg_front_velocity(mass)
        identity_rhs = 3.0 * math.exp(-mu)
        massive_rows.append({
            "mass": mass,
            "critical_strip": mu,
            "correlation_length": 1.0 / mu,
            "front_velocity": v,
            "front_over_massless_front": v / 3.0,
            "identity_rhs": identity_rhs,
            "identity_error": abs(v - identity_rhs),
        })
    max_identity_error = max(row["identity_error"] for row in massive_rows)
    check("BQG_gap_front_identity", max_identity_error < 1.0e-12,
          {"rows": massive_rows, "max_abs_error": max_identity_error},
          "Checks v_front=3 exp(-mu_c), equivalent to the analytic strip/correlation-length relation.")
    check("mass_suppresses_reduced_front", all(
        0.0 < row["front_over_massless_front"] < 1.0 for row in massive_rows
    ), [row["front_over_massless_front"] for row in massive_rows],
          "Massive rows are an internal reduced-model control, not a fit to the external paper.")

    all_passed = all(item["passed"] for item in checks)
    return {
        "schema_version": 1,
        "status": "tested_finite_reference_only" if all_passed else "failed_reference_consistency",
        "passed": all_passed,
        "scientific_scope": {
            "lematre_perche": "QMI task with a massless scalar field on a pre-existing Minkowski background; 1+1 reconstruction slice.",
            "bqg": "Reduced free-fermion propagation model; not the full physical BQG theory.",
            "comparison": "Dimensionless null-slope consistency only; no shared microscopic dynamics is established.",
            "physical_maxwell_kernel": "open",
            "physical_tt_kernel": "open",
            "physical_relational_time": "open",
            "experimental_confirmation": False,
        },
        "external_reference": LEMAITRE,
        "lematre_metric": {
            "tensor": [[gtt, gtx], [gtx, gxx]],
            "determinant": det,
            "null_slopes_dx_dt": list(slopes),
        },
        "bqg_reduced_reference": {
            "dispersion": "E_±(k)=±sqrt(4m^2+9 cos^2(k))",
            "massless_front": v0,
            "normalization_speed": 3.0,
            "dimensionless_massless_front": c_normalized,
            "massive_controls": massive_rows,
        },
        "checks": checks,
        "limitations": [
            "This script does not download or execute the external notebook.",
            "It does not reproduce the five-agent training or recompute the QMI landscape.",
            "The external metric was extracted by a neighbour-ridge mask and an unweighted quadratic fit with g_xx fixed to 1.",
            "The external agent receives physical (x,t) coordinates and a coordinate decoder target.",
            "The external objective uses a normalized-power QMI transform and reward shaping.",
            "Mutual information includes spacelike correlations and is not a pure signalling observable.",
            "The BQG velocity 3 is a reduced-model unit; division by 3 is dimensionless normalization, not physical calibration.",
            "No BQG physical photon kernel, full Lorentzian continuum limit, or common matter/gravity cone is proved here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
