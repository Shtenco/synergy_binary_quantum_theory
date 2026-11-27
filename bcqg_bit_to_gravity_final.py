#!/usr/bin/env python3
"""Canonical binary quantum geometry -> continuum-GR structural aggregator.

A successful ``--strict`` run certifies the declared internal structural
candidate package only.  Physicalization is tracked independently in
``physicalization_gates.json`` and is reported in the same output so that a
structural GREEN certificate cannot be mistaken for a solved physical
projector, interacting graviton kernel or experimental confirmation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from bcqg_observer_smoothing_unified import run as run_bits
from bcqg_global_manifold_gate import run as run_manifold
from bcqg_quantum_hda_killer import run as run_hda

ROOT = Path(__file__).resolve().parent
CORE_ACCEPTED = {"proved", "tested_finite", "conditional"}
REQUIRED_PHYSICAL = {
    "PHYSICAL_PROJECTOR_HISTORY",
    "CONNECTED_INTERBLOCK_HISTORY",
    "PHYSICAL_TT_KERNEL",
    "IR_SIX_VECTOR",
    "COMMON_SCALE_CALIBRATION",
    "DYNAMICAL_MAXWELL_KERNEL",
    "PHYSICAL_BACKGROUND_COSMOLOGY",
    "PHYSICAL_SCALAR_COSMOLOGY",
    "LENSING_DYNAMICS_CLOSURE",
}


def run():
    bits = run_bits()
    manifold = run_manifold(refinements=2)
    hda = run_hda()
    ledger = json.loads((ROOT / "theory_gates.json").read_text(encoding="utf-8"))
    physical = json.loads((ROOT / "physicalization_gates.json").read_text(encoding="utf-8"))
    gates = ledger.get("gates", [])

    bits_pass = bool(bits.get("candidate_all_passed", bits.get("all_passed", False)))
    manifold_pass = bool(manifold.get("passed", False))
    hda_prereq_pass = bool(hda.get("regression_passed", False))
    route_hda_pass = bool(hda.get("route_sector_HDA_principal_symbol_closed", False))

    core = [g for g in gates if g.get("closure_role") == "core"]
    extensions = [g for g in gates if g.get("closure_role") == "extension"]
    experiments = [g for g in gates if g.get("closure_role") == "experiment"]
    negative_controls = [
        {"id": g["id"], "status": g.get("status"), "claim": g.get("claim"), "blocking": False}
        for g in gates if g.get("closure_role") == "negative_control"
    ]

    core_rows = {
        g["id"]: {
            "status": g.get("status"),
            "accepted": g.get("status") in CORE_ACCEPTED,
        }
        for g in core
    }
    ledger_structural_closed = bool(core_rows) and all(x["accepted"] for x in core_rows.values())

    runtime_chain = {
        "binary_route_geometrogenesis": bits_pass,
        "recursive_selected_PL_3manifold": manifold_pass,
        "HDA_prerequisites_and_factorization_no_go": hda_prereq_pass,
        "route_sector_HDA_principal_symbol": route_hda_pass,
    }
    runtime_passed = all(runtime_chain.values())
    structural_closed = bool(runtime_passed and ledger_structural_closed)

    pgates = {
        g["id"]: g
        for g in physical.get("gates", [])
        if isinstance(g, dict) and isinstance(g.get("id"), str)
    }

    def frozen(gate_id: str) -> bool:
        return gate_id in pgates and pgates[gate_id].get("status") == "frozen"

    physical_statuses = {
        gate_id: pgates.get(gate_id, {}).get("status")
        for gate_id in sorted(REQUIRED_PHYSICAL)
    }
    projector_closed = frozen("PHYSICAL_PROJECTOR_HISTORY")
    connected_history_closed = frozen("CONNECTED_INTERBLOCK_HISTORY")
    tt_kernel_frozen = frozen("PHYSICAL_TT_KERNEL")
    six_vector_frozen = frozen("IR_SIX_VECTOR")
    common_scale_calibrated = frozen("COMMON_SCALE_CALIBRATION")
    physicalization_complete = all(frozen(gate_id) for gate_id in REQUIRED_PHYSICAL)

    return {
        "status": (
            "binary quantum geometry -> continuum-GR candidate; "
            "theory-specific physicalization remains fail-closed and separately tracked"
        ),
        "schema_version": ledger.get("schema_version"),
        "candidate_framework": True,
        "structural_candidate_closed": structural_closed,
        "core_candidate_architecture_closed": structural_closed,
        # Backward-compatible alias retained for historical consumers.  Its
        # scope is now machine-readable and deliberately narrower than a
        # statement that physical gravity is complete.
        "core_theory_closed": structural_closed,
        "core_theory_closed_scope": "legacy alias: structural internal candidate only",
        "physicalization": {
            "schema_version": physical.get("schema_version"),
            "physical_projector_history_closed": projector_closed,
            "connected_interblock_history_closed": connected_history_closed,
            "physical_tt_kernel_frozen": tt_kernel_frozen,
            "ir_six_wilson_vector_frozen": six_vector_frozen,
            "common_physical_scale_calibrated": common_scale_calibrated,
            "physicalization_complete": physicalization_complete,
            "required_gate_statuses": physical_statuses,
        },
        "experimentally_confirmed_theory_of_nature": False,
        "experimental_confirmation_is_separate_from_internal_structure_and_physicalization": True,
        "runtime_chain": runtime_chain,
        "core_gate_count": len(core_rows),
        "registered_core_gates": core_rows,
        "negative_controls": negative_controls,
        "non_blocking_extensions": [
            {"id": g["id"], "status": g["status"], "claim": g["claim"]}
            for g in extensions
        ],
        "experimental_tests": [
            {"id": g["id"], "status": g["status"], "claim": g["claim"]}
            for g in experiments
        ],
        "selected_numeric_certificates": {
            "exact_dimension_fixed_point": 3.0,
            "gauss_singlet_weight": 2.0 / 9.0,
            "logical_oriented_volume": "sqrt(3)/4",
            "three_node_joint_HDA_exponent": 1.0064429343878083,
            "euclidean_HH_safe_Jmax": 2.5,
            "declared_lorentzian_HH_safe_Jmax": 6.5,
            "regge_L6_frozen_prediction": 0.11876923193907167,
            "regge_L6_observed": 0.11876075461190198,
            "quartic_TT_dimension": 6,
            "six_wilson_extractor_det": "1/699840000",
            "higher_shell_lambda_min": 10.635759878291307,
            "higher_shell_lambda_max": 15.059927665966466,
        },
        "claim_boundary": {
            "structural": (
                "A passing structural result certifies the declared exact, finite-tested and explicitly conditional candidate architecture."
            ),
            "physicalization": (
                "Finite master-projector and relational/source positive controls do not close the theory-specific refinement/rigging map, connected interblock history, physical TT 1PI kernel, IR six-vector or one-scale calibration."
            ),
            "spectral_variable": (
                "A Hamiltonian-constraint spectral parameter z is not physical omega without an independently derived physical-history/time construction."
            ),
            "experiment": (
                "External observations remain required to establish whether the candidate describes nature."
            ),
        },
        "canonical_status": "THEORY_STATUS.md",
        "canonical_package": "CANONICAL_THEORY_PACKAGE.md",
        "physicalization_ledger": "physicalization_gates.json",
        "experimental_layer": "PREDICTIONS_AND_EXPERIMENTAL_TESTS.md",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    ap.add_argument(
        "--strict",
        action="store_true",
        help="fail unless the declared structural candidate package closes",
    )
    ap.add_argument(
        "--strict-physicalization",
        action="store_true",
        help="fail unless all required physicalization gates are frozen (currently expected to fail)",
    )
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2, ensure_ascii=False)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")

    structural_ok = bool(out["structural_candidate_closed"])
    physical_ok = bool(out["physicalization"]["physicalization_complete"])
    if args.strict and not structural_ok:
        return 2
    if args.strict_physicalization and not physical_ok:
        return 3
    return 0 if structural_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
