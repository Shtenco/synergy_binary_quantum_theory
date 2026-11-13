#!/usr/bin/env python3
"""Validate the fail-closed unified physicalization ledger.

GREEN means truthful status reproduction. It is compatible with a structurally closed candidate
and an explicitly open physical frontier; it is not evidence that the physical theory is complete.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "physicalization_gates.json"
STRUCTURAL_LEDGER = ROOT / "theory_gates.json"
ROLE_STATUSES = {
    "reference": {"proved", "tested_finite"},
    "positive_control": {"tested_finite"},
    "physical": {"open_physical", "frozen"},
    "experiment": {"experimental_test"},
}
REQUIRED_KEYS = {"id", "status", "closure_role", "claim", "evidence", "hard_scope"}
ID_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
EXPECTED_REQUIRED_PHYSICAL = {
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
STRUCTURAL_CORE_ACCEPTED = {"proved", "tested_finite", "conditional"}
EXPECTED_SUPPLEMENTAL_PHYSICAL = {
    "PHYSICAL_ENVIRONMENT_OVERLAP",
    "PHYSICAL_BORN_MEASURE",
    "PHYSICAL_RELATIONAL_CLOCK_RESPONSE",
    "PHYSICAL_COMMON_CAUSAL_CONE",
}


def safe_relative_path(value: str) -> bool:
    path = PurePosixPath(value)
    return bool(value) and not path.is_absolute() and ".." not in path.parts


def load_json(path: Path, errors: list[str]) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"{path.name}: {exc}")
        return {}
    if not isinstance(value, dict):
        errors.append(f"{path.name}: root must be an object")
        return {}
    return value


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    errors: list[str] = []
    data = load_json(LEDGER, errors)
    structural = load_json(STRUCTURAL_LEDGER, errors)

    if data.get("schema_version") != 2:
        errors.append("physicalization schema_version must equal 2")
    expected_policy = {role: sorted(values) for role, values in ROLE_STATUSES.items()}
    policy = data.get("status_policy")
    normalized_policy = {
        role: sorted(values) if isinstance(values, list) else values
        for role, values in policy.items()
    } if isinstance(policy, dict) else None
    if normalized_policy != expected_policy:
        errors.append("physicalization status_policy does not match verifier semantics")
    if data.get("structural_status_source") != "theory_gates.json":
        errors.append("structural_status_source must be theory_gates.json")

    gates = data.get("gates", [])
    if not isinstance(gates, list):
        errors.append("gates must be a list")
        gates = []
    rows: dict[str, dict] = {}
    counts_by_role = {role: 0 for role in ROLE_STATUSES}
    counts_by_status: dict[str, int] = {}

    for index, gate in enumerate(gates):
        label = f"gate[{index}]"
        if not isinstance(gate, dict):
            errors.append(f"{label}: gate must be an object")
            continue
        missing = REQUIRED_KEYS - set(gate)
        if missing:
            errors.append(f"{label}: missing keys {sorted(missing)}")
        gid = gate.get("id")
        if not isinstance(gid, str) or not ID_RE.fullmatch(gid):
            errors.append(f"{label}: invalid id {gid!r}")
            continue
        if gid in rows:
            errors.append(f"duplicate gate id {gid}")
            continue
        role = gate.get("closure_role")
        status = gate.get("status")
        if role not in ROLE_STATUSES:
            errors.append(f"{gid}: invalid closure_role {role!r}")
        else:
            counts_by_role[role] += 1
            if status not in ROLE_STATUSES[role]:
                errors.append(f"{gid}: illegal status {status!r} for role {role!r}")
        counts_by_status[str(status)] = counts_by_status.get(str(status), 0) + 1
        for key in ("claim", "hard_scope"):
            if not isinstance(gate.get(key), str) or not gate[key].strip():
                errors.append(f"{gid}: {key} must be non-empty text")
        evidence = gate.get("evidence", [])
        if not isinstance(evidence, list) or not evidence or not all(isinstance(x, str) for x in evidence):
            errors.append(f"{gid}: evidence must be a non-empty list of paths")
            evidence = []
        if len(evidence) != len(set(evidence)):
            errors.append(f"{gid}: duplicate evidence path")
        for relative in evidence:
            if not safe_relative_path(relative):
                errors.append(f"{gid}: unsafe evidence path {relative!r}")
            elif not (ROOT / relative).is_file():
                errors.append(f"{gid}: missing evidence file {relative}")
        rows[gid] = gate

    declared = data.get("required_physical_gates")
    if not isinstance(declared, list) or set(declared) != EXPECTED_REQUIRED_PHYSICAL:
        errors.append("required_physical_gates must equal the frozen nine-gate unified physicalization frontier")
    missing = sorted(EXPECTED_REQUIRED_PHYSICAL - set(rows))
    if missing:
        errors.append(f"missing required physical gates: {missing}")
    for gid in EXPECTED_REQUIRED_PHYSICAL & set(rows):
        if rows[gid].get("closure_role") != "physical":
            errors.append(f"{gid}: required gate must have role='physical'")
    missing_supplemental = sorted(EXPECTED_SUPPLEMENTAL_PHYSICAL - set(rows))
    if missing_supplemental:
        errors.append(f"missing supplemental physical gates: {missing_supplemental}")
    for gid in EXPECTED_SUPPLEMENTAL_PHYSICAL & set(rows):
        if rows[gid].get("closure_role") != "physical":
            errors.append(f"{gid}: supplemental gate must have role='physical'")

    if "CONSTRAINT_FESHBACH" in rows and "not physical omega" not in rows["CONSTRAINT_FESHBACH"]["hard_scope"].lower():
        errors.append("CONSTRAINT_FESHBACH must state that z is not physical omega")
    if "MASTER_PROJECTOR_FINITE" in rows:
        scope = rows["MASTER_PROJECTOR_FINITE"]["hard_scope"].lower()
        if "finite" not in scope or "rigging" not in scope:
            errors.append("MASTER_PROJECTOR_FINITE must preserve finite/rigging-limit scope")
    if "PHYSICAL_TT_KERNEL" in rows and "not promoted to the physical omega" not in rows["PHYSICAL_TT_KERNEL"]["hard_scope"].lower():
        errors.append("PHYSICAL_TT_KERNEL must reject spectral/history labels as physical omega")
    if "DYNAMICAL_MAXWELL_KERNEL" in rows and "not a physical maxwell derivation" not in rows["DYNAMICAL_MAXWELL_KERNEL"]["hard_scope"].lower():
        errors.append("DYNAMICAL_MAXWELL_KERNEL must reject kinematic U(1) as Maxwell closure")
    if "PHYSICAL_BACKGROUND_COSMOLOGY" in rows and "no bqg background dark component" not in rows["PHYSICAL_BACKGROUND_COSMOLOGY"]["hard_scope"].lower():
        errors.append("PHYSICAL_BACKGROUND_COSMOLOGY must reject fitted/background shortcut")
    if "PHYSICAL_SCALAR_COSMOLOGY" in rows and "promoted to physical dark matter" not in rows["PHYSICAL_SCALAR_COSMOLOGY"]["hard_scope"].lower():
        errors.append("PHYSICAL_SCALAR_COSMOLOGY must reject TT/eigenvalue dark-matter shortcut")
    if "LENSING_DYNAMICS_CLOSURE" in rows and "prerequisites only" not in rows["LENSING_DYNAMICS_CLOSURE"]["hard_scope"].lower():
        errors.append("LENSING_DYNAMICS_CLOSURE must keep finite references distinct from physical closure")

    def frozen(gid: str) -> bool:
        return gid in rows and rows[gid].get("status") == "frozen"

    declared_pairs = {
        "physical_projector_history_closed_declared": frozen("PHYSICAL_PROJECTOR_HISTORY"),
        "physical_tt_kernel_frozen_declared": frozen("PHYSICAL_TT_KERNEL"),
        "ir_six_wilson_vector_frozen_declared": frozen("IR_SIX_VECTOR"),
        "common_physical_scale_calibrated_declared": frozen("COMMON_SCALE_CALIBRATION"),
        "dynamical_maxwell_kernel_frozen_declared": frozen("DYNAMICAL_MAXWELL_KERNEL"),
        "physical_background_cosmology_closed_declared": frozen("PHYSICAL_BACKGROUND_COSMOLOGY"),
        "physical_scalar_cosmology_closed_declared": frozen("PHYSICAL_SCALAR_COSMOLOGY"),
        "lensing_dynamics_closure_closed_declared": frozen("LENSING_DYNAMICS_CLOSURE"),
        "physical_environment_overlap_closed_declared": frozen("PHYSICAL_ENVIRONMENT_OVERLAP"),
        "physical_born_measure_closed_declared": frozen("PHYSICAL_BORN_MEASURE"),
        "physical_relational_clock_response_closed_declared": frozen("PHYSICAL_RELATIONAL_CLOCK_RESPONSE"),
        "physical_common_causal_cone_closed_declared": frozen("PHYSICAL_COMMON_CAUSAL_CONE"),
    }
    for key, derived in declared_pairs.items():
        if data.get(key) is not derived:
            errors.append(f"{key}={data.get(key)!r} disagrees with derived value {derived!r}")

    physicalization_complete = all(frozen(gid) for gid in (EXPECTED_REQUIRED_PHYSICAL | EXPECTED_SUPPLEMENTAL_PHYSICAL))

    structural_gates = structural.get("gates", []) if isinstance(structural.get("gates", []), list) else []
    structural_core = [g for g in structural_gates if isinstance(g, dict) and g.get("closure_role") == "core"]
    structural_candidate_closed = bool(structural_core) and all(
        g.get("status") in STRUCTURAL_CORE_ACCEPTED for g in structural_core
    ) and structural.get("core_theory_closed_declared") is True
    if structural.get("experimentally_confirmed") is not False:
        errors.append("theory_gates.json must keep experimentally_confirmed=false")
    if data.get("experimentally_confirmed") is not False:
        errors.append("physicalization_gates.json must keep experimentally_confirmed=false")

    result = {
        "schema_version": data.get("schema_version"),
        "gate_count": len(gates),
        "counts_by_role": counts_by_role,
        "counts_by_status": counts_by_status,
        "structural_candidate_closed": bool(structural_candidate_closed),
        "physical_projector_history_closed": frozen("PHYSICAL_PROJECTOR_HISTORY"),
        "physical_tt_kernel_frozen": frozen("PHYSICAL_TT_KERNEL"),
        "ir_six_wilson_vector_frozen": frozen("IR_SIX_VECTOR"),
        "common_physical_scale_calibrated": frozen("COMMON_SCALE_CALIBRATION"),
        "dynamical_maxwell_kernel_frozen": frozen("DYNAMICAL_MAXWELL_KERNEL"),
        "physical_background_cosmology_closed": frozen("PHYSICAL_BACKGROUND_COSMOLOGY"),
        "physical_scalar_cosmology_closed": frozen("PHYSICAL_SCALAR_COSMOLOGY"),
        "lensing_dynamics_closure_closed": frozen("LENSING_DYNAMICS_CLOSURE"),
        "physical_environment_overlap": frozen("PHYSICAL_ENVIRONMENT_OVERLAP"),
        "physical_born_measure": frozen("PHYSICAL_BORN_MEASURE"),
        "physical_relational_clock_response": frozen("PHYSICAL_RELATIONAL_CLOCK_RESPONSE"),
        "physical_common_causal_cone": frozen("PHYSICAL_COMMON_CAUSAL_CONE"),
        "supplemental_physical_gate_statuses": {
            gid: rows.get(gid, {}).get("status") for gid in sorted(EXPECTED_SUPPLEMENTAL_PHYSICAL)
        },
        "physicalization_complete": physicalization_complete,
        "experimentally_confirmed": False,
        "required_physical_gate_statuses": {
            gid: rows.get(gid, {}).get("status") for gid in sorted(EXPECTED_REQUIRED_PHYSICAL)
        },
        "scientific_interpretation": (
            "A valid GREEN result certifies truthful separation of structural results and finite references "
            "from the still-open physical history, gravity, Maxwell, background/scalar cosmology, lensing "
            "and one-scale frontier. It is not experimental confirmation."
        ),
        "errors": errors,
        "valid": not errors,
    }
    text = json.dumps(result, ensure_ascii=False, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
