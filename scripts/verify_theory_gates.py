#!/usr/bin/env python3
"""Validate the canonical v2 *structural-candidate* theory ledger.

The ledger deliberately separates three structural-package scopes:

- core: exact / finite-tested / explicitly conditional arrows that define the
  present structural candidate package;
- extension: stronger universality/generalization results that are useful but
  do not reopen that declared structural package;
- experiment: external validation or independent replication.

Physical projector/history, interacting TT-kernel and scale-calibration status
are tracked separately in ``physicalization_gates.json``.  This validator
therefore certifies structural closure only; it does not promote that result to
physical or experimental completion.
"""
from __future__ import annotations

import json
import re
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "theory_gates.json"

ROLE_STATUSES = {
    "core": {"proved", "tested_finite", "conditional"},
    "negative_control": {"no_go"},
    "extension": {"external_extension"},
    "experiment": {"experimental_test"},
}
ALL_STATUSES = set().union(*ROLE_STATUSES.values())
ID_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
REQUIRED_KEYS = {"id", "status", "closure_role", "claim", "evidence"}


def _safe_relative_path(value: str) -> bool:
    path = PurePosixPath(value)
    return bool(value) and not path.is_absolute() and ".." not in path.parts


def main() -> int:
    errors: list[str] = []

    try:
        data = json.loads(LEDGER.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(json.dumps({"valid": False, "errors": [str(exc)]}, indent=2))
        return 1

    if not isinstance(data, dict):
        errors.append("ledger root must be a JSON object")
        data = {}

    if data.get("schema_version") != 2:
        errors.append("schema_version must equal 2")

    policy = data.get("status_policy")
    expected_policy = {role: sorted(values) for role, values in ROLE_STATUSES.items()}
    if not isinstance(policy, dict):
        errors.append("status_policy must be an object")
    else:
        normalized = {}
        for role, values in policy.items():
            normalized[role] = sorted(values) if isinstance(values, list) else values
        if normalized != expected_policy:
            errors.append("status_policy does not match verifier semantics")

    gates = data.get("gates", [])
    if not isinstance(gates, list):
        errors.append("gates must be a list")
        gates = []

    ids: list[str] = []
    counts_by_role = {role: 0 for role in ROLE_STATUSES}
    counts_by_status = {status: 0 for status in sorted(ALL_STATUSES)}
    core_rows: list[dict] = []
    negative_control_rows: list[dict] = []

    for index, gate in enumerate(gates):
        label = f"gate[{index}]"
        if not isinstance(gate, dict):
            errors.append(f"{label}: gate must be an object")
            continue

        missing = REQUIRED_KEYS - set(gate)
        if missing:
            errors.append(f"{label}: missing keys {sorted(missing)}")

        gate_id = gate.get("id", "")
        if not isinstance(gate_id, str) or not ID_RE.fullmatch(gate_id):
            errors.append(f"{label}: invalid id {gate_id!r}")
        else:
            ids.append(gate_id)
            label = gate_id

        role = gate.get("closure_role")
        status = gate.get("status")
        if role not in ROLE_STATUSES:
            errors.append(f"{label}: invalid closure_role {role!r}")
        else:
            counts_by_role[role] += 1
            if status not in ROLE_STATUSES[role]:
                errors.append(
                    f"{label}: status {status!r} is not legal for closure_role {role!r}"
                )

        if status not in ALL_STATUSES:
            errors.append(f"{label}: invalid status {status!r}")
        else:
            counts_by_status[status] += 1

        claim = gate.get("claim")
        if not isinstance(claim, str) or not claim.strip():
            errors.append(f"{label}: claim must be non-empty text")

        evidence = gate.get("evidence", [])
        if not isinstance(evidence, list) or not evidence or not all(isinstance(x, str) for x in evidence):
            errors.append(f"{label}: evidence must be a non-empty list of paths")
            evidence = []

        if len(evidence) != len(set(evidence)):
            errors.append(f"{label}: duplicate evidence path")

        for relative in evidence:
            if not _safe_relative_path(relative):
                errors.append(f"{label}: unsafe evidence path {relative!r}")
                continue
            if not (ROOT / relative).is_file():
                errors.append(f"{label}: missing evidence file {relative}")

        if role == "core":
            core_rows.append({
                "id": gate_id,
                "status": status,
                "accepted_for_structural_closure": status in {"proved", "tested_finite", "conditional"},
            })
        elif role == "negative_control":
            negative_control_rows.append({
                "id": gate_id,
                "status": status,
                "blocking": False,
                "claim": claim,
            })

    if len(ids) != len(set(ids)):
        duplicates = sorted({x for x in ids if ids.count(x) > 1})
        errors.append(f"gate identifiers must be unique; duplicates={duplicates}")

    structural_candidate_closed = bool(core_rows) and all(
        row["accepted_for_structural_closure"] for row in core_rows
    )
    declared = data.get("core_theory_closed_declared")
    if not isinstance(declared, bool):
        errors.append("core_theory_closed_declared must be a Boolean")
    if declared is True and not structural_candidate_closed:
        errors.append("declared structural closure is inconsistent with core gate statuses")

    experimentally_confirmed = data.get("experimentally_confirmed")
    if experimentally_confirmed is not False:
        errors.append("experimentally_confirmed must remain false unless external evidence policy changes")

    closed_and_valid = bool(structural_candidate_closed and not errors)
    result = {
        "schema_version": data.get("schema_version"),
        "gate_count": len(gates),
        "counts_by_role": counts_by_role,
        "counts_by_status": counts_by_status,
        "core_gate_count": len(core_rows),
        "negative_control_count": len(negative_control_rows),
        "negative_controls": negative_control_rows,
        "structural_candidate_closed": closed_and_valid,
        "complete_structural_package": closed_and_valid,
        # Backward-compatible output for historical consumers.  The scope is
        # explicit so it cannot legitimately be read as physical completion.
        "core_theory_closed": closed_and_valid,
        "core_theory_closed_scope": "legacy alias: structural internal candidate only",
        "candidate_framework": True,
        "physicalization_status_source": "physicalization_gates.json",
        "experimentally_confirmed": False,
        "core_gates": core_rows,
        "errors": errors,
        "valid": not errors,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
