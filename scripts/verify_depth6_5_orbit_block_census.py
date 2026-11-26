#!/usr/bin/env python3
"""Fail-closed checks for the recovered [5] orbit-block census artifact."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "verification_results" / "DEPTH6_5_ORBIT_BLOCK_CENSUS_20261126.json"
MANIFEST = ROOT / "research" / "depth6_recovered_source" / "source_manifest.json"


def require(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    errors: list[str] = []
    try:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(json.dumps({"valid": False, "errors": [str(exc)]}, indent=2))
        return 1

    require(errors, data.get("status") == "STRUCTURAL_BLOCK_CENSUS_ONLY", "census status mismatch")
    require(errors, data.get("rank_claim") is None, "census must not claim rank")
    require(errors, data.get("gauss_admissible_spin_assignments") == 264962, "shell assignment count mismatch")
    require(errors, data.get("s5_spin_orbits") == 2757, "S5 orbit count mismatch")
    require(errors, data.get("active_trivial_s5_blocks") == 2757, "active [5] block count mismatch")
    require(errors, data.get("summed_trivial_selector_multiplicity") == 27227, "summed selector multiplicity mismatch")
    require(errors, data.get("max_block_multiplicity") == 180, "max block multiplicity mismatch")
    require(errors, data.get("max_selector_coordinate_dimension") == 3240, "max selector workspace mismatch")
    require(errors, data.get("nonvacuum_active_blocks") == 2756, "nonvacuum block count mismatch")
    require(errors, data.get("nonvacuum_selector_multiplicity") == 27226, "nonvacuum dimension mismatch")
    require(errors, data.get("vacuum_block", {}).get("representative") == [0] * 10, "vacuum block representative mismatch")
    require(errors, data.get("vacuum_block", {}).get("selector_multiplicity") == 1, "vacuum block must be one dimensional")
    require(errors, "floating representation traces rounded" in data.get("multiplicity_method", ""), "trace-rounding limitation missing")
    require(errors, "certified per-block ranks remain open" in data.get("proof_limit", ""), "rank limitation missing")
    require(errors, data.get("source", {}).get("bundle_sha256") == manifest.get("bundle_sha256"), "bundle provenance mismatch")
    require(errors, data.get("source", {}).get("module_sha256") == manifest.get("module_sha256"), "module hash map mismatch")

    blocks = data.get("blocks", [])
    require(errors, len(blocks) == 2757, "per-orbit block list length mismatch")
    ids = [b.get("orbit_id") for b in blocks]
    require(errors, len(set(ids)) == len(ids), "duplicate orbit id")
    require(errors, sorted(ids) == list(range(2757)), "orbit ids are not the complete 0..2756 range")
    require(errors, sum(b.get("selector_multiplicity", 0) for b in blocks) == 27227, "block multiplicities do not sum to 27227")
    require(errors, sum(b.get("selector_multiplicity", 0) for b in blocks if b.get("representative") != [0] * 10) == 27226, "nonvacuum block multiplicities do not sum to 27226")
    require(errors, max((b.get("selector_multiplicity", 0) for b in blocks), default=0) == 180, "block maximum mismatch")
    require(errors, max((b.get("selector_coordinate_dimension", 0) for b in blocks), default=0) == 3240, "coordinate maximum mismatch")
    require(errors, all(len(b.get("representative", [])) == 10 for b in blocks), "malformed spin representative")
    require(errors, all(all(isinstance(s, int) and 0 <= s <= 7 for s in b["representative"]) for b in blocks), "spin digit outside frozen wall")

    root = ROOT / "research" / "depth6_recovered_source"
    actual_hashes = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in manifest.get("module_sha256", {})
    }
    require(errors, actual_hashes == manifest.get("module_sha256"), "recovered source module bytes do not match manifest")
    require(errors, len(manifest.get("module_sha256", {})) == 13, "recovered source bundle module count mismatch")

    result = {
        "valid": not errors,
        "artifact": str(ARTIFACT.relative_to(ROOT)),
        "status": data.get("status"),
        "blocks": len(blocks),
        "selector_dimension": data.get("summed_trivial_selector_multiplicity"),
        "rank_claim": data.get("rank_claim"),
        "errors": errors,
    }
    print(json.dumps(result, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

