#!/usr/bin/env python3
"""Evidence-linked fail-closed certificate for BQG physical spectral history.

This gate exists because a finite spectral quotient and a physical BQG history
are different claims. It never trusts inline `physical_preconditions` from the
moment packet. Instead it binds independently produced artifacts through their
habitat/domain/convention/master hashes.

Closure layers:

1. finite spectral closure on the declared seed;
2. complete finite full-master/HDA certification on the same habitat;
3. refinement/rigging-map history certification;
4. optional source-dressed connected-history certification.

Only (1)+(2)+(3) can close the repository-level PHYSICAL_PROJECTOR_HISTORY
claim. Layer (4) is additionally required before connected correlator claims.

Important status rule: top-level `passed` means the refinement-compatible
PHYSICAL_PROJECTOR_HISTORY certificate passed. A finite linked stage is exposed
separately as `finite_stage_passed` and is NOT a green physical certificate.

The self-test is synthetic and checks positive hash linkage plus mismatch and
missing-refinement negative controls. It is not a BQG physical result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping

SPECTRAL_SCHEMA = "BQG_SPECTRAL_HISTORY_GRAPH_V1"
CERT_SCHEMA = "BQG_PHYSICAL_SPECTRAL_HISTORY_CERTIFICATE_V1"
REFINEMENT_SCHEMA = "BQG_PHYSICAL_REFINEMENT_HISTORY_CERTIFICATE_V1"
SOURCE_SCHEMA = "BQG_SOURCE_DRESSED_HISTORY_CERTIFICATE_V1"

LINK_FIELDS = (
    "habitat_hash",
    "domain_hash",
    "convention_hash",
    "master_pencil_hash",
)


def canonical_sha256(obj: Mapping) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _master_identity(master: Mapping) -> dict:
    return {
        "habitat_hash": str(master.get("habitat_hash", "")),
        "domain_hash": str(master.get("domain_hash", "")),
        "convention_hash": str(master.get("convention_hash", "")),
        "master_pencil_hash": str(master.get("master_pencil_hash", "")),
    }


def _spectral_identity(spectral: Mapping) -> dict:
    p = spectral.get("provenance", {})
    return {name: str(p.get(name, "")) for name in LINK_FIELDS}


def _matching_identity(expected: Mapping, got: Mapping) -> dict:
    return {
        name: bool(expected.get(name)) and str(got.get(name, "")) == str(expected.get(name, ""))
        for name in LINK_FIELDS
    }


def audit_master(master: Mapping) -> dict:
    identity = _master_identity(master)
    hda = master.get("quantum_hda_certificate_audit", {})
    checks = {
        "assembler_passed": bool(master.get("passed", False)),
        "domain_complete": bool(master.get("domain_complete", False)),
        "quantum_hda_closed": bool(master.get("quantum_hda_closed", False)),
        "matching_hda_certificate_valid": bool(hda.get("valid_for_this_master", False)),
        "physical_projector_emitted_by_full_master_assembler": bool(master.get("physical_projector_emitted", False)),
        "identity_hashes_present": all(bool(identity[k]) for k in LINK_FIELDS),
    }
    return {
        "valid": bool(all(checks.values())),
        "checks": checks,
        "identity": identity,
        "status": master.get("status"),
        "spectrum": master.get("spectrum"),
    }


def audit_spectral(spectral: Mapping, expected_identity: Mapping) -> dict:
    identity = _spectral_identity(spectral)
    link = _matching_identity(expected_identity, identity)
    checks = {
        "schema": spectral.get("schema") == SPECTRAL_SCHEMA,
        "finite_spectral_history_closed": bool(spectral.get("finite_spectral_history_closed", False)),
        "spectral_gate_did_not_self_promote_physical_history": not bool(spectral.get("physical_history_closed", False)),
        "spectral_gate_did_not_self_emit_physical_projector": not bool(spectral.get("physical_projector_emitted", False)),
        "provenance_hashes_present": all(bool(identity[k]) for k in LINK_FIELDS),
        "provenance_matches_master": all(link.values()),
    }
    return {
        "valid": bool(all(checks.values())),
        "checks": checks,
        "identity": identity,
        "hash_matches": link,
        "status": spectral.get("status"),
        "seed_label": spectral.get("seed_label"),
        "termination": spectral.get("termination"),
    }


def audit_refinement(refinement: Mapping | None, expected_identity: Mapping) -> dict:
    if refinement is None:
        return {
            "present": False,
            "valid": False,
            "checks": {
                "schema": False,
                "certificate_passed": False,
                "same_anchor_master": False,
                "refinement_family_hash_present": False,
                "low_cluster_scale_separation": False,
                "projector_converged_under_embeddings": False,
                "boundary_history_converged": False,
                "hda_residual_converged": False,
                "rank_rule_preregistered_or_heldout": False,
            },
            "reason": "no refinement/rigging-map certificate supplied",
        }

    same_anchor = (
        str(refinement.get("anchor_habitat_hash", "")) == expected_identity["habitat_hash"]
        and str(refinement.get("anchor_domain_hash", "")) == expected_identity["domain_hash"]
        and str(refinement.get("anchor_convention_hash", "")) == expected_identity["convention_hash"]
        and str(refinement.get("anchor_master_pencil_hash", "")) == expected_identity["master_pencil_hash"]
    )
    checks = {
        "schema": refinement.get("schema") == REFINEMENT_SCHEMA,
        "certificate_passed": bool(refinement.get("passed", False)),
        "same_anchor_master": bool(same_anchor),
        "refinement_family_hash_present": bool(str(refinement.get("refinement_family_hash", ""))),
        "low_cluster_scale_separation": bool(refinement.get("low_cluster_scale_separation", False)),
        "projector_converged_under_embeddings": bool(refinement.get("projector_converged_under_embeddings", False)),
        "boundary_history_converged": bool(refinement.get("boundary_history_converged", False)),
        "hda_residual_converged": bool(refinement.get("hda_residual_converged", False)),
        "rank_rule_preregistered_or_heldout": bool(refinement.get("rank_rule_preregistered_or_heldout", False)),
    }
    return {
        "present": True,
        "valid": bool(all(checks.values())),
        "checks": checks,
        "refinement_family_hash": refinement.get("refinement_family_hash"),
        "anchor": {
            "habitat_hash": refinement.get("anchor_habitat_hash"),
            "domain_hash": refinement.get("anchor_domain_hash"),
            "convention_hash": refinement.get("anchor_convention_hash"),
            "master_pencil_hash": refinement.get("anchor_master_pencil_hash"),
        },
        "reason": "matched refinement/rigging-map certificate" if all(checks.values()) else "refinement certificate incomplete or mismatched",
    }


def audit_source(source: Mapping | None, expected_identity: Mapping) -> dict:
    if source is None:
        return {
            "present": False,
            "valid": False,
            "checks": {
                "schema": False,
                "certificate_passed": False,
                "same_master": False,
                "source_operator_set_complete": False,
                "source_dressed_history_converged": False,
                "connected_W_not_Z_used": False,
            },
            "reason": "no source-dressed history certificate supplied",
        }
    identity = {name: str(source.get(name, "")) for name in LINK_FIELDS}
    link = _matching_identity(expected_identity, identity)
    checks = {
        "schema": source.get("schema") == SOURCE_SCHEMA,
        "certificate_passed": bool(source.get("passed", False)),
        "same_master": all(link.values()),
        "source_operator_set_complete": bool(source.get("source_operator_set_complete", False)),
        "source_dressed_history_converged": bool(source.get("source_dressed_history_converged", False)),
        "connected_W_not_Z_used": bool(source.get("connected_W_not_Z_used", False)),
    }
    return {
        "present": True,
        "valid": bool(all(checks.values())),
        "checks": checks,
        "identity": identity,
        "hash_matches": link,
        "reason": "matched source-dressed history certificate" if all(checks.values()) else "source certificate incomplete or mismatched",
    }


def certify(
    spectral: Mapping,
    master: Mapping,
    refinement: Mapping | None = None,
    source: Mapping | None = None,
) -> dict:
    ma = audit_master(master)
    sa = audit_spectral(spectral, ma["identity"])
    ra = audit_refinement(refinement, ma["identity"])
    soa = audit_source(source, ma["identity"])

    finite_evidence_linked = bool(ma["valid"] and sa["valid"])
    physical_projector_history_closed = bool(finite_evidence_linked and ra["valid"])
    connected_source_history_closed = bool(physical_projector_history_closed and soa["valid"])

    if connected_source_history_closed:
        status = "CONNECTED_SOURCE_DRESSED_BQG_HISTORY_CERTIFIED"
    elif physical_projector_history_closed:
        status = "PHYSICAL_BQG_PROJECTOR_HISTORY_CERTIFIED"
    elif finite_evidence_linked:
        status = "FINITE_FULL_MASTER_SPECTRAL_HISTORY_CERTIFIED_REFINEMENT_OPEN"
    else:
        status = "PHYSICAL_SPECTRAL_HISTORY_OPEN"

    return {
        "schema": CERT_SCHEMA,
        "status": status,
        "passed": physical_projector_history_closed,
        "finite_stage_passed": finite_evidence_linked,
        "artifact_hashes": {
            "spectral_result_sha256": canonical_sha256(spectral),
            "master_result_sha256": canonical_sha256(master),
            "refinement_result_sha256": canonical_sha256(refinement) if refinement is not None else None,
            "source_result_sha256": canonical_sha256(source) if source is not None else None,
        },
        "master_audit": ma,
        "spectral_audit": sa,
        "refinement_audit": ra,
        "source_audit": soa,
        "finite_full_master_spectral_history_certified": finite_evidence_linked,
        "physical_projector_history_closed": physical_projector_history_closed,
        "connected_source_history_closed": connected_source_history_closed,
        "physicalization_gate_update_allowed": physical_projector_history_closed,
        "forbidden_promotions": {
            "constraint_sigma_is_physical_time": False,
            "constraint_spectral_dimension_is_spacetime_dimension": False,
            "master_eigenvalue_is_particle_mass": False,
            "finite_history_without_refinement_is_continuum_history": False,
        },
        "claim_boundary": (
            "Finite full-master spectral history requires a complete HDA-certified master plus a matching finite spectral closure. "
            "Repository-level PHYSICAL_PROJECTOR_HISTORY remains open until an independent refinement/rigging-map certificate also passes. "
            "Connected correlator claims additionally require a matching source-dressed history certificate."
        ),
    }


def self_test() -> dict:
    identity = {
        "habitat_hash": "habitat:test",
        "domain_hash": "domain:test",
        "convention_hash": "convention:test",
        "master_pencil_hash": "master:test",
    }
    master = {
        "passed": True,
        "status": "COMPLETE_FINITE_FULL_DIRAC_MASTER_QUANTUM_HDA_CERTIFIED",
        "domain_complete": True,
        "quantum_hda_closed": True,
        "physical_projector_emitted": True,
        "quantum_hda_certificate_audit": {"valid_for_this_master": True},
        **identity,
        "spectrum": {"nullity": 2},
    }
    spectral = {
        "schema": SPECTRAL_SCHEMA,
        "status": "FINITE_SPECTRAL_HISTORY_CLOSED_PHYSICAL_CERTIFICATE_REQUIRED",
        "seed_label": "synthetic-evidence-link-control",
        "finite_spectral_history_closed": True,
        "physical_history_closed": False,
        "physical_projector_emitted": False,
        "provenance": dict(identity),
        "termination": {"mode": "direct_block_residual", "residual_norm": 0.0, "upstream_certified": True},
    }
    refinement = {
        "schema": REFINEMENT_SCHEMA,
        "passed": True,
        "refinement_family_hash": "refinement-family:test",
        "anchor_habitat_hash": identity["habitat_hash"],
        "anchor_domain_hash": identity["domain_hash"],
        "anchor_convention_hash": identity["convention_hash"],
        "anchor_master_pencil_hash": identity["master_pencil_hash"],
        "low_cluster_scale_separation": True,
        "projector_converged_under_embeddings": True,
        "boundary_history_converged": True,
        "hda_residual_converged": True,
        "rank_rule_preregistered_or_heldout": True,
    }
    source = {
        "schema": SOURCE_SCHEMA,
        "passed": True,
        **identity,
        "source_operator_set_complete": True,
        "source_dressed_history_converged": True,
        "connected_W_not_Z_used": True,
    }

    pos = certify(spectral, master, refinement, source)
    if not pos["passed"] or not pos["physical_projector_history_closed"] or not pos["connected_source_history_closed"]:
        raise AssertionError(pos)

    no_ref = certify(spectral, master)
    if not no_ref["finite_stage_passed"] or not no_ref["finite_full_master_spectral_history_certified"]:
        raise AssertionError("finite linked certificate unexpectedly failed")
    if no_ref["passed"] or no_ref["physical_projector_history_closed"]:
        raise AssertionError("missing refinement incorrectly closed physical history")

    bad_spec = json.loads(json.dumps(spectral))
    bad_spec["provenance"]["master_pencil_hash"] = "wrong-master"
    mismatch = certify(bad_spec, master, refinement, source)
    if mismatch["finite_stage_passed"] or mismatch["finite_full_master_spectral_history_certified"]:
        raise AssertionError("master-hash mismatch did not fail closed")

    bad_master = json.loads(json.dumps(master))
    bad_master["quantum_hda_certificate_audit"]["valid_for_this_master"] = False
    bad_hda = certify(spectral, bad_master, refinement, source)
    if bad_hda["finite_stage_passed"] or bad_hda["finite_full_master_spectral_history_certified"]:
        raise AssertionError("invalid HDA-master link did not fail closed")

    forged = json.loads(json.dumps(spectral))
    forged["physical_history_closed"] = True
    forged["physical_projector_emitted"] = True
    forged_result = certify(forged, master, refinement, source)
    if forged_result["finite_stage_passed"] or forged_result["finite_full_master_spectral_history_certified"]:
        raise AssertionError("self-promoted spectral result was accepted")

    return {
        "passed": True,
        "positive_status": pos["status"],
        "missing_refinement_status": no_ref["status"],
        "missing_refinement_top_level_passed": no_ref["passed"],
        "hash_mismatch_status": mismatch["status"],
        "hda_mismatch_status": bad_hda["status"],
        "self_promotion_rejected": True,
    }


def load_json(path: Path | None):
    if path is None:
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spectral-result", type=Path)
    ap.add_argument("--master-result", type=Path)
    ap.add_argument("--refinement-result", type=Path)
    ap.add_argument("--source-result", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        out = self_test()
    else:
        if args.spectral_result is None or args.master_result is None:
            ap.error("--spectral-result and --master-result are required unless --self-test")
        out = certify(
            load_json(args.spectral_result),
            load_json(args.master_result),
            load_json(args.refinement_result),
            load_json(args.source_result),
        )

    text = json.dumps(out, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0 if out.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
