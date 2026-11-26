#!/usr/bin/env python3
"""Fail-closed verifier for the corrected finite depth-6 BQG frontier.

This validator deliberately distinguishes:
- finite numerical master closure,
- structural branch-sum closure,
- pending independent numerical witnesses,
- the full finite depth-6 theorem,
- continuum/physical completion.

It MUST fail if the canonical ledger reopens already closed structural work,
reverts [2,1,1,1] to ACTIVE, or promotes the full theorem beyond preserved
evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "depth6_frontier.json"
CERT_2111 = ROOT / "BQG_DEPTH6_2111_MASTER_CLOSED_2026-10-01.json"
C5_CENSUS = ROOT / "verification_results" / "DEPTH6_5_ORBIT_BLOCK_CENSUS_20261126.json"
C5_SOURCE_BRIDGE = ROOT / "verification_results" / "DEPTH6_5_VACUUM_SOURCE_BRIDGE_20261126.json"
C5_EXACT_BRIDGE = ROOT / "verification_results" / "DEPTH6_5_VACUUM_COLLAPSED_EXACT_EQUIVALENCE_20261126.json"
C5_SELECTOR_PROBE = ROOT / "verification_results" / "DEPTH6_5_SELECTOR_VACUUM_GENERIC_ROW_20261126.json"

EXPECTED_MULT = {
    "[5]": 27227,
    "[4,1]": 104146,
    "[3,2]": 130903,
    "[3,1,1]": 153455,
    "[2,2,1]": 130503,
    "[2,1,1,1]": 103318,
    "[1^5]": 26794,
}
S5_DIMS = {
    "[5]": 1,
    "[4,1]": 4,
    "[3,2]": 5,
    "[3,1,1]": 6,
    "[2,2,1]": 5,
    "[2,1,1,1]": 4,
    "[1^5]": 1,
}
EXPECTED_STRUCTURAL = {
    "[3,2]": ("2755/2755", "130903/130903", 0),
    "[3,1,1]": ("2719/2719", "153455/153455", 0),
    "[2,2,1]": ("2749/2749", "130503/130503", 0),
    "[2,1,1,1]": ("2712/2712", "103318/103318", 0),
}


def require(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    errors: list[str] = []
    data = json.loads(LEDGER.read_text(encoding="utf-8"))
    cert = json.loads(CERT_2111.read_text(encoding="utf-8"))
    census = json.loads(C5_CENSUS.read_text(encoding="utf-8"))
    source_bridge = json.loads(C5_SOURCE_BRIDGE.read_text(encoding="utf-8"))
    exact_bridge = json.loads(C5_EXACT_BRIDGE.read_text(encoding="utf-8"))
    selector_probe = json.loads(C5_SELECTOR_PROBE.read_text(encoding="utf-8"))

    require(errors, data.get("schema_version") == 2, "ledger schema_version must be 2")
    require(errors, cert.get("schema_version") == 1, "2111 certificate schema_version must be 1")
    require(errors, data.get("status") == "NUMERICAL_WITNESSES_PENDING", "frontier status mismatch")

    shell = data.get("shell", {})
    require(errors, shell.get("gauss_admissible_spin_assignments") == 264962, "Gauss shell count mismatch")
    require(errors, shell.get("hilbert_dimension") == 3111637, "Hilbert dimension mismatch")
    require(errors, shell.get("s5_spin_orbits") == 2757, "S5 orbit count mismatch")

    mult = data.get("multiplicities", {})
    require(errors, mult == EXPECTED_MULT, "S5 multiplicity table mismatch")
    reconstructed = sum(EXPECTED_MULT[k] * S5_DIMS[k] for k in EXPECTED_MULT)
    require(errors, reconstructed == 3111637, f"internal dimension reconstruction mismatch: {reconstructed}")

    # Structural support/capacity is already closed for all seven S5 irreps.
    require(
        errors,
        data.get("structural_depth6_status") == "CLOSED_FOR_ALL_7_S5_IRREPS",
        "all-seven-irrep structural closure must remain frozen",
    )
    support = data.get("structural_support", {})
    for key, (blocks, columns, remainder) in EXPECTED_STRUCTURAL.items():
        s = support.get(key, {})
        require(errors, s.get("status") == "CLOSED", f"{key} structural status must remain CLOSED")
        require(errors, s.get("blocks") == blocks, f"{key} structural block ledger mismatch")
        require(errors, s.get("columns") == columns, f"{key} structural column ledger mismatch")
        require(errors, s.get("remainder") == remainder, f"{key} structural remainder must remain zero")

    sector_evidence = data.get("sector_kernel_evidence", {})
    required_sectors = {"[1^5]", "[5]", "[4,1]", "[2,1,1,1]"}
    require(
        errors,
        set(sector_evidence) == required_sectors,
        "sector_kernel_evidence must describe the four previously claimed sectors",
    )
    # [1^5] is a finite numerical consequence of the saved whole S4-sign result,
    # not an exact-arithmetic proof and not an independent sector certificate.
    c1f5 = sector_evidence.get("[1^5]", {})
    require(errors, c1f5.get("status") == "CLOSED_FINITE_NUMERICAL", "[1^5] must remain finite numerical only")
    require(errors, c1f5.get("rank") == 26794, "[1^5] numerical rank mismatch")
    require(errors, c1f5.get("dimension") == 26794, "[1^5] dimension mismatch")
    require(errors, c1f5.get("kernel") == "empty", "[1^5] numerical kernel mismatch")
    require(errors, c1f5.get("certificate") == CERT_2111.name, "[1^5] S4-sign certificate pointer mismatch")
    require(errors, "whole S4-sign" in c1f5.get("derivation", ""), "[1^5] numerical derivation must identify its scope")

    # [5] now has an exact vacuum-inclusion witness for the recovered node map,
    # but still lacks both a nonvacuum rank certificate and identity with the
    # historical branch-sum operator. Preserve the original claim as unverified.
    c5 = sector_evidence.get("[5]", {})
    require(errors, c5.get("status") == "UNVERIFIED_LEDGER_CLAIM", "[5] unsupported CLOSED claim must remain unverified")
    require(errors, c5.get("ledger_claim", {}).get("dimension") == 27227, "[5] original dimension claim missing")
    require(errors, c5.get("ledger_claim", {}).get("nonvacuum_rank") == 27226, "[5] original rank claim missing")
    require(errors, c5.get("evidence") is None, "[5] must not imply a missing certificate exists")
    census_ref = c5.get("orbit_block_census", {})
    require(errors, census_ref.get("status") == "STRUCTURAL_BLOCK_CENSUS_ONLY", "[5] orbit census scope missing")
    require(errors, census_ref.get("artifact") == "verification_results/DEPTH6_5_ORBIT_BLOCK_CENSUS_20261126.json", "[5] orbit census artifact pointer mismatch")
    require(errors, census_ref.get("active_s5_orbit_blocks") == 2757, "[5] orbit census block count mismatch")
    require(errors, census_ref.get("summed_trivial_selector_multiplicity") == 27227, "[5] orbit census multiplicity mismatch")
    require(errors, census_ref.get("vacuum_block_multiplicity") == 1, "[5] vacuum block dimension mismatch")
    require(errors, census_ref.get("nonvacuum_active_s5_orbit_blocks") == 2756, "[5] nonvacuum block count mismatch")
    require(errors, census_ref.get("nonvacuum_selector_multiplicity") == 27226, "[5] nonvacuum multiplicity mismatch")
    require(errors, census_ref.get("rank_claim") is None, "[5] orbit census must not claim rank")
    require(errors, census.get("status") == "STRUCTURAL_BLOCK_CENSUS_ONLY", "[5] census artifact missing or overstated")
    require(errors, census.get("active_trivial_s5_blocks") == 2757, "[5] census artifact block count mismatch")
    require(errors, census.get("summed_trivial_selector_multiplicity") == 27227, "[5] census artifact multiplicity mismatch")
    require(errors, census.get("nonvacuum_active_blocks") == 2756, "[5] nonvacuum census block count mismatch")
    require(errors, census.get("nonvacuum_selector_multiplicity") == 27226, "[5] nonvacuum census dimension mismatch")
    require(errors, census.get("rank_claim") is None, "[5] census artifact must not claim rank")
    selector = c5.get("trivial_sector_selector", {})
    require(errors, selector.get("status") == "EXACT_ALGEBRAIC_GROUP_ALGEBRA_IDENTITY", "[5] exact sector selector identity missing")
    require(errors, selector.get("projector") == "((I+J5)/5)((I+J4)/4)e3 = e5", "[5] exact coset-projector formula mismatch")
    require(errors, selector.get("verifier") == "scripts/verify_s5_trivial_projector.py", "[5] exact selector verifier pointer mismatch")
    require(errors, "no master-map rank" in selector.get("scope_limit", ""), "[5] selector must not overstate master-map rank")
    vacuum = c5.get("reconstructed_operator_vacuum_inclusion", {})
    require(errors, vacuum.get("status") == "EXACT_FINITE_FOR_RECONSTRUCTED_NODE_FORMULA", "[5] reconstructed exact vacuum result missing or overstated")
    require(errors, vacuum.get("node_constraint") == "H_v |0> = 0 for each v in {0,1,2,3,4} in the exact reduced-basis reconstruction", "[5] exact node-vacuum scope mismatch")
    require(errors, vacuum.get("master_consequence") == "<0|sum_v H_v^dagger H_v|0> = 0 for the positive master assembled from the reconstructed node maps", "[5] exact vacuum master consequence mismatch")
    require(errors, len(vacuum.get("artifacts", [])) == 5, "[5] exact vacuum artifact list incomplete")
    require(errors, "source Python execution itself remains floating-point" in vacuum.get("limitation", ""), "[5] source runtime precision limitation missing")
    require(errors, "source selector CFG has no [5] entry" in vacuum.get("limitation", ""), "[5] source selector limitation missing")
    require(errors, "historical branch-sum identities remain unproved" in vacuum.get("limitation", ""), "[5] historical-operator limitation missing")
    require(errors, "Exact finite equality on the vacuum column" in vacuum.get("node_map_identity_to_source", ""), "[5] exact vacuum source node-map bridge missing")
    require(errors, "Vacuum-only inclusion is now established" in vacuum.get("master_map_assembly", ""), "[5] recovered generic-row vacuum inclusion missing")
    selector_ref = vacuum.get("selector_vacuum_probe", {})
    require(errors, selector_ref.get("status") == "EXACT_TRIVIAL_SELECTOR_VACUUM_ROW_ZERO", "[5] selector vacuum probe status mismatch")
    require(errors, selector_ref.get("script") == "scripts/audit_depth6_5_selector_vacuum.py", "[5] selector vacuum script pointer mismatch")
    require(errors, selector_ref.get("artifact") == "verification_results/DEPTH6_5_SELECTOR_VACUUM_GENERIC_ROW_20261126.json", "[5] selector vacuum artifact pointer mismatch")
    require(errors, selector_ref.get("selector_dimension") == 1 and selector_ref.get("generic_master_row_output_support") == 0, "[5] selector vacuum dimensions/output mismatch")
    require(errors, selector_probe.get("status") == "EXACT_TRIVIAL_SELECTOR_VACUUM_ROW_ZERO", "[5] selector vacuum artifact missing or overstated")
    require(errors, selector_probe.get("selector", {}).get("temporary_cfg") == ["triv", 3, [-1, 3], 4, [-1, 4]], "[5] trivial selector Jucys config mismatch")
    require(errors, selector_probe.get("selector", {}).get("dimension") == 1 and selector_probe.get("selector", {}).get("max_deviation_from_identity") == 0.0, "[5] vacuum selector identity mismatch")
    require(errors, selector_probe.get("generic_master_row", {}).get("output_support") == 0 and selector_probe.get("generic_master_row", {}).get("output_norm") == 0.0, "[5] generic master-row vacuum output mismatch")
    require(errors, selector_probe.get("source", {}).get("bundle_sha256") == "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32", "[5] selector probe source bundle hash mismatch")
    require(errors, len(selector_probe.get("source", {}).get("module_sha256", {})) == 13, "[5] selector probe module hash list incomplete")
    require(errors, selector_probe.get("source", {}).get("module_sha256") == census.get("source", {}).get("module_sha256"), "[5] selector probe source module hashes differ from census")
    require(errors, "historical B_[5] operator identity" in selector_probe.get("claim_limit", ""), "[5] selector probe historical-identity limitation missing")
    require(errors, "non-vacuum rank remain open" in selector_probe.get("claim_limit", ""), "[5] selector probe rank limitation missing")
    bridge_ref = vacuum.get("source_bridge_crosscheck", {})
    require(errors, bridge_ref.get("status") == "NUMERICAL_SOURCE_BRIDGE_CROSSCHECK", "[5] source bridge scope mismatch")
    require(errors, bridge_ref.get("artifact") == "verification_results/DEPTH6_5_VACUUM_SOURCE_BRIDGE_20261126.json", "[5] source bridge artifact pointer mismatch")
    require(errors, bridge_ref.get("state_resolved_channel_matrices") == 60, "[5] source bridge matrix coverage mismatch")
    require(errors, bridge_ref.get("channel_entries_compared") == 320, "[5] source bridge entry coverage mismatch")
    require(errors, source_bridge.get("status") == "NUMERICAL_SOURCE_BRIDGE_CROSSCHECK", "[5] numerical source bridge artifact missing or overstated")
    require(errors, source_bridge.get("state_resolved_channel_matrices") == 60, "[5] numerical source bridge matrix count mismatch")
    require(errors, source_bridge.get("channel_entries_compared") == 320, "[5] numerical source bridge entry count mismatch")
    require(errors, source_bridge.get("max_absolute_error", 1.0) <= 7e-16, "[5] numerical source bridge absolute error exceeds recorded bound")
    require(errors, source_bridge.get("max_relative_error_on_nonzero_exact_entries", 1.0) <= 2e-15, "[5] numerical source bridge relative error exceeds recorded bound")
    require(errors, source_bridge.get("max_source_leakage_at_exact_zero_entries", 1.0) <= 8e-17, "[5] numerical source zero-entry leakage exceeds recorded bound")
    require(errors, "not an exact-arithmetic proof" in source_bridge.get("claim_limit", ""), "[5] numerical source bridge limitation missing")
    exact_ref = vacuum.get("exact_collapsed_equivalence", {})
    require(errors, exact_ref.get("status") == "EXACT_FINITE_COLLAPSED_EQUIVALENCE_ON_VACUUM", "[5] exact collapsed-stage status mismatch")
    require(errors, exact_ref.get("script") == "scripts/verify_depth6_vacuum_collapsed_exact_equivalence.py", "[5] exact collapsed-stage script pointer mismatch")
    require(errors, exact_ref.get("artifact") == "verification_results/DEPTH6_5_VACUUM_COLLAPSED_EXACT_EQUIVALENCE_20261126.json", "[5] exact collapsed-stage artifact pointer mismatch")
    require(errors, exact_ref.get("channel_state_entries_compared") == 240, "[5] exact collapsed-stage entry count mismatch")
    require(errors, exact_bridge.get("status") == "EXACT_FINITE_COLLAPSED_EQUIVALENCE_ON_VACUUM", "[5] exact collapsed-stage artifact missing or overstated")
    require(errors, exact_bridge.get("source", {}).get("bundle_sha256") == "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32", "[5] exact collapsed-stage source bundle hash mismatch")
    require(errors, exact_bridge.get("source", {}).get("module_sha256") == census.get("source", {}).get("module_sha256"), "[5] exact collapsed-stage source module hashes differ from census")
    require(errors, exact_bridge.get("exact_identity_checks", {}).get("state_channel_entries") == 240, "[5] exact collapsed-stage artifact entry count mismatch")
    require(errors, len(exact_bridge.get("exact_identity_checks", {}).get("per_node", [])) == 5, "[5] exact collapsed-stage node coverage mismatch")
    require(errors, all(row.get("oriented_terms") == 12 and row.get("triangle_exact_state_keys_compared") == 48 and row.get("radial_exact_state_keys_compared") == 48 and row.get("channel_state_entries_compared") == 48 for row in exact_bridge.get("exact_identity_checks", {}).get("per_node", [])), "[5] exact collapsed-stage per-node coverage mismatch")
    exact_checks = exact_bridge.get("exact_identity_checks", {})
    require(errors, "symbolically equal" in exact_checks.get("expanded_vs_collapsed_triangle_adjoint", ""), "[5] exact triangle identity claim missing")
    require(errors, "symbolically equal" in exact_checks.get("expanded_vs_collapsed_radial_adjoint", ""), "[5] exact radial identity claim missing")
    exact_limit = exact_bridge.get("claim_limit", "")
    require(errors, "historical B_[5] identity remain open" in exact_limit, "[5] exact collapsed-stage historical identity limitation missing")
    require(errors, "no non-vacuum rank is implied" in exact_limit, "[5] exact collapsed-stage rank limitation missing")
    for node in range(5):
        artifact_path = ROOT / "verification_results" / f"DEPTH6_5_VACUUM_ADJOINT_EXACT_NODE{node}_20261126.json"
        try:
            artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"[5] exact vacuum node {node} artifact unavailable/invalid: {exc}")
            continue
        require(errors, artifact.get("node") == node, f"[5] exact vacuum node {node} label mismatch")
        require(errors, artifact.get("exact_cg_arithmetic") is True, f"[5] node {node} is not exact-CG arithmetic")
        require(errors, artifact.get("node_constraint_annihilates_vacuum") is True, f"[5] node {node} vacuum image is nonzero")
        require(errors, artifact.get("node_constraint_pre_gauss_output") == {}, f"[5] node {node} pre-Gauss exact output is nonzero")
        require(errors, artifact.get("source_gauss_projection_applied") is True, f"[5] node {node} source Gauss projection missing")
        terms = artifact.get("terms", [])
        require(errors, len(terms) == 12, f"[5] node {node} oriented-term coverage mismatch")
        require(errors, all(term.get("term_annihilates_vacuum") is True for term in terms), f"[5] node {node} has a nonzero term output")
        allowed_channel_forms = {
            "Matrix([[3**(1/4)/6,3**(1/4)/3],[-3**(1/4)/3,-3**(1/4)/6]])",
            "Matrix([[-3**(1/4)/6,-3**(1/4)/3],[3**(1/4)/3,3**(1/4)/6]])",
        }
        channel_matrices = [
            "".join(matrix.split())
            for term in terms
            for matrix in term.get("adjoint_channel_matrices_by_output_state", {}).values()
        ]
        require(errors, len(channel_matrices) == 12, f"[5] node {node} state-resolved channel matrix count mismatch")
        require(errors, all(matrix in allowed_channel_forms for matrix in channel_matrices), f"[5] node {node} exact channel matrix form mismatch")
        provenance = artifact.get("source_provenance", {})
        require(errors, provenance.get("commit") == "f6bfcd775b8c3cd33c6c235d48293e40657baf9d", f"[5] node {node} source commit mismatch")
        require(errors, provenance.get("bundle_sha256") == "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32", f"[5] node {node} bundle hash mismatch")
        require(errors, "final local-J=0 Gauss projection" in provenance.get("comparison_scope", ""), f"[5] node {node} source Gauss projection provenance missing")
        require(errors, provenance.get("identity_with_historical_branch_sum") == "not established", f"[5] node {node} historical identity limitation missing")
        require(errors, len(provenance.get("source_module_sha256", {})) == 5, f"[5] node {node} source module hash list missing")
    c41 = sector_evidence.get("[4,1]", {})
    require(errors, c41.get("status") == "UNVERIFIED_LEDGER_CLAIM", "[4,1] unsupported CLOSED claim must remain unverified")
    require(errors, c41.get("ledger_claim", {}).get("rank") == 104146, "[4,1] original rank claim missing")
    require(errors, c41.get("evidence") is None, "[4,1] must not imply a missing certificate exists")

    c2111 = sector_evidence.get("[2,1,1,1]", {})
    require(errors, c2111.get("status") == "CLOSED_FINITE_NUMERICAL", "[2,1,1,1] must be CLOSED_FINITE_NUMERICAL")
    require(errors, c2111.get("rank") == 103318, "[2,1,1,1] rank mismatch")
    require(errors, c2111.get("dimension") == 103318, "[2,1,1,1] dimension mismatch")
    require(errors, c2111.get("kernel") == "empty", "[2,1,1,1] kernel must be empty")
    require(errors, c2111.get("certificate") == CERT_2111.name, "[2,1,1,1] certificate pointer mismatch")

    # Verify the preserved [2,1,1,1] master certificate rather than trusting only the ledger.
    require(errors, cert.get("status") == "PASS", "2111 master certificate must be PASS")
    require(errors, cert.get("input_dimension") == 130112, "2111 S4-sign input dimension mismatch")
    h0 = cert.get("h0", {})
    zero = h0.get("exact_zero_subspace", {})
    require(errors, zero.get("dimension") == 16, "2111 exact H0 kernel dimension mismatch")
    require(errors, h0.get("complement_dimension") == 130096, "2111 H0 complement dimension mismatch")
    main_q = h0.get("main_unique_q_certificate", {})
    residual = h0.get("residual_sparse_certificate", {})
    require(errors, main_q.get("columns") == 123744, "2111 main H0 column count mismatch")
    require(errors, residual.get("columns") == 6352, "2111 residual H0 column count mismatch")
    require(errors, main_q.get("columns", 0) + residual.get("columns", 0) == 130096, "2111 H0 complement partition mismatch")
    require(errors, float(main_q.get("selected_sigma_min", 0.0)) > 1e-2, "2111 main H0 sigma_min below certified threshold")
    require(errors, float(residual.get("sigma_min", 0.0)) > 0.3, "2111 residual H0 sigma_min unexpectedly weak")
    require(errors, bool(h0.get("decomposition_row_disjoint")), "2111 H0 decomposition must remain row-disjoint")

    h1 = cert.get("h1_on_exact_h0_kernel", {})
    require(errors, h1.get("dimension") == 16 and h1.get("rank") == 16, "2111 H1 lift must remain rank 16/16")
    require(errors, float(h1.get("sigma_min", 0.0)) > 1.0, "2111 H1 sigma_min unexpectedly weak")
    require(errors, bool(h1.get("fresh_recompute")), "2111 H1 lift must come from fresh recompute")

    master = cert.get("master", {})
    require(errors, master.get("s4_sign_kernel_dimension") == 0, "2111 S4-sign master kernel must be empty")
    require(errors, master.get("s4_sign_status") == "CLOSED", "2111 S4-sign certificate not CLOSED")
    consequence = cert.get("irrep_consequence", {}).get("[2,1,1,1]", {})
    require(errors, consequence.get("dimension") == 103318, "2111 irrep consequence dimension mismatch")
    require(errors, consequence.get("rank") == 103318, "2111 irrep consequence rank mismatch")
    require(errors, consequence.get("kernel") == "empty", "2111 irrep consequence kernel mismatch")
    require(errors, consequence.get("status") == "CLOSED", "2111 irrep consequence status mismatch")

    s4 = data.get("s4_sign_frontier", {})
    require(errors, s4.get("total_dimension") == 130112, "S4-sign total dimension mismatch")
    require(errors, s4.get("exact_h0_null_directions") == 16, "S4-sign exact H0-null dimension mismatch")
    require(errors, s4.get("h0_complement_dimension") == 130096, "S4-sign H0 complement dimension mismatch")
    require(errors, s4.get("h0_complement_rank") == 130096, "S4-sign H0 complement rank mismatch")
    require(errors, s4.get("master_kernel_dimension") == 0, "S4-sign master kernel dimension mismatch")
    require(errors, s4.get("sector_status") == "CLOSED_FINITE_NUMERICAL", "[2,1,1,1] frontier must not regress to ACTIVE")
    require(errors, s4.get("obsolete_gate", {}).get("status") == "SUPERSEDED_DO_NOT_USE_AS_BLOCKER", "obsolete 130007 minor must not remain a blocking gate")

    # [3,2] is structurally closed; only its independent numerical witness remains.
    witnesses = data.get("numerical_master_witnesses", {})
    w32 = witnesses.get("[3,2]", {})
    require(errors, w32.get("structural_status") == "CLOSED", "[3,2] structural proof was incorrectly reopened")
    require(errors, w32.get("dimension") == 130903, "[3,2] numerical witness dimension mismatch")
    require(errors, w32.get("blocks") == 2755, "[3,2] numerical witness block count mismatch")
    require(errors, w32.get("status") == "RECOVER_OR_FINISH", "[3,2] must be recovery-first, not structural NEXT")
    require(errors, data.get("current_numerical_witness") == "[3,2]", "current numerical witness target mismatch")
    require(
        errors,
        data.get("remaining_numerical_irreps") == ["[3,2]", "[3,1,1]", "[2,2,1]"],
        "remaining numerical irrep frontier mismatch",
    )

    do_not = data.get("do_not_recompute_for_32", [])
    require(errors, any("130903/130903" in x for x in do_not), "[3,2] structural 130903/130903 must be frozen as do-not-recompute")
    require(errors, any("Jucys" in x for x in do_not), "[3,2] Jucys selector work must be frozen as do-not-recompute")

    require(errors, data.get("finite_depth6_theorem_status") == "NOT_YET_PROVED", "full finite depth-6 theorem must remain NOT_YET_PROVED")
    require(errors, "5 A_[5],[4]" in data.get("branch_sum", {}).get("[5]", ""), "[5] unique-branch conditional formula missing")
    require(errors, "operator identity with recovered generic row remains unproved" in data.get("branch_sum", {}).get("[5]", ""), "[5] historical operator identity limitation missing")
    hist5 = data.get("historical_5_claim_audit", {})
    require(errors, hist5.get("status") == "HISTORICAL_LABEL_FOUND_NO_RANK_WITNESS_VALIDATED", "[5] historical audit status mismatch")
    require(errors, hist5.get("source_branch") == "archive/w5-v09-immutable-6ae145a", "[5] historical audit branch mismatch")
    require(errors, hist5.get("source_frontier_blob_sha") == "e9ca9dc8ffa69df908f7ae1348701a5f738e1202", "[5] archived frontier blob pin mismatch")
    require(errors, hist5.get("source_verifier_blob_sha") == "7f9cf62b8386bb464ecdf62efac9d28a58a16375", "[5] archived verifier blob pin mismatch")
    archive_audit = hist5.get("retained_branch_archive_corpus", {})
    require(errors, archive_audit.get("archive_count") == 179, "[5] archive corpus size mismatch")
    require(errors, archive_audit.get("copies_containing_this_claim") == 42, "[5] repeated historical claim copy count mismatch")
    require(errors, archive_audit.get("distinct_frontier_content_count") == 1, "[5] historical frontier copies are not one identical snapshot")
    require(errors, archive_audit.get("frontier_content_sha256") == "ed33ad07523803ec3b979cbae1f1755f5b3d67e12491dc127b4d9eb29e81286c", "[5] historical frontier content hash mismatch")
    require(errors, "not 42 independent rank checks" in archive_audit.get("interpretation", ""), "[5] repeated copies must not be counted as independent evidence")
    old5 = hist5.get("archived_claim", {})
    require(errors, old5.get("global_status") == "NUMERICAL_WITNESSES_PENDING", "[5] historical global pending status mismatch")
    require(errors, old5.get("sector_status") == "CLOSED" and old5.get("dimension") == 27227, "[5] historical closed sector claim mismatch")
    require(errors, old5.get("nonvacuum_rank") == 27226 and old5.get("certificate_pointer") is None, "[5] historical rank/certificate provenance mismatch")
    require(errors, "[5] nonvacuum rank" in hist5.get("verifier_does_not_check", []), "[5] historical verifier rank gap missing")
    require(errors, "[5] certificate or numerical witness" in hist5.get("verifier_does_not_check", []), "[5] historical verifier witness gap missing")
    require(errors, "not a disproof" in hist5.get("interpretation", ""), "[5] historical audit must not overstate the result")
    require(errors, cert.get("finite_depth6_full_theorem_status") == "NOT_YET_PROVED", "2111 certificate must not promote full depth-6 theorem")

    result = {
        "valid": not errors,
        "status": data.get("status"),
        "structural_depth6_status": data.get("structural_depth6_status"),
        "finite_depth6_theorem_status": data.get("finite_depth6_theorem_status"),
        "sector_kernel_evidence": {
            key: value.get("status") for key, value in sorted(sector_evidence.items())
        },
        "s4_sign_status": s4.get("sector_status"),
        "s4_sign_master_kernel_dimension": s4.get("master_kernel_dimension"),
        "current_numerical_witness": data.get("current_numerical_witness"),
        "32_structural_status": w32.get("structural_status"),
        "32_numerical_status": w32.get("status"),
        "errors": errors,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
