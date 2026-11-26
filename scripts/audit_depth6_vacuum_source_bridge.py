#!/usr/bin/env python3
"""Numerically cross-check the exact vacuum evaluator against recovered source stages.

This is a convention/provenance audit, not an exact-arithmetic source identity:
the recovered source stage still uses floating CG/recoupling coefficients.
"""
from __future__ import annotations

import ast
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "depth6_recovered_source"
sys.path.insert(0, str(SOURCE))

import bqg_block_apply_invariants as A  # noqa: E402

NODE_ARTIFACTS = ROOT / "verification_results"
EXPECTED_BUNDLE_SHA256 = "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32"


def verify_source() -> dict[str, str]:
    manifest = json.loads((SOURCE / "source_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("bundle_sha256") != EXPECTED_BUNDLE_SHA256:
        raise RuntimeError("recovered source bundle digest mismatch")
    actual = {
        name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
        for name in manifest["module_sha256"]
    }
    if actual != manifest["module_sha256"]:
        raise RuntimeError("recovered source module digest mismatch")
    return actual


def run() -> dict[str, object]:
    hashes = verify_source()
    vacuum_spins = (0,) * 10
    initial = A.initial_reduced_block(vacuum_spins, np.eye(1, dtype=complex))
    node_summaries = []
    global_max_abs = 0.0
    global_max_rel_nonzero = 0.0
    global_max_zero_leakage = 0.0
    state_channel_count = 0
    compared_channels = 0

    for node in range(5):
        artifact_path = NODE_ARTIFACTS / f"DEPTH6_5_VACUUM_ADJOINT_EXACT_NODE{node}_20261126.json"
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        if artifact.get("node") != node or artifact.get("exact_cg_arithmetic") is not True:
            raise RuntimeError(f"invalid exact-CG node artifact: {artifact_path.name}")
        specs = A.R.oriented_specs(node)
        if len(specs) != 12 or len(artifact.get("terms", [])) != 12:
            raise RuntimeError(f"expected 12 oriented terms for node {node}")

        node_max_abs = 0.0
        node_max_zero_leakage = 0.0
        node_compared = 0
        node_state_channels = 0
        for (sign, spec), term in zip(specs, artifact["terms"], strict=True):
            v, a, b, c = spec
            if term.get("spec") != [v, a, b, c] or term.get("orientation_sign") != sign:
                raise RuntimeError(f"oriented-term ordering mismatch at node={node}, spec={spec}")
            exact_matrices = {
                ast.literal_eval(state_key): sp.sympify(matrix_text, locals={"Matrix": sp.Matrix})
                for state_key, matrix_text in term["adjoint_channel_matrices_by_output_state"].items()
            }
            node_state_channels += len(exact_matrices)

            for i in range(2):
                for j in range(2):
                    # The source's collapsed triangle and radial stages implement
                    # the same adjoint ordering as the exact expanded evaluator.
                    loop = A.triangle_adjoint_collapsed(initial, v, a, b, i, j, 1, tol=1e-30)
                    source_output = A.apply_radial_collapsed_block(
                        loop, v, c, i, j, adj=True, tol=1e-30
                    )
                    keys = set(exact_matrices) | set(source_output)
                    for key in keys:
                        exact_matrix = exact_matrices.get(key)
                        exact_value = complex(sp.N(exact_matrix[i, j], 40)) if exact_matrix is not None else 0j
                        source_value = complex(np.asarray(source_output.get(key, np.zeros(1)))[0])
                        abs_error = abs(source_value - exact_value)
                        node_max_abs = max(node_max_abs, abs_error)
                        global_max_abs = max(global_max_abs, abs_error)
                        if exact_value == 0:
                            node_max_zero_leakage = max(node_max_zero_leakage, abs_error)
                            global_max_zero_leakage = max(global_max_zero_leakage, abs_error)
                        else:
                            global_max_rel_nonzero = max(
                                global_max_rel_nonzero, abs_error / abs(exact_value)
                            )
                        node_compared += 1
                        compared_channels += 1

        node_summaries.append({
            "node": node,
            "oriented_terms": len(specs),
            "state_resolved_channel_matrices": node_state_channels,
            "source_vs_exact_channel_entries_compared": node_compared,
            "max_absolute_error": node_max_abs,
            "max_leakage_at_exact_zero_entries": node_max_zero_leakage,
        })
        state_channel_count += node_state_channels

    return {
        "status": "NUMERICAL_SOURCE_BRIDGE_CROSSCHECK",
        "source": {
            "commit": "f6bfcd775b8c3cd33c6c235d48293e40657baf9d",
            "bundle_sha256": EXPECTED_BUNDLE_SHA256,
            "module_sha256": hashes,
            "backend": "bqg_block_apply_invariants.triangle_adjoint_collapsed + apply_radial_collapsed_block",
        },
        "input": "all-zero doubled-spin vacuum, 5 nodes, 12 oriented terms/node, four SU(2) color pairs/term",
        "exact_reference": "DEPTH6_5_VACUUM_ADJOINT_EXACT_NODE{0..4}_20261126.json",
        "state_resolved_channel_matrices": state_channel_count,
        "channel_entries_compared": compared_channels,
        "max_absolute_error": global_max_abs,
        "max_relative_error_on_nonzero_exact_entries": global_max_rel_nonzero,
        "max_source_leakage_at_exact_zero_entries": global_max_zero_leakage,
        "node_summaries": node_summaries,
        "claim_limit": "Numerical agreement checks conventions and optimized source stages on this vacuum input. It is not an exact-arithmetic proof of source/reconstruction identity, does not add a [5] selector to CFG, and does not establish identity with historical B_[5].",
    }


if __name__ == "__main__":
    result = run()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "verification_results" / "DEPTH6_5_VACUUM_SOURCE_BRIDGE_20261126.json",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "channel_entries_compared": result["channel_entries_compared"],
        "state_resolved_channel_matrices": result["state_resolved_channel_matrices"],
        "max_absolute_error": result["max_absolute_error"],
        "max_relative_error_on_nonzero_exact_entries": result["max_relative_error_on_nonzero_exact_entries"],
        "max_source_leakage_at_exact_zero_entries": result["max_source_leakage_at_exact_zero_entries"],
        "claim_limit": result["claim_limit"],
        "artifact": str(args.output),
    }, indent=2))

