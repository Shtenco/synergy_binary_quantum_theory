#!/usr/bin/env python3
"""Audit the exact trivial [5] selector path through the recovered vacuum row.

The pinned source CFG omits [5].  This non-mutating wrapper adds its analytically
derived Jucys data in memory, constructs the one-dimensional vacuum selector,
and calls the recovered generic master row with that selector.  It does not
certify any non-vacuum selector multiplicity or rank.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "depth6_recovered_source"
sys.path.insert(0, str(SOURCE))

import bqg_depth6_generic_jucys_selector as J  # noqa: E402
import bqg_depth6_generic_master_row as M  # noqa: E402

EXPECTED_BUNDLE = "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32"
TRIVIAL_S5_CFG = ("triv", 3, (-1, 3), 4, (-1, 4))


def verify_source_hashes() -> dict[str, str]:
    manifest = json.loads((SOURCE / "source_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("bundle_sha256") != EXPECTED_BUNDLE:
        raise RuntimeError("recovered source bundle hash mismatch")
    hashes = {
        name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
        for name in manifest["module_sha256"]
    }
    if hashes != manifest["module_sha256"]:
        raise RuntimeError("recovered source module hashes mismatch")
    return hashes


def run_probe() -> dict[str, object]:
    hashes = verify_source_hashes()
    if "5" in J.CFG:
        raise RuntimeError("pinned source unexpectedly already defines the [5] CFG")
    # J4 is the S4 transposition sum: eigenvalues -1,3; J5 has -1,4.
    # The selected targets 3 and 4 are the trivial branches at both steps.
    J.CFG["5"] = TRIVIAL_S5_CFG
    vacuum = (0,) * 10
    seed = J.CFG["5"][0]
    _blocks, _slices, selector_dim = J.layout(vacuum, seed)
    if selector_dim != 1:
        raise RuntimeError(("vacuum selector dimension", selector_dim))
    selector = J.apply_selector_projector(vacuum, "5", np.eye(selector_dim, dtype=complex))
    selector_error = float(np.max(np.abs(selector - np.eye(selector_dim))))
    if selector_error != 0.0:
        raise ArithmeticError(("vacuum selector is not identity", selector_error, selector))

    output, _ = M.master_map(vacuum, "5", W=selector)
    output_norm = float(np.sqrt(sum(np.linalg.norm(matrix) ** 2 for matrix in output.values())))
    if output or output_norm != 0.0:
        raise AssertionError(("generic [5] vacuum row is nonzero", len(output), output_norm))

    return {
        "status": "EXACT_TRIVIAL_SELECTOR_VACUUM_ROW_ZERO",
        "source": {
            "commit": "f6bfcd775b8c3cd33c6c235d48293e40657baf9d",
            "bundle_sha256": EXPECTED_BUNDLE,
            "module_sha256": hashes,
            "source_cfg_unchanged_on_disk": True,
        },
        "selector": {
            "temporary_cfg": list(TRIVIAL_S5_CFG),
            "derived_projector": "((J4+I)/4)((J5+I)/5)e3=e5",
            "dimension": selector_dim,
            "max_deviation_from_identity": selector_error,
            "proof_basis": "exact coset-average identity verified by scripts/verify_s5_trivial_projector.py",
        },
        "generic_master_row": {
            "input": "zero-spin vacuum coordinate",
            "output_support": len(output),
            "output_norm": output_norm,
            "evaluated_by": "recovered bqg_depth6_generic_master_row.master_map with explicit one-dimensional selector W",
        },
        "claim_limit": (
            "This pins the [5] trivial-selector and generic-row vacuum path only. The recovered source CFG still omits [5] on disk; "
            "the generic row runtime is floating-point; the historical B_[5] operator identity and non-vacuum rank remain open."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(run_probe(), indent=2, ensure_ascii=False))
