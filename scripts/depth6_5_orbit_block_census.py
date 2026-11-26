#!/usr/bin/env python3
"""Rebuild the recovered depth-6 [5] orbit-block census.

This is a structural/numerical census, not a master-map rank certificate.
The recovered source's ``exact_mult`` evaluates permutation traces in floating
arithmetic and rounds them to the nearest integer under its stated tolerance.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "depth6_recovered_source"
sys.path.insert(0, str(SOURCE))

import bqg_depth6_41 as Z  # noqa: E402
import bqg_depth6_generic_jucys_selector as J  # noqa: E402
import bqg_radial_return_reachability_gate as Reach  # noqa: E402


EXPECTED_MODULE_HASHES = {
    "bqg_block_apply_invariants.py": "199a73d45645981f6a6c195bce8b595519172a71a0555004d14b5b7a4afb39f5",
    "bqg_covariant_volume_recoupling.py": "45ab0a5befc51ad8323e3c9d658662318bcb2d58460387e8fe0b8b1f6e1afba4",
    "bqg_depth6_generic_master_row.py": "0435668bb05595598a020fbb1b4dde471b6cb09b634fbe078a5976336f8b505e",
    "bqg_reduced_transfer_hamiltonian.py": "7179b328431a31ffec5ce1e66ec94e97d02f78ec49323b9d4e87ab6cb881fe8f",
    "bqg_safe_boundary_master_gate.py": "9570b696e241ab1b01e5e327f86798ad6099f90b39d8db5cd5ceeac04a268c05",
}
EXPECTED_SHELL_STATES = 264_962
EXPECTED_S5_ORBITS = 2_757
EXPECTED_TRIVIAL_MULTIPLICITY = 27_227


def verify_source_manifest() -> dict[str, str]:
    manifest_path = SOURCE / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("source_commit") != "f6bfcd775b8c3cd33c6c235d48293e40657baf9d":
        raise RuntimeError("recovered source commit mismatch")
    if manifest.get("bundle_sha256") != "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32":
        raise RuntimeError("recovered source bundle hash mismatch")
    module_names = manifest.get("module_sha256", {})
    hashes = {
        name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
        for name in module_names
    }
    if len(hashes) != manifest.get("module_count"):
        raise RuntimeError("recovered source module count mismatch")
    if hashes != manifest.get("module_sha256"):
        raise RuntimeError("recovered source module hash mismatch")
    for name, digest in EXPECTED_MODULE_HASHES.items():
        if hashes.get(name) != digest:
            raise RuntimeError(f"artifact-pinned source hash mismatch: {name}")
    return hashes


def build_shell(depth: int) -> set[tuple[int, ...]]:
    current = {(1,) * 10}
    for level in range(1, depth + 1):
        nxt: set[tuple[int, ...]] = set()
        for vertex in Reach.V:
            nxt |= Reach.apply_H(current, vertex)
        _, current = Reach.summarize(nxt)
        print(f"shell {level}/{depth}: {len(current)} Gauss-admissible states", flush=True)
    return current


def canonical_s5_representatives(states: set[tuple[int, ...]]) -> list[tuple[int, ...]]:
    """Vectorized equivalent of source ``canon5`` for the bounded spin digits."""
    ordered = sorted(states)
    matrix = np.asarray(ordered, dtype=np.int16)
    if matrix.shape != (len(states), 10) or int(matrix.max(initial=0)) > 7:
        raise RuntimeError("unexpected K5 spin-state shape or digit bound")

    # Base 8 preserves lexicographic order for ten doubled-spin digits in 0..7.
    weights = np.asarray([8 ** p for p in range(9, -1, -1)], dtype=np.int64)
    best = np.full(len(ordered), np.iinfo(np.int64).max, dtype=np.int64)
    edge_index = Reach.EI
    for permutation in Z.S5:
        input_at_output = [0] * 10
        for old, (u, v) in enumerate(Reach.EDGES):
            new_edge = edge_index[tuple(sorted((permutation[u], permutation[v])))]
            input_at_output[new_edge] = old
        candidate = matrix[:, input_at_output] @ weights
        np.minimum(best, candidate, out=best)

    codes = np.unique(best)
    representatives = []
    for code in codes.tolist():
        digits = [0] * 10
        value = int(code)
        for i in range(9, -1, -1):
            digits[i] = value % 8
            value //= 8
        rep = tuple(digits)
        # Spot-check each unique code against the recovered implementation.
        if Z.canon5(rep)[0] != rep:
            raise RuntimeError(f"vectorized/source canonicalizer disagreement: {rep}")
        representatives.append(rep)
    return sorted(representatives)


def build_census(depth: int = 6) -> dict[str, object]:
    started = time.time()
    source_hashes = verify_source_manifest()
    states = build_shell(depth)
    if len(states) != EXPECTED_SHELL_STATES:
        raise RuntimeError(f"shell size mismatch: {len(states)}")

    representatives = canonical_s5_representatives(states)
    if len(representatives) != EXPECTED_S5_ORBITS:
        raise RuntimeError(f"S5 orbit count mismatch: {len(representatives)}")

    blocks = []
    multiplicity_total = 0
    multiplicity_histogram: collections.Counter[int] = collections.Counter()
    coordinate_histogram: collections.Counter[int] = collections.Counter()
    for orbit_id, rep in enumerate(representatives):
        multiplicity = J.exact_mult(tuple(rep), "5")
        if not multiplicity:
            continue
        coordinate_dim = J.layout(tuple(rep), "triv")[2]
        multiplicity_total += multiplicity
        multiplicity_histogram[multiplicity] += 1
        coordinate_histogram[coordinate_dim] += 1
        blocks.append(
            {
                "orbit_id": orbit_id,
                "representative": list(rep),
                "selector_multiplicity": multiplicity,
                "selector_coordinate_dimension": coordinate_dim,
                "work_proxy_m_times_d": multiplicity * coordinate_dim,
            }
        )

    if multiplicity_total != EXPECTED_TRIVIAL_MULTIPLICITY:
        raise RuntimeError(f"[5] multiplicity mismatch: {multiplicity_total}")
    if len(blocks) != EXPECTED_S5_ORBITS:
        raise RuntimeError(f"active orbit mismatch: {len(blocks)}")
    vacuum_block = next((b for b in blocks if b["representative"] == [0] * 10), None)
    if vacuum_block is None or vacuum_block["selector_multiplicity"] != 1:
        raise RuntimeError("all-zero vacuum orbit must be a one-dimensional [5] block")
    nonvacuum_blocks = [b for b in blocks if b is not vacuum_block]
    nonvacuum_dimension = sum(b["selector_multiplicity"] for b in nonvacuum_blocks)
    if len(nonvacuum_blocks) != 2756 or nonvacuum_dimension != 27226:
        raise RuntimeError("nonvacuum complement block census mismatch")

    return {
        "status": "STRUCTURAL_BLOCK_CENSUS_ONLY",
        "scope": "Recovered finite depth-6 spin shell; selector multiplicities are source-rounded character traces.",
        "source": {
            "repository": "Shtenco/synergy_binary_quantum_theory",
            "commit": "f6bfcd775b8c3cd33c6c235d48293e40657baf9d",
            "bundle_path": "mixed_parallel/bqg_mixed_master_bundle.tar.gz.b64",
            "bundle_sha256": "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32",
            "module_sha256": source_hashes,
        },
        "depth": depth,
        "gauss_admissible_spin_assignments": len(states),
        "s5_spin_orbits": len(representatives),
        "active_trivial_s5_blocks": len(blocks),
        "summed_trivial_selector_multiplicity": multiplicity_total,
        "multiplicity_method": "bqg_depth6_generic_jucys_selector.exact_mult; floating representation traces rounded to integer within source tolerance",
        "multiplicity_histogram": {str(k): multiplicity_histogram[k] for k in sorted(multiplicity_histogram)},
        "selector_coordinate_dimension_histogram": {str(k): coordinate_histogram[k] for k in sorted(coordinate_histogram)},
        "max_block_multiplicity": max(b["selector_multiplicity"] for b in blocks),
        "max_selector_coordinate_dimension": max(b["selector_coordinate_dimension"] for b in blocks),
        "vacuum_block": vacuum_block,
        "nonvacuum_active_blocks": len(nonvacuum_blocks),
        "nonvacuum_selector_multiplicity": nonvacuum_dimension,
        "largest_work_proxy_blocks": sorted(blocks, key=lambda b: (-b["work_proxy_m_times_d"], b["orbit_id"]))[:25],
        "blocks": blocks,
        "elapsed_seconds": round(time.time() - started, 3),
        "rank_claim": None,
        "proof_limit": "Block decomposition does not establish constraint-map injectivity on any block. Exact source operator identity and certified per-block ranks remain open.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "verification_results" / "DEPTH6_5_ORBIT_BLOCK_CENSUS_20261126.json",
    )
    args = parser.parse_args()
    result = build_census()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": result["status"],
                "shell": result["gauss_admissible_spin_assignments"],
                "orbits": result["s5_spin_orbits"],
                "active_blocks": result["active_trivial_s5_blocks"],
                "selector_dimension": result["summed_trivial_selector_multiplicity"],
                "max_block_multiplicity": result["max_block_multiplicity"],
                "max_coordinate_dim": result["max_selector_coordinate_dimension"],
                "elapsed_seconds": result["elapsed_seconds"],
                "rank_claim": result["rank_claim"],
                "artifact": str(args.output),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

