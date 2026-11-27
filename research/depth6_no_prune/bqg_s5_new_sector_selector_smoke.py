#!/usr/bin/env python3
"""Validate the three newly censused S5 selector sectors before full sharding."""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research" / "depth6_recovered_source"))
sys.path.insert(0, str(ROOT / "research" / "depth6_no_prune"))

from bqg_s5_full_irrep_census import shell
import bqg_depth6_41 as Z
import bqg_depth6_generic_jucys_selector as J

EXPECTED = {"5": (2757, 27227), "41": (2755, 104146), "11111": (2656, 26794)}
SEEDS = (17, 12345, 20261127)
SAMPLE_BLOCKS = 12


def main() -> int:
    states = shell(6)
    if len(states) != 264962:
        raise RuntimeError(f"depth-6 shell mismatch: {len(states)}")
    reps = sorted({Z.canon5(tuple(state))[0] for state in states})
    if len(reps) != 2757:
        raise RuntimeError(f"S5 orbit mismatch: {len(reps)}")
    result = {
        "status": "passed",
        "scope": "complete multiplicity census plus deterministic multi-seed selector smoke for [5], [4,1], [1^5]",
        "shell_states": len(states),
        "s5_orbits": len(reps),
        "seeds": list(SEEDS),
        "sample_blocks_per_sector": SAMPLE_BLOCKS,
        "sectors": {},
    }
    for key, (expected_blocks, expected_columns) in EXPECTED.items():
        active, columns = [], 0
        for orbit_index, rep in enumerate(reps):
            multiplicity = J.exact_mult(tuple(rep), key)
            if multiplicity:
                active.append((orbit_index, tuple(rep), multiplicity))
                columns += multiplicity
        if len(active) != expected_blocks or columns != expected_columns:
            raise RuntimeError(
                f"{key} census mismatch: blocks={len(active)} columns={columns}; "
                f"expected {expected_blocks}/{expected_columns}"
            )
        tested = []
        for orbit_index, rep, multiplicity in active[:SAMPLE_BLOCKS]:
            dim = J.layout(rep, J.CFG[key][0])[2]
            for seednum in SEEDS:
                Q = J.fast_selector(rep, key, extra=2, seednum=seednum)
                if Q.shape != (dim, multiplicity):
                    raise RuntimeError(
                        f"{key} shape mismatch at orbit {orbit_index}: {Q.shape}, "
                        f"expected {(dim, multiplicity)}"
                    )
            tested.append({
                "orbit_index": orbit_index,
                "representative": list(rep),
                "multiplicity": multiplicity,
                "coordinate_dimension": dim,
                "seeds_passed": list(SEEDS),
            })
        result["sectors"][key] = {
            "active_blocks": len(active),
            "multiplicity_columns": columns,
            "sampled_selector_blocks": len(tested),
            "selector_seed_runs": len(tested) * len(SEEDS),
            "samples": tested,
        }
        print(f"PASS {key}: {len(active)} blocks / {columns} columns; "
              f"{len(tested) * len(SEEDS)} selector runs", flush=True)
    output = Path("verification_results/S5_NEW_SECTOR_SELECTOR_SMOKE.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
