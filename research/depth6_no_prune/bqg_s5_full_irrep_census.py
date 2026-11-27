#!/usr/bin/env python3
"""Independent seven-sector S5 multiplicity census on the complete depth-6 shell."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "research" / "depth6_recovered_source"
sys.path.insert(0, str(SOURCE_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import bqg_radial_return_reachability_gate as Reach  # noqa: E402
import bqg_depth6_41 as Z  # noqa: E402
import bqg_depth6_generic_jucys_selector as J  # noqa: E402
import bqg_s5_irrep_decomposition_gate as D  # noqa: E402

EXPECTED_SHELL_STATES = 264_962
EXPECTED_S5_ORBITS = 2_757
EXPECTED_HILBERT_DIMENSION = 3_111_637


def shell(depth: int = 6) -> tuple[tuple[int, ...], ...]:
    current = {(1,) * 10}
    for level in range(depth):
        nxt = set()
        for vertex in Reach.V:
            nxt |= Reach.apply_H(current, vertex)
        _, current = Reach.summarize(nxt)
        print(f"SHELL {level + 1} {len(current)}", flush=True)
    return tuple(sorted(current))


def run() -> dict:
    states = shell(6)
    errors = []
    if len(states) != EXPECTED_SHELL_STATES:
        errors.append(
            f"depth-6 shell mismatch: {len(states)} != {EXPECTED_SHELL_STATES}"
        )

    representatives = sorted({Z.canon5(tuple(state))[0] for state in states})
    if len(representatives) != EXPECTED_S5_ORBITS:
        errors.append(
            f"S5 orbit mismatch: {len(representatives)} != {EXPECTED_S5_ORBITS}"
        )

    dimensions = {
        irrep: int(D.irrep_characters(tuple(range(5)))[irrep])
        for irrep in D.IRREPS
    }
    rows = {
        irrep: {"active_blocks": 0, "multiplicity_columns": 0}
        for irrep in D.IRREPS
    }

    for orbit_index, representative in enumerate(representatives):
        for irrep in D.IRREPS:
            multiplicity = J.exact_mult(tuple(representative), irrep)
            if multiplicity < 0:
                errors.append(
                    f"negative multiplicity: orbit={orbit_index} irrep={irrep} m={multiplicity}"
                )
            if multiplicity:
                rows[irrep]["active_blocks"] += 1
                rows[irrep]["multiplicity_columns"] += multiplicity
        # Keep only one orbit's dense recoupling matrices resident at a time.
        J.X._transform_matrix.cache_clear()
        if (orbit_index + 1) % 100 == 0:
            print(
                f"CENSUS {orbit_index + 1}/{len(representatives)}",
                flush=True,
            )

    unweighted_columns = sum(
        item["multiplicity_columns"] for item in rows.values()
    )
    weighted_dimension = sum(
        dimensions[irrep] * rows[irrep]["multiplicity_columns"]
        for irrep in D.IRREPS
    )
    if weighted_dimension != EXPECTED_HILBERT_DIMENSION:
        errors.append(
            "seven-sector weighted dimension mismatch: "
            f"{weighted_dimension} != {EXPECTED_HILBERT_DIMENSION}"
        )

    all_seven_present = all(
        rows[irrep]["active_blocks"] > 0
        and rows[irrep]["multiplicity_columns"] > 0
        for irrep in D.IRREPS
    )
    if not all_seven_present:
        errors.append("one or more S5 irreps have no active blocks")

    return {
        "status": "passed" if not errors else "failed",
        "source_commit": os.environ.get("GITHUB_SHA"),
        "errors": errors,
        "scope": "all seven S5 irreps over every canonical depth-6 orbit",
        "shell_states": len(states),
        "s5_orbits": len(representatives),
        "irrep_dimensions": dimensions,
        "sectors": rows,
        "unweighted_multiplicity_columns": unweighted_columns,
        "weighted_hilbert_dimension": weighted_dimension,
        "expected_weighted_hilbert_dimension": EXPECTED_HILBERT_DIMENSION,
        "all_seven_sectors_present": all_seven_present,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
