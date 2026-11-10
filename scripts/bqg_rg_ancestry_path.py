#!/usr/bin/env python3
"""Find and validate the maximum-fidelity ancestry path through BQG channels."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np


def load(paths, stochastic_tol=1e-8):
    mats = []
    for path in paths:
        item = json.loads(Path(path).read_text(encoding="utf-8"))
        mats.append(item)
    mats.sort(key=lambda item: item["source_s2"])

    for item in mats:
        F = np.asarray(item["channel_flow_matrix"], dtype=float)
        ns = len(item["source_master_eigenvalues"])
        nt = len(item["target_master_eigenvalues"])
        if F.shape != (ns, nt):
            raise ValueError(
                f"channel matrix shape {F.shape} does not match "
                f"source/target multiplicities {(ns, nt)}"
            )
        if not np.isfinite(F).all():
            raise ValueError("channel-flow matrix contains non-finite entries")
        if np.min(F) < -1e-10 or np.max(F) > 1.0 + 1e-10:
            raise ValueError("channel-flow entries lie outside [0,1]")
        row_sums = F.sum(axis=1)
        defect = float(np.max(np.abs(row_sums - 1.0)))
        if defect > stochastic_tol:
            raise ValueError(
                f"channel-flow matrix is not row-stochastic: defect={defect}"
            )
        if item.get("passed") is False:
            raise ValueError("a channel-flow shard did not pass its own gate")

    for left, right in zip(mats, mats[1:]):
        if left["target_s2"] != right["source_s2"]:
            raise ValueError("nonconsecutive flow matrices")
    return mats


def solve(mats):
    if not mats:
        raise ValueError("at least one channel-flow matrix is required")

    n0 = len(mats[0]["source_master_eigenvalues"])
    # Neutral prior if the initial source itself has multiplicity.
    score = np.full(n0, -math.log(n0), dtype=float)
    paths = [[i] for i in range(n0)]
    stages = []

    for item in mats:
        F = np.asarray(item["channel_flow_matrix"], dtype=float)
        if F.shape[0] != len(score):
            raise ValueError("matrix/source score shape mismatch")

        new_score = np.full(F.shape[1], -np.inf, dtype=float)
        new_paths = [None] * F.shape[1]
        parent = [None] * F.shape[1]

        for target in range(F.shape[1]):
            candidates = [
                score[source] + math.log(F[source, target])
                if F[source, target] > 0.0 else -np.inf
                for source in range(F.shape[0])
            ]
            source = int(np.argmax(candidates))
            if not math.isfinite(candidates[source]):
                raise ValueError(f"target channel {target} has no reachable ancestry")
            new_score[target] = candidates[source]
            new_paths[target] = paths[source] + [target]
            parent[target] = source

        stages.append({
            "source_j": item["source_j"],
            "target_j": item["target_j"],
            "flow_matrix": F.tolist(),
            "best_parent_for_target": parent,
            "best_log_fidelity_to_target": [float(x) for x in new_score],
            "row_sums": [float(x) for x in F.sum(axis=1)],
        })
        score, paths = new_score, new_paths

    end = int(np.argmax(score))
    path = paths[end]
    best_log = float(score[end])
    best_product = float(math.exp(best_log))

    steps = []
    for index, item in enumerate(mats):
        source = path[index]
        target = path[index + 1]
        F = np.asarray(item["channel_flow_matrix"], dtype=float)
        overlap = float(F[source, target])
        steps.append({
            "source_j": item["source_j"],
            "target_j": item["target_j"],
            "source_channel": int(source),
            "target_channel": int(target),
            "source_master_eigenvalue": float(item["source_master_eigenvalues"][source]),
            "target_master_eigenvalue": float(item["target_master_eigenvalues"][target]),
            "overlap": overlap,
            "mismatch": float(math.sqrt(max(0.0, 1.0 - overlap))),
        })

    # Explicit baseline for the rejected shortcut: always choose channel index 0.
    lowest_path = [0] * (len(mats) + 1)
    lowest_steps = []
    lowest_product = 1.0
    for index, item in enumerate(mats):
        F = np.asarray(item["channel_flow_matrix"], dtype=float)
        if F.shape[0] == 0 or F.shape[1] == 0:
            raise ValueError("empty channel multiplicity")
        overlap = float(F[0, 0])
        lowest_product *= overlap
        lowest_steps.append({
            "source_j": item["source_j"],
            "target_j": item["target_j"],
            "source_channel": 0,
            "target_channel": 0,
            "overlap": overlap,
        })

    return {
        "status": "maximum-fidelity ancestry path through BQG multiplicity channels",
        "channel_path": [int(x) for x in path],
        "best_log_fidelity": best_log,
        "best_cumulative_fidelity": best_product,
        "steps": steps,
        "naive_lowest_channel_path": lowest_path,
        "naive_lowest_cumulative_fidelity": float(lowest_product),
        "ancestry_to_lowest_product_ratio": (
            float(best_product / lowest_product) if lowest_product > 0 else None
        ),
        "naive_lowest_steps": lowest_steps,
        "stages": stages,
        "claim_boundary": (
            "This optimizes representation-RG projector overlap, not physical "
            "probability, physical time, or action. The result is a finite "
            "microscopic ancestry criterion and must be followed by the "
            "graph-changing habitat residual and physical-projector tests."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(load(args.inputs))
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
