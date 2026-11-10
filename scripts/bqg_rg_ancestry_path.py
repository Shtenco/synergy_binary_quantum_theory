#!/usr/bin/env python3
"""Find and validate the maximum-fidelity ancestry path through BQG channels."""
from __future__ import annotations
import argparse
import json
import math
import hashlib
from pathlib import Path

import numpy as np


def validate_refinement_spectrum(item, context):
    """Cross-check singular and Gram spectra of the actual induced tensor."""
    diag = item.get("refinement_tensor_diagnostics")
    if not isinstance(diag, dict):
        raise ValueError(f"{context}: missing refinement_tensor_diagnostics")
    shape = diag.get("shape")
    if (not isinstance(shape, list) or len(shape) != 2
            or any(not isinstance(x, int) or isinstance(x, bool) or x <= 0 for x in shape)):
        raise ValueError(f"{context}: malformed refinement tensor shape")
    def numeric_vector(field):
        raw = diag.get(field)
        if (not isinstance(raw, list)
                or any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in raw)):
            raise ValueError(f"{context}: {field} must be a JSON numeric array")
        return np.asarray(raw, dtype=float)

    singular = numeric_vector("singular_values")
    row_gram = numeric_vector("row_gram_eigenvalues")
    col_gram = numeric_vector("column_gram_nonzero_eigenvalues")
    rank_tol_raw = diag.get("rank_tolerance")
    frobenius_raw = diag.get("frobenius_norm")
    if (isinstance(rank_tol_raw, bool) or not isinstance(rank_tol_raw, (int, float))
            or isinstance(frobenius_raw, bool) or not isinstance(frobenius_raw, (int, float))):
        raise ValueError(f"{context}: rank tolerance and Frobenius norm must be numeric")
    rank_tol = float(rank_tol_raw)
    rank = diag.get("numerical_rank")
    frobenius = float(frobenius_raw)
    row_identity_defect = float(diag.get("row_gram_identity_defect", math.nan))
    column_idempotence_defect = float(diag.get("column_gram_idempotence_defect", math.nan))
    if (singular.ndim != 1 or singular.size != min(shape)
            or row_gram.shape != (shape[0],) or col_gram.ndim != 1
            or not np.isfinite(singular).all() or not np.isfinite(row_gram).all()
            or not np.isfinite(col_gram).all() or not math.isfinite(rank_tol)
            or rank_tol <= 0.0 or not math.isfinite(frobenius) or frobenius < 0.0
            or not math.isfinite(row_identity_defect) or row_identity_defect < 0.0
            or not math.isfinite(column_idempotence_defect) or column_idempotence_defect < 0.0):
        raise ValueError(f"{context}: malformed/non-finite refinement spectra")
    if np.any(singular < 0.0) or np.any(row_gram < -1e-10) or np.any(col_gram < -1e-10):
        raise ValueError(f"{context}: negative singular/Gram spectrum")
    if singular.size > 1 and np.any(np.diff(singular) > 1e-12 * max(float(singular[0]), 1.0)):
        raise ValueError(f"{context}: singular values are not descending")
    if frobenius <= 0.0 or float(np.sum(singular * singular)) <= 0.0:
        raise ValueError(f"{context}: refinement tensor has zero spectral mass")
    expected_rank = int(np.sum(singular > rank_tol))
    if expected_rank == 0:
        raise ValueError(f"{context}: induced refinement tensor has zero numerical rank")
    if isinstance(rank, bool) or not isinstance(rank, int) or rank != expected_rank:
        raise ValueError(f"{context}: numerical rank does not match singular spectrum")
    if float(np.sum(singular * singular)) <= 0.0:
        raise ValueError(f"{context}: induced refinement tensor has zero spectral mass")
    expected_row = np.sort(np.concatenate((np.zeros(max(0, shape[0] - singular.size)),
                                           singular * singular)))
    if expected_row.shape != row_gram.shape:
        raise ValueError(f"{context}: row Gram spectrum dimension mismatch")
    scale = max(float(np.max(np.abs(expected_row))), float(np.max(np.abs(row_gram))), 1.0)
    row_defect = float(np.max(np.abs(expected_row - row_gram)))
    expected_col = np.sort((singular[singular > rank_tol]) ** 2)
    if expected_col.shape != col_gram.shape:
        raise ValueError(f"{context}: column Gram nonzero spectrum dimension mismatch")
    col_scale = max(float(np.max(np.abs(expected_col))) if expected_col.size else 0.0,
                    float(np.max(np.abs(col_gram))) if col_gram.size else 0.0, 1.0)
    col_defect = (float(np.max(np.abs(expected_col - col_gram)))
                  if expected_col.size else 0.0)
    frob_defect = abs(frobenius * frobenius - float(np.sum(singular * singular)))
    if row_defect > 1e-8 * scale or col_defect > 1e-8 * col_scale:
        raise ValueError(f"{context}: Gram spectra disagree with singular values")
    if frob_defect > 1e-8 * max(frobenius * frobenius, 1.0):
        raise ValueError(f"{context}: Frobenius norm disagrees with singular spectrum")
    return {
        "shape": shape,
        "rank": rank,
        "singular_values": singular.tolist(),
        "normalized_squared_singular_values": (
            (singular[singular > rank_tol] ** 2
             / np.sum(singular[singular > rank_tol] ** 2)).tolist()
            if expected_rank > 0 else []
        ),
        "row_gram_singular_spectrum_max_defect": row_defect,
        "column_gram_singular_spectrum_max_defect": col_defect,
        "frobenius_singular_spectrum_defect": frob_defect,
    }


def load(paths, stochastic_tol=1e-8, require_refinement_spectra=False):
    mats = []
    if not paths:
        raise ValueError("at least one channel-flow file is required")
    for path in paths:
        source_path = Path(path)
        raw = source_path.read_bytes()
        item = json.loads(raw.decode("utf-8"))
        item["_source_path"] = str(source_path)
        item["_source_sha256"] = hashlib.sha256(raw).hexdigest()
        mats.append(item)
    mats.sort(key=lambda item: item["source_s2"])

    if not math.isfinite(stochastic_tol) or stochastic_tol < 0:
        raise ValueError("stochastic_tol must be finite and nonnegative")

    for item in mats:
        source_s2 = item.get("source_s2")
        target_s2 = item.get("target_s2")
        if (not isinstance(source_s2, int) or isinstance(source_s2, bool)
                or not isinstance(target_s2, int) or isinstance(target_s2, bool)):
            raise ValueError("source_s2 and target_s2 must be integer spin labels")
        if target_s2 != source_s2 + 1:
            raise ValueError(
                f"each refinement step must advance s2 by one: "
                f"{source_s2} -> {target_s2}"
            )
        if item.get("passed") is not True:
            raise ValueError("a channel-flow shard must explicitly pass its own gate")
        run_id = item.get("workflow_run_id")
        commit = item.get("workflow_commit")
        if not isinstance(run_id, int) or isinstance(run_id, bool) or run_id <= 0:
            raise ValueError("workflow_run_id must be a positive integer")
        if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit.lower()):
            raise ValueError("workflow_commit must be a full 40-character hexadecimal SHA")
        source_eigenvalues = np.asarray(item.get("source_master_eigenvalues"), dtype=float)
        target_eigenvalues = np.asarray(item.get("target_master_eigenvalues"), dtype=float)
        if (source_eigenvalues.ndim != 1 or target_eigenvalues.ndim != 1
                or source_eigenvalues.size == 0 or target_eigenvalues.size == 0):
            raise ValueError("source and target master spectra must be nonempty vectors")
        if not np.isfinite(source_eigenvalues).all() or not np.isfinite(target_eigenvalues).all():
            raise ValueError("source or target master spectrum contains non-finite values")
        if not np.all(np.diff(source_eigenvalues) > 0.0) or not np.all(np.diff(target_eigenvalues) > 0.0):
            raise ValueError("source and target eigenvalues must be strictly sorted ascending")
        if item.get("refinement_tensor_diagnostics") is None:
            if require_refinement_spectra:
                raise ValueError("missing refinement_tensor_diagnostics")
            item["_refinement_spectrum_validation"] = None
        else:
            item["_refinement_spectrum_validation"] = validate_refinement_spectrum(
                item, f"s2={source_s2}"
            )
        if not math.isclose(float(item.get("source_j", math.nan)), source_s2 / 2.0,
                            rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("source_j does not match source_s2")
        if not math.isclose(float(item.get("target_j", math.nan)), target_s2 / 2.0,
                            rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("target_j does not match target_s2")

        F = np.asarray(item["channel_flow_matrix"], dtype=float)
        ns = len(source_eigenvalues)
        nt = len(target_eigenvalues)
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
    for left, right in zip(mats, mats[1:]):
        if left["target_s2"] != right["source_s2"]:
            raise ValueError("nonconsecutive flow matrices")
        if left["workflow_run_id"] != right["workflow_run_id"]:
            raise ValueError("workflow_run_id changes within the ancestry chain")
        if left["workflow_commit"] != right["workflow_commit"]:
            raise ValueError("workflow_commit changes within the ancestry chain")
        left_spectrum = np.asarray(left["target_master_eigenvalues"], dtype=float)
        right_spectrum = np.asarray(right["source_master_eigenvalues"], dtype=float)
        if left_spectrum.shape != right_spectrum.shape:
            raise ValueError("adjacent target/source spectrum dimensions do not match")
        scale = max(float(np.max(np.abs(left_spectrum))), float(np.max(np.abs(right_spectrum))), 1.0)
        spectral_defect = float(np.max(np.abs(left_spectrum - right_spectrum)))
        if spectral_defect > 1e-10 * scale:
            raise ValueError(f"adjacent target/source spectra are discontinuous: defect={spectral_defect}")
        left["_spectral_stitch_defect_to_next"] = spectral_defect
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

    refinement_evolution = []
    for index, item in enumerate(mats):
        current = item.get("_refinement_spectrum_validation")
        if current is None:
            refinement_evolution.append({
                "source_s2": item["source_s2"],
                "target_s2": item["target_s2"],
                "spectrum_validation": "NOT_AVAILABLE_IN_LEGACY_SHARD",
            })
            continue
        entry = {
            "source_s2": item["source_s2"],
            "target_s2": item["target_s2"],
            **current,
        }
        if index + 1 < len(mats):
            nxt = mats[index + 1].get("_refinement_spectrum_validation")
            if nxt is not None:
                a = np.asarray(current["normalized_squared_singular_values"], dtype=float)
                b = np.asarray(nxt["normalized_squared_singular_values"], dtype=float)
                width = max(a.size, b.size)
                ap = np.pad(a, (0, width - a.size))
                bp = np.pad(b, (0, width - b.size))
                entry["next_scale_normalized_singular_mass_l1"] = float(np.sum(np.abs(ap - bp)))
                entry["next_scale_normalized_singular_mass_hellinger"] = float(
                    np.linalg.norm(np.sqrt(ap) - np.sqrt(bp)) / math.sqrt(2.0)
                )
        refinement_evolution.append(entry)

    selected_path_overlap_product = float(math.prod(step["overlap"] for step in steps))
    local_maximum_per_step = [
        float(np.max(np.asarray(item["channel_flow_matrix"], dtype=float)))
        for item in mats
    ]
    independent_local_maxima_product = float(math.prod(local_maximum_per_step))
    initial_uniform_prior_weight = 1.0 / n0

    return {
        "status": "maximum-fidelity ancestry path through BQG multiplicity channels",
        "provenance": {
            "workflow_run_id": mats[0]["workflow_run_id"],
            "workflow_commit": mats[0]["workflow_commit"],
            "input_files": [{"path": x["_source_path"], "sha256": x["_source_sha256"]} for x in mats],
            "spectral_stitch_defects": [x.get("_spectral_stitch_defect_to_next", 0.0) for x in mats[:-1]],
        },
        "refinement_spectrum_validation_complete": all(x.get("_refinement_spectrum_validation") is not None for x in mats),
        "refinement_spectral_evolution": refinement_evolution,
        "channel_path": [int(x) for x in path],
        "best_log_fidelity": best_log,
        "best_cumulative_fidelity": best_product,
        "best_stepwise_ancestry_path_score": best_product,
        "selected_path_local_projector_overlap_product": selected_path_overlap_product,
        "initial_uniform_prior_weight": initial_uniform_prior_weight,
        "selected_path_joint_weight_with_prior": float(
            initial_uniform_prior_weight * selected_path_overlap_product
        ),
        "independent_local_maximum_per_step": local_maximum_per_step,
        "product_of_independent_local_maxima_not_a_single_path": independent_local_maxima_product,
        "local_one_step_max_overlap_by_source": [
            [float(x) for x in np.max(np.asarray(item["channel_flow_matrix"], dtype=float), axis=1)]
            for item in mats
        ],
        "steps": steps,
        "global_microscopic_endpoint_overlap": "NOT_COMPUTED_BY_STEPWISE_PATH_DP",
        "overlap_semantics": (
            "Each steps[].overlap is a local one-step representation-RG projector overlap. "
            "best_stepwise_ancestry_path_score is the maximum product along a discrete "
            "sequence of these local overlaps. It is not a directly contracted global "
            "microscopic endpoint overlap; that requires composing the microscopic "
            "refinement tensors and retaining intermediate carrier coherence."
        ),
        "naive_lowest_channel_path": lowest_path,
        "naive_lowest_cumulative_fidelity": float(lowest_product),
        "best_path_score_to_fixed_channel0_score_ratio": (
            float(best_product / lowest_product) if lowest_product > 0 else None
        ),
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
    parser.add_argument("--require-refinement-spectra", action="store_true")
    args = parser.parse_args()
    result = solve(load(args.inputs, require_refinement_spectra=args.require_refinement_spectra))
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
