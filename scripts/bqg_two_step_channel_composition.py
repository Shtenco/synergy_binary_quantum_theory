#!/usr/bin/env python3
"""Direct two-step BQG channel composition versus scalar RG-matrix product.

The script composes the actual generated rank-2 [2,2] carrier subspace through
two consecutive microscopic binary-blocking tensors. It compares the direct
endpoint channel overlap with F_j @ F_{j+1}, which only retains one-step scalar
channel weights. No operator, cutoff, channel ordering, or tolerance is fitted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import bqg_full_channel_flow_matrix as FLOW


def projectors(master):
    return [channel @ channel.conj().T for channel in master["channels"]]


def flow_matrix(T, source, target):
    target_projectors = projectors(target)
    rows = []
    singular_values = []
    ranks = []
    for source_channel in source["channels"]:
        Pblock, svals, rank = FLOW.generated_copy(
            T, source_channel, target["P22_basis"]
        )
        singular_values.append([float(value) for value in svals])
        ranks.append(int(rank))
        rows.append([
            float(np.trace(Pblock @ Ptarget).real / 2.0)
            for Ptarget in target_projectors
        ])
    return np.asarray(rows, dtype=float), singular_values, ranks


def orthonormal_range(projector, expected_rank, context):
    values, vectors = np.linalg.eigh((projector + projector.conj().T) / 2)
    basis = vectors[:, values > 0.5]
    if basis.shape[1] != expected_rank:
        raise RuntimeError(
            f"{context}: expected rank {expected_rank}, got {basis.shape[1]}"
        )
    defect = float(np.linalg.norm(basis.conj().T @ basis - np.eye(expected_rank)))
    if defect >= 1e-10:
        raise RuntimeError(f"{context}: extracted range is not orthonormal: {defect}")
    return basis


def block_diagnostics(Pimage, target_projectors):
    off_diagonal = []
    internal_anisotropy = []
    weights = []
    for b, Pb in enumerate(target_projectors):
        block = Pb @ Pimage @ Pb
        weight = float(np.trace(block).real / 2.0)
        weights.append(weight)
        scalar_block = weight * Pb
        internal_anisotropy.append(
            float(np.linalg.norm(block - scalar_block, ord="fro"))
        )
        for bp, Pbp in enumerate(target_projectors):
            if b != bp:
                off_diagonal.append({
                    "b": b,
                    "bp": bp,
                    "frobenius_norm": float(
                        np.linalg.norm(Pb @ Pimage @ Pbp, ord="fro")
                    ),
                })
    return {
        "weights": weights,
        "off_diagonal_blocks": off_diagonal,
        "max_off_diagonal_frobenius": max(
            (item["frobenius_norm"] for item in off_diagonal), default=0.0
        ),
        "internal_channel_anisotropy_frobenius": internal_anisotropy,
        "max_internal_channel_anisotropy": max(internal_anisotropy, default=0.0),
    }


def run_pair(source_s2):
    if source_s2 not in (3, 4):
        raise ValueError("supported two-step source s2 values are 3 and 4")

    FLOW.ZVM.patch_and_clear()
    source = FLOW.master_channels(source_s2)
    middle = FLOW.master_channels(source_s2 + 1)
    target = FLOW.master_channels(source_s2 + 2)

    T1, leakage1 = FLOW.refinement_tensor(source_s2)
    T2, leakage2 = FLOW.refinement_tensor(source_s2 + 1)

    F1, first_singulars, first_ranks = flow_matrix(T1, source, middle)
    F2, second_singulars, second_ranks = flow_matrix(T2, middle, target)
    scalar_product = F1 @ F2
    max_product_path = np.zeros_like(scalar_product)
    max_product_middle_channel = np.zeros(scalar_product.shape, dtype=int)
    for source_index in range(F1.shape[0]):
        for target_index in range(F2.shape[1]):
            path_products = F1[source_index, :] * F2[:, target_index]
            middle_index = int(np.argmax(path_products))
            max_product_middle_channel[source_index, target_index] = middle_index
            max_product_path[source_index, target_index] = path_products[middle_index]

    saved_first_path = (
        Path(__file__).resolve().parents[1]
        / "docs/numerical_continuation"
        / f"BQG_FULL_CHANNEL_FLOW_S2_{source_s2}.json"
    )
    saved_second_path = (
        Path(__file__).resolve().parents[1]
        / "docs/numerical_continuation"
        / f"BQG_FULL_CHANNEL_FLOW_S2_{source_s2 + 1}.json"
    )
    saved_first = json.loads(saved_first_path.read_text(encoding="utf-8"))
    saved_second = json.loads(saved_second_path.read_text(encoding="utf-8"))
    if not saved_first.get("passed") or not saved_second.get("passed"):
        raise RuntimeError("one-step reference matrix is not explicitly passed")
    saved_F1 = np.asarray(saved_first["channel_flow_matrix"], dtype=float)
    saved_F2 = np.asarray(saved_second["channel_flow_matrix"], dtype=float)
    if saved_F1.shape != F1.shape or np.max(np.abs(saved_F1 - F1)) >= 1e-8:
        raise RuntimeError("recomputed first one-step matrix differs from persisted reference")
    if saved_F2.shape != F2.shape or np.max(np.abs(saved_F2 - F2)) >= 1e-8:
        raise RuntimeError("recomputed second one-step matrix differs from persisted reference")

    middle_projectors = projectors(middle)
    target_projectors = projectors(target)

    direct_rows = []
    direct_first_singulars = []
    direct_second_singulars = []
    direct_first_ranks = []
    direct_second_ranks = []
    channel_blocks = []

    for source_index, source_channel in enumerate(source["channels"]):
        P1, s1, rank1 = FLOW.generated_copy(
            T1, source_channel, middle["P22_basis"]
        )
        Q1 = orthonormal_range(P1, 2, f"source channel {source_index}, first step")

        P2, s2, rank2 = FLOW.generated_copy(T2, Q1, target["P22_basis"])
        row = [
            float(np.trace(P2 @ Ptarget).real / 2.0)
            for Ptarget in target_projectors
        ]
        direct_rows.append(row)
        direct_first_singulars.append([float(x) for x in s1])
        direct_second_singulars.append([float(x) for x in s2])
        direct_first_ranks.append(int(rank1))
        direct_second_ranks.append(int(rank2))
        channel_blocks.append(
            block_diagnostics(P1, middle_projectors)
        )

    direct = np.asarray(direct_rows, dtype=float)
    defect = direct - scalar_product

    def row_error(matrix):
        return float(np.max(np.abs(np.sum(matrix, axis=1) - 1.0)))

    finite = all(
        np.isfinite(matrix).all()
        for matrix in (F1, F2, scalar_product, direct, defect)
    )
    structural_checks_passed = bool(
        finite
        and leakage1 < 1e-8
        and leakage2 < 1e-8
        and row_error(F1) < 1e-8
        and row_error(F2) < 1e-8
        and row_error(direct) < 1e-8
        and all(rank == 2 for rank in first_ranks + second_ranks)
        and all(rank == 2 for rank in direct_first_ranks + direct_second_ranks)
    )
    factorization_tolerance = 1e-8
    scalar_factorization_passed = bool(
        float(np.max(np.abs(defect))) < factorization_tolerance
    )

    result = {
        "status": "COMPUTED" if structural_checks_passed else "STRUCTURAL_CHECK_FAILURE",
        "structural_checks_passed": structural_checks_passed,
        "scalar_factorization_passed": scalar_factorization_passed,
        "scalar_factorization_tolerance": factorization_tolerance,
        "source_one_step_reference_sha256": hashlib.sha256(
            saved_first_path.read_bytes()
        ).hexdigest(),
        "middle_one_step_reference_sha256": hashlib.sha256(
            saved_second_path.read_bytes()
        ).hexdigest(),
        "source_s2": source_s2,
        "middle_s2": source_s2 + 1,
        "target_s2": source_s2 + 2,
        "source_master_eigenvalues": source["eigenvalues"],
        "middle_master_eigenvalues": middle["eigenvalues"],
        "target_master_eigenvalues": target["eigenvalues"],
        "refinement_tensor_leakage_first": float(leakage1),
        "refinement_tensor_leakage_second": float(leakage2),
        "F_first": F1.tolist(),
        "F_second": F2.tolist(),
        "F_scalar_product": scalar_product.tolist(),
        "max_product_local_path_score_matrix": max_product_path.tolist(),
        "max_product_intermediate_channel_index": max_product_middle_channel.tolist(),
        "G_direct_composed": direct.tolist(),
        "composition_defect": defect.tolist(),
        "composition_defect_max_abs": float(np.max(np.abs(defect))),
        "composition_defect_frobenius": float(np.linalg.norm(defect, ord="fro")),
        "composition_defect_relative_frobenius": float(
            np.linalg.norm(defect, ord="fro")
            / max(np.linalg.norm(direct, ord="fro"), 1e-300)
        ),
        "row_stochasticity_error_F_first": row_error(F1),
        "row_stochasticity_error_F_second": row_error(F2),
        "row_stochasticity_error_G_direct": row_error(direct),
        "first_step_generated_singular_values": first_singulars,
        "second_step_generated_singular_values": second_singulars,
        "direct_first_step_singular_values": direct_first_singulars,
        "direct_second_step_singular_values": direct_second_singulars,
        "first_step_generated_ranks": first_ranks,
        "second_step_generated_ranks": second_ranks,
        "direct_first_step_generated_ranks": direct_first_ranks,
        "direct_second_step_generated_ranks": direct_second_ranks,
        "intermediate_channel_block_diagnostics": channel_blocks,
        "interpretation": (
            "G_direct composes the actual generated rank-2 carrier subspace "
            "through both microscopic refinement tensors. F_first @ F_second "
            "is the scalar channel-matrix composition. Their difference "
            "measures information discarded by replacing the intermediate "
            "carrier with channel weights."
        ),
        "claim_boundary": (
            "Finite two-step representation-RG composition test only. "
            "Not a physical transition probability or continuum theorem."
        ),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-s2", type=int, choices=(3, 4), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = run_pair(args.source_s2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["structural_checks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
