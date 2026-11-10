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
import math
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


def direct_composition_matrix(T1, T2, source, middle, target):
    """Recompute direct endpoint projector overlaps without scalar-channel averaging."""
    target_projectors = projectors(target)
    rows = []
    for source_index, source_channel in enumerate(source["channels"]):
        P1, _, rank1 = FLOW.generated_copy(T1, source_channel, middle["P22_basis"])
        if rank1 != 2:
            raise RuntimeError(f"source channel {source_index}: first image rank is not two")
        Q1 = orthonormal_range(P1, 2, f"perturbed source channel {source_index}")
        P2, _, rank2 = FLOW.generated_copy(T2, Q1, target["P22_basis"])
        if rank2 != 2:
            raise RuntimeError(f"source channel {source_index}: endpoint image rank is not two")
        rows.append([float(np.trace(P2 @ Ptar).real / 2.0) for Ptar in target_projectors])
    return np.asarray(rows, dtype=float)


def compose_refinement_tensors(T1, T2, P_middle):
    """Return T2 P_middle T1 with both binary ancillas retained."""
    T1 = np.asarray(T1, dtype=complex)
    T2 = np.asarray(T2, dtype=complex)
    P_middle = np.asarray(P_middle, dtype=complex)
    if T1.ndim != 3 or T2.ndim != 3 or P_middle.ndim != 2:
        raise ValueError("refinement tensors must be rank-3 and middle projector rank-2")
    if T1.shape[0] != P_middle.shape[0] or P_middle.shape[0] != P_middle.shape[1]:
        raise ValueError("middle projector and first tensor output dimensions do not match")
    if T2.shape[1] != P_middle.shape[0]:
        raise ValueError("middle projector and second tensor input dimensions do not match")
    if not all(np.isfinite(x).all() for x in (T1.real, T1.imag, T2.real, T2.imag, P_middle.real, P_middle.imag)):
        raise ValueError("non-finite tensor or intermediate projector")
    herm = float(np.linalg.norm(P_middle - P_middle.conj().T, ord="fro"))
    idem = float(np.linalg.norm(P_middle @ P_middle - P_middle, ord="fro"))
    if herm > 1e-10 or idem > 1e-10:
        raise ValueError("intermediate carrier operator is not an orthogonal projector")
    return np.einsum("rmb,mn,nca->rcab", T2, P_middle, T1, optimize=True)


def direct_microscopic_endpoint_channel_weights(T_composed, source, target, context):
    """Contract both microscopic tensors before channel normalization.

    Returns norm fractions of the composed source-channel carrier in each
    target master channel, after the explicitly declared intermediate P22
    projection. These are finite two-step Hilbert-space norm weights, not
    physical-time probabilities and not a full six-level endpoint overlap.
    """
    target_projectors = projectors(target)
    target_basis = target["P22_basis"]
    rows = []
    total_norms = []
    completeness_defects = []
    for source_index, source_channel in enumerate(source["channels"]):
        # T_composed[r,c,a,b], source_channel[c,u] -> endpoint carrier
        # G[r,u,a,b], retaining both ancilla indices and source multiplicity.
        G = np.einsum("rcab,cu->ruab", T_composed, source_channel, optimize=True)
        G = G.reshape(G.shape[0], -1)
        G22 = target_basis @ (target_basis.conj().T @ G)
        total = float(np.linalg.norm(G22, ord="fro") ** 2)
        if not math.isfinite(total) or total <= 1e-300:
            raise RuntimeError(f"{context}: source channel {source_index} has zero/nonfinite endpoint norm")
        weights = []
        for target_index, Ptarget in enumerate(target_projectors):
            projected = Ptarget @ G22
            weight = float(np.linalg.norm(projected, ord="fro") ** 2 / total)
            if not math.isfinite(weight) or weight < -1e-10 or weight > 1.0 + 1e-8:
                raise RuntimeError(f"{context}: invalid endpoint channel weight")
            weights.append(weight)
        row_sum = float(sum(weights))
        defect = abs(row_sum - 1.0)
        rows.append(weights)
        total_norms.append(total)
        completeness_defects.append(defect)
    return {
        "weights": np.asarray(rows, dtype=float),
        "endpoint_norm_squared_by_source": total_norms,
        "row_sums": [float(np.sum(row)) for row in rows],
        "row_stochasticity_error": max(completeness_defects, default=0.0),
        "semantics": (
            "Two-step amplitude-level tensor contraction with intermediate P22 projection; "
            "norm fractions over target master channels. Not physical-time probability "
            "and not the full multi-level/global microscopic endpoint overlap."
        ),
    }


def normalized_refinement_spectral_distance(left_diag, right_diag):
    left_sv = np.asarray(left_diag["singular_values"], dtype=float)
    right_sv = np.asarray(right_diag["singular_values"], dtype=float)
    left_tol = float(left_diag["rank_tolerance"])
    right_tol = float(right_diag["rank_tolerance"])
    left = left_sv[left_sv > left_tol] ** 2
    right = right_sv[right_sv > right_tol] ** 2
    left_total = float(np.sum(left))
    right_total = float(np.sum(right))
    if (not math.isfinite(left_total) or not math.isfinite(right_total)
            or left_total <= 0.0 or right_total <= 0.0):
        raise RuntimeError("refinement tensor has zero or invalid spectral mass")
    left /= left_total
    right /= right_total
    width = max(left.size, right.size)
    left = np.pad(left, (0, width - left.size))
    right = np.pad(right, (0, width - right.size))
    return {
        "normalized_squared_singular_mass_l1": float(np.sum(np.abs(left - right))),
        "normalized_squared_singular_mass_hellinger": float(
            np.linalg.norm(np.sqrt(left) - np.sqrt(right)) / math.sqrt(2.0)
        ),
        "left_rank": int(left_diag["numerical_rank"]),
        "right_rank": int(right_diag["numerical_rank"]),
    }


def refinement_perturbation_sensitivity(T1, T2, source, middle, target, baseline, source_s2):
    rng = np.random.default_rng(20261110 + int(source_s2))
    rows = []
    for relative_eps in (1e-14, 1e-13, 1e-12):
        perturbed = []
        for T in (T1, T2):
            noise = rng.normal(size=T.shape) + 1j * rng.normal(size=T.shape)
            noise_norm = float(np.linalg.norm(noise))
            noise *= relative_eps * max(float(np.linalg.norm(T)), 1.0) / max(noise_norm, 1e-300)
            perturbed.append(T + noise)
        direct_p = direct_composition_matrix(perturbed[0], perturbed[1], source, middle, target)
        defect = direct_p - baseline
        rows.append({
            "relative_tensor_perturbation": relative_eps,
            "max_direct_overlap_change": float(np.max(np.abs(defect))),
            "frobenius_direct_overlap_change": float(np.linalg.norm(defect, ord="fro")),
            "finite": bool(np.isfinite(direct_p).all()),
        })
    return rows


def run_pair(source_s2, reference_dir=None):
    if source_s2 not in (3, 4, 5, 6, 7):
        raise ValueError("supported two-step source s2 values are 3 through 7")

    FLOW.ZVM.patch_and_clear()
    source = FLOW.master_channels(source_s2)
    middle = FLOW.master_channels(source_s2 + 1)
    target = FLOW.master_channels(source_s2 + 2)

    T1, leakage1 = FLOW.refinement_tensor(source_s2)
    T2, leakage2 = FLOW.refinement_tensor(source_s2 + 1)
    refinement_diag1 = FLOW.refinement_tensor_diagnostics(T1, f"two-step first s2={source_s2}")
    refinement_diag2 = FLOW.refinement_tensor_diagnostics(T2, f"two-step second s2={source_s2 + 1}")
    # Compose the microscopic tensors with the intermediate [2,2] carrier
    # projector inserted. This is a raw linear composition; it is NOT the
    # normalized-range map used by direct_composition_matrix below.
    # T1[middle, source, ancilla_1], T2[target, middle, ancilla_2].
    P_middle = middle["P22_basis"] @ middle["P22_basis"].conj().T
    T_composed = compose_refinement_tensors(T1, T2, P_middle)
    refinement_diag_composed = FLOW.refinement_tensor_diagnostics(
        T_composed, f"raw projected two-step refinement s2={source_s2}->{source_s2 + 2}"
    )
    refinement_spectral_step_change = normalized_refinement_spectral_distance(
        refinement_diag1, refinement_diag2
    )
    refinement_spectral_composition_change = normalized_refinement_spectral_distance(
        refinement_diag1, refinement_diag_composed
    )
    if not math.isfinite(float(leakage1)) or not math.isfinite(float(leakage2)):
        raise RuntimeError("non-finite refinement leakage in direct two-step composition")

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

    reference_root = (
        Path(reference_dir) if reference_dir is not None
        else Path(__file__).resolve().parents[1] / "docs/numerical_continuation"
    )
    saved_first_path = reference_root / f"BQG_FULL_CHANNEL_FLOW_S2_{source_s2}.json"
    saved_second_path = reference_root / f"BQG_FULL_CHANNEL_FLOW_S2_{source_s2 + 1}.json"
    saved_first = json.loads(saved_first_path.read_text(encoding="utf-8"))
    saved_second = json.loads(saved_second_path.read_text(encoding="utf-8"))
    if saved_first.get("passed") is not True or saved_second.get("passed") is not True:
        raise RuntimeError("one-step reference matrix must explicitly pass with boolean true")
    expected_scales = (
        (saved_first, source_s2, source_s2 + 1, "first"),
        (saved_second, source_s2 + 1, source_s2 + 2, "second"),
    )
    for item, expected_source, expected_target, label in expected_scales:
        if item.get("source_s2") != expected_source or item.get("target_s2") != expected_target:
            raise RuntimeError(f"{label} reference scale labels do not match requested two-step chain")
        run_id = item.get("workflow_run_id")
        commit = item.get("workflow_commit")
        if not isinstance(run_id, int) or isinstance(run_id, bool) or run_id <= 0:
            raise RuntimeError(f"{label} reference has invalid workflow_run_id")
        if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit.lower()):
            raise RuntimeError(f"{label} reference has invalid workflow_commit")
    if saved_first["workflow_run_id"] != saved_second["workflow_run_id"]:
        raise RuntimeError("one-step references come from different workflow runs")
    if saved_first["workflow_commit"] != saved_second["workflow_commit"]:
        raise RuntimeError("one-step references come from different source commits")
    saved_first_target = np.asarray(saved_first["target_master_eigenvalues"], dtype=float)
    saved_second_source = np.asarray(saved_second["source_master_eigenvalues"], dtype=float)
    if saved_first_target.shape != saved_second_source.shape or not np.isfinite(saved_first_target).all() or not np.isfinite(saved_second_source).all():
        raise RuntimeError("adjacent one-step reference spectra are malformed")
    if (not np.all(np.diff(saved_first_target) > 0.0)
        or not np.all(np.diff(saved_second_source) > 0.0)):
        raise RuntimeError("adjacent reference spectra must be strictly sorted ascending")
    spectral_pairs = (
        (np.asarray(saved_first["source_master_eigenvalues"], dtype=float), np.asarray(source["eigenvalues"], dtype=float), "first source"),
        (saved_first_target, np.asarray(middle["eigenvalues"], dtype=float), "middle"),
        (np.asarray(saved_second["target_master_eigenvalues"], dtype=float), np.asarray(target["eigenvalues"], dtype=float), "final target"),
    )
    for persisted, recomputed, label in spectral_pairs:
        if persisted.shape != recomputed.shape or not np.isfinite(persisted).all() or not np.isfinite(recomputed).all():
            raise RuntimeError(f"{label} persisted/recomputed spectrum shape or finiteness mismatch")
        local_scale = max(float(np.max(np.abs(persisted))), float(np.max(np.abs(recomputed))), 1.0)
        if float(np.max(np.abs(persisted - recomputed))) > 1e-10 * local_scale:
            raise RuntimeError(f"{label} persisted/recomputed spectrum mismatch")
    spectral_stitch_defect = float(np.max(np.abs(saved_first_target - saved_second_source)))
    spectral_scale = max(float(np.max(np.abs(saved_first_target))), float(np.max(np.abs(saved_second_source))), 1.0)
    if spectral_stitch_defect > 1e-10 * spectral_scale:
        raise RuntimeError(f"adjacent one-step reference spectra do not stitch: {spectral_stitch_defect}")
    saved_F1 = np.asarray(saved_first["channel_flow_matrix"], dtype=float)
    saved_F2 = np.asarray(saved_second["channel_flow_matrix"], dtype=float)
    if not np.isfinite(saved_F1).all() or not np.isfinite(saved_F2).all():
        raise RuntimeError("persisted one-step matrix contains NaN or infinity")
    if saved_F1.shape != F1.shape or float(np.max(np.abs(saved_F1 - F1))) >= 1e-8:
        raise RuntimeError("recomputed first one-step matrix differs from persisted reference")
    if saved_F2.shape != F2.shape or float(np.max(np.abs(saved_F2 - F2))) >= 1e-8:
        raise RuntimeError("recomputed second one-step matrix differs from persisted reference")

    middle_projectors = projectors(middle)
    target_projectors = projectors(target)

    direct_rows = []
    direct_first_singulars = []
    direct_second_singulars = []
    direct_first_ranks = []
    direct_second_ranks = []
    channel_blocks = []
    endpoint_channel_blocks = []
    basis_rotation_defects = []
    endpoint_projector_defects = []

    for source_index, source_channel in enumerate(source["channels"]):
        P1, s1, rank1 = FLOW.generated_copy(
            T1, source_channel, middle["P22_basis"]
        )
        Q1 = orthonormal_range(P1, 2, f"source channel {source_index}, first step")

        P2, s2, rank2 = FLOW.generated_copy(T2, Q1, target["P22_basis"])

        # The direct endpoint carrier must not depend on the arbitrary U(2)
        # basis chosen inside the intermediate rank-2 image subspace.
        theta, phi = 0.371, 0.193
        rotation = np.array([
            [math.cos(theta), -np.exp(1j * phi) * math.sin(theta)],
            [np.exp(-1j * phi) * math.sin(theta), math.cos(theta)],
        ], dtype=complex)
        P2_rotated, _, _ = FLOW.generated_copy(
            T2, Q1 @ rotation, target["P22_basis"]
        )
        basis_rotation_defects.append(
            float(np.linalg.norm(P2_rotated - P2, ord="fro"))
        )

        p2_hermiticity = float(np.linalg.norm(P2 - P2.conj().T, ord="fro"))
        p2_idempotence = float(np.linalg.norm(P2 @ P2 - P2, ord="fro"))
        endpoint_projector_defects.append({"hermiticity": p2_hermiticity, "idempotence": p2_idempotence})
        endpoint_channel_blocks.append(block_diagnostics(P2, target_projectors))
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
    amplitude_endpoint = direct_microscopic_endpoint_channel_weights(
        T_composed, source, target, f"two-step s2={source_s2}"
    )
    amplitude_endpoint_matrix = amplitude_endpoint["weights"]
    sensitivity = refinement_perturbation_sensitivity(T1, T2, source, middle, target, direct, source_s2)
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
        and amplitude_endpoint["row_stochasticity_error"] < 1e-8
        and np.isfinite(amplitude_endpoint_matrix).all()
        and all(rank == 2 for rank in first_ranks + second_ranks)
        and all(rank == 2 for rank in direct_first_ranks + direct_second_ranks)
        and max(basis_rotation_defects, default=0.0) < 1e-8
        and max((x["hermiticity"] for x in endpoint_projector_defects), default=0.0) < 1e-10
        and max((x["idempotence"] for x in endpoint_projector_defects), default=0.0) < 1e-10
        and all(x["finite"] for x in sensitivity)
        and max((x["max_direct_overlap_change"] for x in sensitivity), default=0.0) < 1e-8
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
        "workflow_run_id": saved_first["workflow_run_id"],
        "workflow_commit": saved_first["workflow_commit"],
        "adjacent_reference_spectrum_stitch_defect": spectral_stitch_defect,
        "flow_source_sha256": hashlib.sha256(Path(FLOW.__file__).read_bytes()).hexdigest(),
        "two_step_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_one_step_reference_sha256": hashlib.sha256(
            saved_first_path.read_bytes()
        ).hexdigest(),
        "middle_one_step_reference_sha256": hashlib.sha256(
            saved_second_path.read_bytes()
        ).hexdigest(),
        "source_s2": source_s2,
        "reference_directory": str(reference_root),
        "middle_s2": source_s2 + 1,
        "target_s2": source_s2 + 2,
        "source_master_eigenvalues": source["eigenvalues"],
        "middle_master_eigenvalues": middle["eigenvalues"],
        "target_master_eigenvalues": target["eigenvalues"],
        "refinement_tensor_leakage_first": float(leakage1),
        "refinement_tensor_leakage_second": float(leakage2),
        "first_refinement_tensor_diagnostics": refinement_diag1,
        "second_refinement_tensor_diagnostics": refinement_diag2,
        "direct_composed_refinement_tensor_shape": list(T_composed.reshape(T_composed.shape[0], -1).shape),
        "direct_composed_refinement_tensor_diagnostics": refinement_diag_composed,
        "consecutive_refinement_spectral_change": refinement_spectral_step_change,
        "one_step_to_direct_two_step_spectral_change": refinement_spectral_composition_change,
        "direct_composition_refinement_perturbation_sensitivity": sensitivity,
        "max_direct_composition_perturbation_sensitivity": max(
            (x["max_direct_overlap_change"] for x in sensitivity), default=0.0
        ),
        "F_first": F1.tolist(),
        "F_second": F2.tolist(),
        "F_scalar_product": scalar_product.tolist(),
        "max_product_local_path_score_matrix": max_product_path.tolist(),
        "max_product_intermediate_channel_index": max_product_middle_channel.tolist(),
        "G_direct_composed": direct.tolist(),
        "direct_two_step_endpoint_projector_overlap_matrix": direct.tolist(),
        "direct_two_step_amplitude_level_endpoint_channel_weights": amplitude_endpoint_matrix.tolist(),
        "direct_two_step_amplitude_level_endpoint_norm_squared_by_source": amplitude_endpoint["endpoint_norm_squared_by_source"],
        "direct_two_step_amplitude_level_endpoint_row_sums": amplitude_endpoint["row_sums"],
        "direct_two_step_amplitude_level_endpoint_row_stochasticity_error": amplitude_endpoint["row_stochasticity_error"],
        "direct_two_step_amplitude_level_endpoint_semantics": amplitude_endpoint["semantics"],
        "global_microscopic_endpoint_overlap_matrix_legacy_alias": direct.tolist(),
        "direct_endpoint_overlap_semantics": (
            "Direct two-step endpoint projector overlap obtained by composing the two "
            "microscopic refinement maps through the intermediate rank-2 carrier. "
            "This is not the six-level/global microscopic overlap and is not a "
            "physical-time probability."
        ),
        "max_product_intermediate_channel_scope": (
            "Pairwise maximum over the intermediate channel for each fixed source/target "
            "pair; not a globally consistent multi-step ancestry path."
        ),
        "scalar_markov_composition_matrix": scalar_product.tolist(),
        "local_one_step_overlap_matrices": {"source_to_middle": F1.tolist(), "middle_to_target": F2.tolist()},
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
        "intermediate_basis_rotation_invariance_defects": basis_rotation_defects,
        "max_intermediate_basis_rotation_invariance_defect": max(
            basis_rotation_defects, default=0.0
        ),
        "intermediate_channel_block_diagnostics": channel_blocks,
        "direct_endpoint_channel_block_diagnostics": endpoint_channel_blocks,
        "direct_endpoint_projector_defects": endpoint_projector_defects,
        "max_direct_endpoint_off_diagonal_frobenius": max(
            (x["max_off_diagonal_frobenius"] for x in endpoint_channel_blocks), default=0.0
        ),
        "max_direct_endpoint_internal_channel_anisotropy": max(
            (x["max_internal_channel_anisotropy"] for x in endpoint_channel_blocks), default=0.0
        ),
        "interpretation": (
            "G_direct_composed is the directly contracted two-step microscopic "
            "endpoint projector overlap, computed by propagating the actual "
            "generated rank-2 carrier subspace through both refinement tensors. "
            "F_first and F_second are local one-step channel overlaps; their "
            "matrix product is a scalar Markov-style composition, not the "
            "microscopic endpoint overlap. max_product_local_path_score is "
            "the best single intermediate-channel path and is distinct from "
            "both the sum over intermediate channels and G_direct_composed."
        ),
        "overlap_semantics": {
            "F_first_F_second": "local one-step overlaps and their scalar matrix product",
            "max_product_local_path_score_matrix": "maximum product over one intermediate channel; not a summed or direct microscopic overlap",
            "G_direct_composed": "direct microscopic two-step endpoint projector overlap",
            "raw_projected_refinement_tensor_spectrum": "diagnostic of T2 P_middle T1; not the normalized-range carrier map used to form G_direct_composed",
        },
        "claim_boundary": (
            "Finite two-step representation-RG composition test only. "
            "Not a physical transition probability or continuum theorem."
        ),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-s2", type=int, choices=(3, 4, 5, 6, 7), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path)
    args = parser.parse_args()

    result = run_pair(args.source_s2, args.reference_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["structural_checks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
