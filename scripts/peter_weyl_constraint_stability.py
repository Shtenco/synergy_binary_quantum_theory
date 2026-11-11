#!/usr/bin/env python3
"""Compare legacy and no-prune Peter-Weyl constraint maps on one carrier.

The map C is the direct sum over constraint nodes:
    C |a> = (H_0|a>, ..., H_{m-1}|a>).
The legacy Gram diagnostics are retained for comparison, while the original
constraint map is also analyzed by streaming Householder TSQR + direct SVD.
The output Hilbert-space matrix is never materialized densely.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import heapq
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_master_projector_pilot as PILOT
import peter_weyl_streaming_svd as STREAM_SVD


def stack_columns(columns, nodes, dim):
    """Combine node outputs into a direct-sum output space for each input."""
    stacked = {}
    supports = {}
    for a in range(dim):
        out = {}
        per_node = {}
        for v in nodes:
            col = columns[(v, a)]
            per_node[str(v)] = len(col)
            for key, value in col.items():
                out[(v, key)] = complex(value)
        stacked[a] = out
        supports[str(a)] = per_node
    return stacked, supports


def subtract_stacked(left, right, dim):
    """Return right-left, retaining every exactly nonzero difference."""
    delta = {}
    for a in range(dim):
        keys = set(left[a]) | set(right[a])
        col = {}
        for key in keys:
            value = right[a].get(key, 0j) - left[a].get(key, 0j)
            if value != 0j:
                col[key] = value
        delta[a] = col
    return delta


def gram_from_columns(columns, dim):
    G = np.empty((dim, dim), dtype=complex)
    for a in range(dim):
        for b in range(dim):
            G[a, b] = PILOT.inner_sparse(columns[a], columns[b])
    return G


def spectrum_of_gram(G, gram_eigenvalue_rtol):
    """Gram-only diagnostic; its tolerance is relative to eigenvalues, not sigmas."""
    if not math.isfinite(gram_eigenvalue_rtol) or gram_eigenvalue_rtol < 0:
        raise ValueError("gram_eigenvalue_rtol must be finite and nonnegative")
    if not np.all(np.isfinite(G)):
        raise ValueError("Gram matrix contains non-finite entries")
    raw_hermiticity = PILOT.stable_complex_norm((G - G.conj().T).ravel())
    H = 0.5 * G + 0.5 * G.conj().T
    evals = np.linalg.eigvalsh(H)
    if not np.all(np.isfinite(evals)):
        raise FloatingPointError("Gram eigensolver returned non-finite eigenvalues")
    scale = float(np.max(np.abs(evals))) if len(evals) else 0.0
    normalized = np.zeros_like(evals) if scale == 0.0 else evals / scale
    psd_rtol = 100.0 * np.finfo(float).eps
    if len(evals) and float(np.min(normalized)) < -psd_rtol:
        raise ValueError(f"Gram matrix is not PSD: min normalized eigenvalue {float(np.min(normalized))}")
    # Classify in normalized units: rtol*scale can underflow to zero for
    # subnormal Gram matrices, even though the relative spectral information
    # remains representable.
    threshold = gram_eigenvalue_rtol * scale
    clipped = np.maximum(evals, 0.0)
    singular = np.sqrt(clipped)
    positive = singular[normalized > gram_eigenvalue_rtol]
    psd_tol = psd_rtol * scale
    return {
        "gram_hermiticity_defect_frobenius": raw_hermiticity,
        "gram_min_eigenvalue": float(np.min(evals)) if len(evals) else 0.0,
        "gram_max_eigenvalue": float(np.max(evals)) if len(evals) else 0.0,
        "gram_eigenvalue_scale": scale,
        "normalized_gram_eigenvalues": normalized.tolist(),
        "rank_threshold_on_gram_eigenvalues": threshold,
        "normalized_rank_threshold_on_gram_eigenvalues": float(gram_eigenvalue_rtol),
        "rank_threshold_resolves_roundoff": bool(gram_eigenvalue_rtol > psd_rtol),
        "estimated_rank": int(np.count_nonzero(normalized > gram_eigenvalue_rtol)),
        "minimum_singular_value_including_nullity": float(singular[0]) if len(singular) else 0.0,
        "smallest_positive_singular_value": float(np.min(positive)) if len(positive) else None,
        "largest_singular_value": float(singular[-1]) if len(singular) else 0.0,
        "singular_values": [float(x) for x in singular],
        "psd_tolerance": psd_tol,
        "rank_tolerance_semantics": "relative Gram eigenvalue threshold; equivalent sigma threshold is sqrt(rtol)",
        "singular_value_method": "sqrt(eigvalsh(C^dagger C)); Gram formation squares condition number",
        "gram_spectral_resolution_floor": float(psd_tol),
    }


def iter_operator_row_blocks(columns, dim, block_rows=4096):
    """Yield deterministic dense row blocks from sparse columns, never full C.

    Sparse keys are (node, output-state). Sorted key references are merged
    deterministically; only one block_rows x dim dense numerical buffer is
    allocated at a time (the sparse input dictionaries remain resident).
    """
    if not isinstance(block_rows, int) or block_rows <= 0:
        raise ValueError("block_rows must be a positive integer")
    key_streams = [iter(sorted(columns[a], key=repr)) for a in range(dim)]
    merged = heapq.merge(*key_streams, key=repr)
    block = np.empty((block_rows, dim), dtype=np.complex128)
    used = 0
    previous = None
    have_previous = False
    for key in merged:
        if have_previous and key == previous:
            continue
        block[used, :] = [columns[a].get(key, 0j) for a in range(dim)]
        used += 1
        previous = key
        have_previous = True
        if used == block_rows:
            yield block.copy()
            used = 0
    if used:
        yield block[:used].copy()
    elif not have_previous:
        # An empty sparse support is the zero operator. Represent it by one
        # explicit zero row so the SVD routine can return its full nullspace.
        yield np.zeros((1, dim), dtype=np.complex128)


def direct_svd_from_columns(columns, dim, singular_rtol=1e-10, block_rows=4096):
    result = STREAM_SVD.streaming_svd(
        iter_operator_row_blocks(columns, dim, block_rows), dim, singular_rtol
    )
    return {
        "rank": result.rank,
        "singular_value_rtol": result.rank_rtol,
        "singular_values": [float(x) for x in result.singular_values],
        "right_singular_vectors_real": result.right_singular_vectors.real.tolist(),
        "right_singular_vectors_imag": result.right_singular_vectors.imag.tolist(),
        "right_singular_vectors_encoding": "separate real and imaginary arrays; columns are right singular vectors",
        "row_count": result.row_count,
        "block_count": result.block_count,
        "input_scale": result.input_scale,
        "max_abs_component": result.max_abs_component,
        "scale_separation_loss": result.scale_separation_loss,
        "method": result.method,
        "exact_rank_certified": False,
        "classification": (
            "INCONCLUSIVE_SCALE_SEPARATION_UNDERFLOW"
            if result.scale_separation_loss
            else "NUMERICAL_ONLY_NO_UPSTREAM_ERROR_BOUND"
        ),
    }



def kernel_projector_from_direct_svd(direct, dim):
    """Build a kernel projector only from a self-consistent full-V SVD record."""
    if not isinstance(dim, int) or isinstance(dim, bool) or dim <= 0:
        raise ValueError("dim must be a positive integer")
    rank = direct.get("rank")
    if not isinstance(rank, int) or isinstance(rank, bool) or rank < 0 or rank > dim:
        raise ValueError("direct-SVD rank is outside the input dimension")
    singular = np.asarray(direct.get("singular_values"), dtype=float)
    if (singular.shape != (dim,) or not np.all(np.isfinite(singular))
            or np.any(singular < 0.0)):
        raise ValueError("direct singular values must be a finite nonnegative dim-vector")
    if np.any(singular[:-1] < singular[1:]):
        raise ValueError("direct singular values must be sorted descending")
    rtol = direct.get("singular_value_rtol")
    if (isinstance(rtol, (bool, np.bool_)) or not isinstance(rtol, (int, float, np.integer, np.floating))
            or not math.isfinite(float(rtol)) or float(rtol) < 0.0):
        raise ValueError("direct SVD record must contain a finite nonnegative singular_value_rtol")
    sigma_max = float(singular[0])
    expected_rank = (
        int(np.count_nonzero(singular / sigma_max > float(rtol)))
        if sigma_max > 0.0 else 0
    )
    if rank != expected_rank:
        raise ValueError(
            "direct-SVD rank is inconsistent with singular values and singular_value_rtol"
        )
    real = np.asarray(direct.get("right_singular_vectors_real"), dtype=float)
    imag = np.asarray(direct.get("right_singular_vectors_imag"), dtype=float)
    if real.shape != (dim, dim) or imag.shape != (dim, dim):
        raise ValueError("full right-singular-vector matrices must be dim x dim")
    vectors = real + 1j * imag
    if not np.all(np.isfinite(vectors)):
        raise ValueError("right singular vectors contain non-finite entries")
    unitarity_defect = float(np.linalg.norm(
        vectors.conj().T @ vectors - np.eye(dim), ord="fro"
    ))
    unitarity_tolerance = 1e-8 * max(1.0, math.sqrt(dim))
    if (not math.isfinite(unitarity_defect)
            or unitarity_defect > unitarity_tolerance):
        raise ValueError(
            "right singular vectors are not orthonormal: "
            f"defect={unitarity_defect}, tolerance={unitarity_tolerance}"
        )
    kernel_basis = vectors[:, rank:]
    nullity = dim - rank
    orthogonality = (
        float(np.linalg.norm(kernel_basis.conj().T @ kernel_basis
                             - np.eye(nullity), ord="fro"))
        if nullity else 0.0
    )
    projector = kernel_basis @ kernel_basis.conj().T
    return {
        "rank": rank, "nullity": nullity, "kernel_basis": kernel_basis,
        "projector": projector,
        "right_singular_vector_unitarity_defect_frobenius": unitarity_defect,
        "kernel_basis_orthogonality_defect_frobenius": orthogonality,
        "projector_hermiticity_defect_frobenius":
            float(np.linalg.norm(projector - projector.conj().T, ord="fro")),
        "projector_idempotence_defect_frobenius":
            float(np.linalg.norm(projector @ projector - projector, ord="fro")),
    }


def compare_kernel_projectors(legacy_direct, no_prune_direct, dim):
    """Compare finite right-kernel subspaces, not merely their dimensions."""
    scale_flags = (
        legacy_direct.get("scale_separation_loss", True),
        no_prune_direct.get("scale_separation_loss", True),
    )
    if any(not isinstance(flag, (bool, np.bool_)) for flag in scale_flags):
        raise ValueError("scale_separation_loss metadata must be boolean")
    if any(scale_flags):
        # Physical-scale singular values may have underflowed after the
        # normalized SVD. Do not construct a kernel from that lossy spectrum.
        return {
            "status": "INCONCLUSIVE_SCALE_SEPARATION",
            "legacy_rank": None, "no_prune_rank": None,
            "legacy_nullity": None, "no_prune_nullity": None,
            "same_nullity": None,
            "projector_difference_spectral_norm": None,
            "max_principal_angle_sine": None,
            "scale_separation_warning": True,
            "projector_geometry_checks_passed": False,
            "interpretation": (
                "Scale separation may have erased nonzero singular modes; "
                "kernel ranks and projectors were not constructed."
            ),
        }
    left = kernel_projector_from_direct_svd(legacy_direct, dim)
    right = kernel_projector_from_direct_svd(no_prune_direct, dim)
    same_nullity = left["nullity"] == right["nullity"]
    projector_distance = float(np.linalg.norm(
        left["projector"] - right["projector"], ord=2
    ))
    scale_warning = bool(
        legacy_direct.get("scale_separation_loss", True)
        or no_prune_direct.get("scale_separation_loss", True)
    )
    geometry_valid = all(
        item[key] <= 1e-8
        for item in (left, right)
        for key in (
            "kernel_basis_orthogonality_defect_frobenius",
            "projector_hermiticity_defect_frobenius",
            "projector_idempotence_defect_frobenius",
        )
    )
    if not same_nullity:
        status, angle = "NULLITY_MISMATCH", None
    elif left["nullity"] == 0:
        status, angle = "TRIVIAL_ZERO_KERNEL", 0.0
    else:
        status = "SAME_DIMENSION_SUBSPACE_DISTANCE"
        angle = min(1.0, max(0.0, projector_distance))
    if scale_warning:
        status, angle = "INCONCLUSIVE_SCALE_SEPARATION", None
    elif not geometry_valid:
        status, angle = "INVALID_PROJECTOR_GEOMETRY", None
    return {
        "status": status,
        "legacy_rank": left["rank"], "no_prune_rank": right["rank"],
        "legacy_nullity": left["nullity"], "no_prune_nullity": right["nullity"],
        "same_nullity": same_nullity,
        "projector_difference_spectral_norm": projector_distance,
        "max_principal_angle_sine": angle,
        "scale_separation_warning": scale_warning,
        "projector_geometry_checks_passed": geometry_valid,
        "legacy_kernel_basis_orthogonality_defect_frobenius":
            left["kernel_basis_orthogonality_defect_frobenius"],
        "no_prune_kernel_basis_orthogonality_defect_frobenius":
            right["kernel_basis_orthogonality_defect_frobenius"],
        "legacy_projector_idempotence_defect_frobenius":
            left["projector_idempotence_defect_frobenius"],
        "no_prune_projector_idempotence_defect_frobenius":
            right["projector_idempotence_defect_frobenius"],
        "interpretation": (
            "For equal nonzero nullities, ||P_L-P_N||_2 equals sin of the "
            "largest principal angle. For unequal nullities, no principal-angle "
            "comparison is reported. Scale-separation warnings block certification."
        ),
    }


def singular_threshold_sweep(direct, dim, rtol_values):
    """Audit rank/kernel sensitivity across a declared relative-SVD threshold grid."""
    if not isinstance(dim, int) or isinstance(dim, bool) or dim <= 0:
        raise ValueError("dim must be a positive integer")
    try:
        raw_values = list(rtol_values)
    except TypeError as exc:
        raise ValueError("threshold grid must be an iterable of finite nonnegative values") from exc
    if any(isinstance(x, (bool, np.bool_)) for x in raw_values):
        raise ValueError("boolean values are not valid singular-value thresholds")
    try:
        values = [float(x) for x in raw_values]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("threshold grid must contain finite numeric values") from exc
    if not values or any(not math.isfinite(x) or x < 0.0 for x in values):
        raise ValueError("threshold grid must be nonempty, finite, and nonnegative")
    if values != sorted(set(values)):
        raise ValueError("threshold grid must be strictly increasing")
    scale_flag = direct.get("scale_separation_loss", True)
    if not isinstance(scale_flag, (bool, np.bool_)):
        raise ValueError("scale_separation_loss metadata must be boolean")
    if scale_flag:
        return {
            "spectral_scale": None,
            "right_singular_vector_unitarity_defect_frobenius": None,
            "right_singular_vector_unitarity_tolerance": None,
            "scale_separation_warning": True,
            "status": "INCONCLUSIVE_SCALE_SEPARATION",
            "rank_plateaus_on_declared_grid": [],
            "threshold_results": [],
            "thresholds_not_evaluated": values,
        }
    singular = np.asarray(direct.get("singular_values"), dtype=float)
    if singular.shape != (dim,) or not np.all(np.isfinite(singular)) or np.any(singular < 0.0):
        raise ValueError("direct singular values must be a finite nonnegative dim-vector")
    if np.any(singular[:-1] < singular[1:]):
        raise ValueError("direct singular values must be sorted descending")
    real = np.asarray(direct.get("right_singular_vectors_real"), dtype=float)
    imag = np.asarray(direct.get("right_singular_vectors_imag"), dtype=float)
    if real.shape != (dim, dim) or imag.shape != (dim, dim):
        raise ValueError("full right-singular-vector matrices must be dim x dim")
    vectors = real + 1j * imag
    if not np.all(np.isfinite(vectors)):
        raise ValueError("right singular vectors contain non-finite entries")
    # Kernel projectors and principal-angle identities require a unitary
    # full right-singular-vector matrix, not merely finite entries and shape.
    unitarity_defect = float(np.linalg.norm(
        vectors.conj().T @ vectors - np.eye(dim), ord="fro"
    ))
    unitarity_tolerance = 1e-8 * max(1.0, math.sqrt(dim))
    if (not math.isfinite(unitarity_defect)
            or unitarity_defect > unitarity_tolerance):
        raise ValueError(
            "right singular vectors are not orthonormal: "
            f"defect={unitarity_defect}, tolerance={unitarity_tolerance}"
        )
    sigma_max = float(singular[0]) if dim else 0.0
    normalized = singular / sigma_max if sigma_max > 0.0 else np.zeros(dim)
    rows = []
    for rtol in values:
        rank = int(np.count_nonzero(normalized > rtol)) if sigma_max > 0.0 else 0
        basis = vectors[:, rank:]
        projector = basis @ basis.conj().T
        rows.append({
            "singular_value_rtol": rtol,
            "rank": rank,
            "nullity": dim - rank,
            "largest_retained_normalized_sigma":
                float(normalized[rank - 1]) if rank else None,
            "largest_discarded_normalized_sigma":
                float(normalized[rank]) if rank < dim else None,
            "threshold_gap_below":
                float(normalized[rank - 1] - rtol) if rank else None,
            "threshold_gap_above":
                float(rtol - normalized[rank]) if rank < dim else None,
            "kernel_projector": projector,
        })
    plateaus = []
    for row in rows:
        if not plateaus or plateaus[-1]["rank"] != row["rank"]:
            plateaus.append({
                "rank": row["rank"], "nullity": row["nullity"],
                "rtol_min": row["singular_value_rtol"],
                "rtol_max": row["singular_value_rtol"],
                "grid_points": 1,
            })
        else:
            plateaus[-1]["rtol_max"] = row["singular_value_rtol"]
            plateaus[-1]["grid_points"] += 1
    return {
        "spectral_scale": sigma_max,
        "right_singular_vector_unitarity_defect_frobenius": unitarity_defect,
        "right_singular_vector_unitarity_tolerance": unitarity_tolerance,
        "scale_separation_warning": bool(direct.get("scale_separation_loss", True)),
        "rank_plateaus_on_declared_grid": plateaus,
        "threshold_results": rows,
    }


def compare_threshold_sweeps(legacy_direct, no_prune_direct, dim, rtol_values):
    """Compare ranks and kernel geometry at each identical threshold."""
    left = singular_threshold_sweep(legacy_direct, dim, rtol_values)
    right = singular_threshold_sweep(no_prune_direct, dim, rtol_values)
    scale_warning = (
        left["scale_separation_warning"] or right["scale_separation_warning"]
    )
    if scale_warning:
        return {
            "threshold_semantics": "relative threshold on singular values of C, not eigenvalues of C^dagger C",
            "legacy_rank_plateaus": [],
            "no_prune_rank_plateaus": [],
            "comparisons": [],
            "scale_separation_warning": True,
            "status": "INCONCLUSIVE_SCALE_SEPARATION",
            "interpretation": (
                "At least one physical-scale singular spectrum may have lost "
                "nonzero modes; no threshold rank or kernel geometry is reported."
            ),
        }
    comparisons = []
    for a, b in zip(left["threshold_results"], right["threshold_results"]):
        same_nullity = a["nullity"] == b["nullity"]
        distance = float(np.linalg.norm(
            a["kernel_projector"] - b["kernel_projector"], ord=2
        ))
        if not same_nullity:
            status, angle = "NULLITY_MISMATCH", None
        elif a["nullity"] == 0:
            status, angle = "TRIVIAL_ZERO_KERNEL", 0.0
        else:
            status = "SAME_DIMENSION_SUBSPACE_DISTANCE"
            angle = min(1.0, max(0.0, distance))
        if scale_warning:
            status, angle = "INCONCLUSIVE_SCALE_SEPARATION", None
        comparisons.append({
            "singular_value_rtol": a["singular_value_rtol"],
            "legacy_rank": a["rank"], "no_prune_rank": b["rank"],
            "legacy_nullity": a["nullity"], "no_prune_nullity": b["nullity"],
            "same_nullity": same_nullity,
            "projector_difference_spectral_norm": distance,
            "max_principal_angle_sine": angle,
            "status": status,
        })
    return {
        "threshold_semantics": "relative threshold on singular values of C, not eigenvalues of C^dagger C",
        "legacy_rank_plateaus": left["rank_plateaus_on_declared_grid"],
        "no_prune_rank_plateaus": right["rank_plateaus_on_declared_grid"],
        "comparisons": comparisons,
        "scale_separation_warning": scale_warning,
        "interpretation": (
            "Grid-based sensitivity diagnostic only. A rank plateau over sampled "
            "thresholds is not a proof of exact rank; thresholds between grid points "
            "and upstream operator-construction error remain unbounded."
        ),
    }

def column_digest(columns):
    records = []
    for a in sorted(columns):
        for (v, key), value in sorted(columns[a].items(), key=lambda item: repr(item[0])):
            z = complex(value)
            records.append((a, v, repr(key), z.real.hex(), z.imag.hex()))
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def run(domain="low", jmax2=5, rank_rtol=1e-10,
        singular_rtol=1e-10, block_rows=4096):
    if not isinstance(jmax2, int) or isinstance(jmax2, bool) or jmax2 < 1:
        raise ValueError("jmax2 must be a positive doubled-spin cutoff covering the j=1/2 carrier")
    if (isinstance(rank_rtol, (bool, np.bool_))
            or not isinstance(rank_rtol, (int, float, np.integer, np.floating))
            or not math.isfinite(float(rank_rtol)) or float(rank_rtol) < 0.0):
        raise ValueError("rank_rtol must be finite and nonnegative")
    if (isinstance(singular_rtol, (bool, np.bool_))
            or not isinstance(singular_rtol, (int, float, np.integer, np.floating))
            or not math.isfinite(float(singular_rtol)) or float(singular_rtol) < 0.0):
        raise ValueError("singular_rtol must be finite and nonnegative")
    if not isinstance(block_rows, int) or isinstance(block_rows, bool) or block_rows <= 0:
        raise ValueError("block_rows must be a positive integer")
    nodes = tuple(PILOT.PW.VERT)
    basis = PILOT.carrier_basis(domain)
    dim = len(basis)

    legacy_raw, legacy_supports = PILOT.constraint_columns(
        basis, nodes, jmax2, mode="legacy"
    )
    legacy, legacy_node_supports = stack_columns(legacy_raw, nodes, dim)
    legacy_hash = column_digest(legacy)
    legacy_gram = gram_from_columns(legacy, dim)
    del legacy_raw, legacy_node_supports

    # The producer cache is keyed by pruning mode; clear it between modes so
    # the second calculation cannot accidentally reuse legacy cache entries.
    if hasattr(PILOT.PW.T_cached, "cache_clear"):
        PILOT.PW.T_cached.cache_clear()
    gc.collect()

    no_prune_raw, no_prune_supports = PILOT.constraint_columns(
        basis, nodes, jmax2, mode="no_prune"
    )
    no_prune, no_prune_node_supports = stack_columns(no_prune_raw, nodes, dim)
    no_prune_hash = column_digest(no_prune)
    no_prune_gram = gram_from_columns(no_prune, dim)
    del no_prune_raw, no_prune_node_supports

    delta = subtract_stacked(legacy, no_prune, dim)
    delta_hash = column_digest(delta)
    delta_gram = gram_from_columns(delta, dim)

    # Direct singular values of C, not square roots of Gram eigenvalues.
    legacy_direct = direct_svd_from_columns(legacy, dim, singular_rtol, block_rows)
    no_prune_direct = direct_svd_from_columns(no_prune, dim, singular_rtol, block_rows)
    delta_direct = direct_svd_from_columns(delta, dim, singular_rtol, block_rows)
    kernel_comparison = compare_kernel_projectors(
        legacy_direct, no_prune_direct, dim
    )
    threshold_grid = (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4)
    threshold_sensitivity = compare_threshold_sweeps(
        legacy_direct, no_prune_direct, dim, threshold_grid
    )

    legacy_spec = spectrum_of_gram(legacy_gram, rank_rtol)
    no_prune_spec = spectrum_of_gram(no_prune_gram, rank_rtol)
    delta_spec = spectrum_of_gram(delta_gram, rank_rtol)
    sigma_min = legacy_direct["singular_values"][-1]
    delta_norm = delta_direct["singular_values"][0]
    scale_loss = bool(
        legacy_direct["scale_separation_loss"]
        or no_prune_direct["scale_separation_loss"]
        or delta_direct["scale_separation_loss"]
    )
    stable = bool(
        not scale_loss
        and legacy_direct["rank"] == dim
        and no_prune_direct["rank"] == dim
        and delta_norm < sigma_min
    )
    if scale_loss:
        sigma_min = None
        delta_norm = None
    return {
        "status": "legacy versus no-prune finite Peter-Weyl constraint-map comparison",
        "domain": domain,
        "domain_dimension": dim,
        "nodes": list(nodes),
        "Jmax": jmax2 / 2,
        "gram_eigenvalue_rtol": rank_rtol,
        "singular_value_rtol": singular_rtol,
        "streaming_block_rows": block_rows,
        "row_space_definition": "direct sum over node labels; row keys are (node, output Peter-Weyl basis key)",
        "legacy": {**legacy_spec, "direct_svd": legacy_direct,
                   "stacked_column_sha256": legacy_hash,
                   "node_column_supports": legacy_supports},
        "no_prune": {**no_prune_spec, "direct_svd": no_prune_direct,
                     "stacked_column_sha256": no_prune_hash,
                     "node_column_supports": no_prune_supports},
        "difference_no_prune_minus_legacy": {
            **delta_spec, "direct_svd": delta_direct,
            "stacked_column_sha256": delta_hash,
            "interpretation": "operator norm of C_no_prune-C_legacy",
        },
        "physical_kernel_subspace_comparison": kernel_comparison,
        "spectral_threshold_sensitivity": threshold_sensitivity,
        "weyl_rank_stability_sufficient_condition": {
            "delta_operator_norm": delta_norm,
            "legacy_minimum_singular_value": sigma_min,
            "strict_inequality_holds": (
                None if scale_loss else bool(delta_norm < sigma_min)
            ),
            "legacy_full_column_rank_at_declared_tolerance": (
                None if scale_loss else bool(legacy_direct["rank"] == dim)
            ),
            "no_prune_full_column_rank_at_declared_tolerance": (
                None if scale_loss else bool(no_prune_direct["rank"] == dim)
            ),
            "scale_separation_warning": scale_loss,
            "rank_source": "direct streaming SVD of C",
            "sufficient_rank_stability_for_computed_float_maps": stable,
            "exact_operator_rank_certified": False,
            "scope": "represented float64 maps only; upstream operator-construction and subtraction error bounds are not included",
            "note": "Weyl sufficient condition on the represented maps only; failure to satisfy it does not prove rank loss.",
        },
        "claim_boundary": "Finite carrier and finite Jmax only; float64 numerical result, not exact arithmetic or a continuum theorem.",
    }


def cutoff_convergence_sweep(domain="low", jmax2_values=(1, 3, 5),
                               rank_rtol=1e-10, singular_rtol=1e-10,
                               block_rows=4096):
    """Compare finite-cutoff kernel geometry and normalized spectra.

    This is a regulator-sensitivity diagnostic. It does not estimate the
    operator-norm difference between C_J and C_{J+1}, because this runner
    retains spectra/projectors rather than the complete sparse maps.
    """
    try:
        cutoffs = list(jmax2_values)
    except TypeError as exc:
        raise ValueError("jmax2_values must be an iterable of increasing integers") from exc
    if (not cutoffs or any(
            not isinstance(x, int) or isinstance(x, bool) or x < 1
            for x in cutoffs)):
        raise ValueError("every doubled-spin cutoff must be an integer >= 1")
    if cutoffs != sorted(set(cutoffs)):
        raise ValueError("doubled-spin cutoffs must be strictly increasing")

    results = []
    summaries = []
    for cutoff in cutoffs:
        result = run(domain, cutoff, rank_rtol, singular_rtol, block_rows)
        results.append(result)
        np_direct = result["no_prune"]["direct_svd"]
        leg_direct = result["legacy"]["direct_svd"]

        def spectrum_summary(direct):
            if direct["scale_separation_loss"]:
                return {
                    "rank": None, "nullity": None, "sigma_max": None,
                    "sigma_min_including_numerical_kernel": None,
                    "smallest_retained_sigma": None,
                    "normalized_singular_values": None,
                    "scale_separation_warning": True,
                    "stacked_column_sha256": None,
                    "classification": "INCONCLUSIVE_SCALE_SEPARATION",
                }
            singular = np.asarray(direct["singular_values"], dtype=float)
            scale = float(singular[0]) if singular.size else 0.0
            normalized = singular / scale if scale > 0.0 else np.zeros_like(singular)
            rank = int(direct["rank"])
            positive = singular[:rank]
            return {
                "rank": rank,
                "nullity": int(len(singular) - rank),
                "sigma_max": scale,
                "sigma_min_including_numerical_kernel": (
                    float(singular[-1]) if singular.size else 0.0
                ),
                "smallest_retained_sigma": (
                    float(positive[-1]) if positive.size else None
                ),
                "normalized_singular_values": normalized.tolist(),
                "scale_separation_warning": bool(direct["scale_separation_loss"]),
                "stacked_column_sha256": None,
            }

        no_prune_summary = spectrum_summary(np_direct)
        legacy_summary = spectrum_summary(leg_direct)
        no_prune_summary["stacked_column_sha256"] = result["no_prune"]["stacked_column_sha256"]
        legacy_summary["stacked_column_sha256"] = result["legacy"]["stacked_column_sha256"]
        summaries.append({
            "jmax2": cutoff,
            "Jmax": cutoff / 2.0,
            "domain_dimension": result["domain_dimension"],
            "legacy": legacy_summary,
            "no_prune": no_prune_summary,
        })

        if hasattr(PILOT.PW.T_cached, "cache_clear"):
            PILOT.PW.T_cached.cache_clear()
        gc.collect()

    adjacent = []
    for index in range(1, len(results)):
        previous = results[index - 1]
        current = results[index]
        previous_direct = previous["no_prune"]["direct_svd"]
        current_direct = current["no_prune"]["direct_svd"]
        if previous_direct["scale_separation_loss"] or current_direct["scale_separation_loss"]:
            adjacent.append({
                "from_jmax2": cutoffs[index - 1],
                "to_jmax2": cutoffs[index],
                "no_prune_kernel_subspace": {
                    "status": "INCONCLUSIVE_SCALE_SEPARATION",
                    "max_principal_angle_sine": None,
                },
                "max_abs_difference_of_normalized_singular_values": None,
                "normalized_spectrum_comparison_scope": (
                    "not evaluated because a physical-scale singular spectrum "
                    "may have lost nonzero modes"
                ),
            })
            continue
        previous_sigma = np.asarray(previous_direct["singular_values"], dtype=float)
        current_sigma = np.asarray(current_direct["singular_values"], dtype=float)
        previous_scale = float(previous_sigma[0]) if previous_sigma.size else 0.0
        current_scale = float(current_sigma[0]) if current_sigma.size else 0.0
        previous_normalized = (
            previous_sigma / previous_scale if previous_scale > 0.0
            else np.zeros_like(previous_sigma)
        )
        current_normalized = (
            current_sigma / current_scale if current_scale > 0.0
            else np.zeros_like(current_sigma)
        )
        kernel = compare_kernel_projectors(
            previous_direct, current_direct, previous["domain_dimension"]
        )
        adjacent.append({
            "from_jmax2": cutoffs[index - 1],
            "to_jmax2": cutoffs[index],
            "no_prune_kernel_subspace": kernel,
            "max_abs_difference_of_normalized_singular_values": float(
                np.max(np.abs(current_normalized - previous_normalized))
            ),
            "normalized_spectrum_comparison_scope": (
                "shape-only diagnostic after independent sigma_max normalization; "
                "not an operator-norm convergence bound"
            ),
        })

    return {
        "status": "finite Peter-Weyl regulator sensitivity sweep",
        "domain": domain,
        "domain_dimension": summaries[0]["domain_dimension"],
        "doubled_spin_cutoffs": cutoffs,
        "points": summaries,
        "adjacent_cutoff_comparisons": adjacent,
        "claim_boundary": (
            "Compares finite-cutoff numerical ranks, normalized singular spectra, "
            "and right-kernel projectors only. It does not compute ||C_J-C_Jnext||, "
            "prove nested-operator convergence, establish regulator independence, "
            "or imply a continuum theorem."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", choices=("low", "full"), default="low")
    ap.add_argument("--jmax2", type=int, default=5)
    ap.add_argument("--rank-rtol", type=float, default=1e-10,
                    help="relative Gram-eigenvalue tolerance (not singular-value tolerance)")
    ap.add_argument("--singular-rtol", type=float, default=1e-10,
                    help="relative singular-value tolerance for direct TSQR/SVD")
    ap.add_argument("--block-rows", type=int, default=4096)
    ap.add_argument(
        "--jmax2-sweep",
        help="optional comma-separated increasing doubled-spin cutoffs, e.g. 1,3,5",
    )
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.jmax2_sweep:
        try:
            cutoffs = [int(value.strip()) for value in args.jmax2_sweep.split(",")]
        except ValueError as exc:
            ap.error(f"invalid --jmax2-sweep: {exc}")
        result = cutoff_convergence_sweep(
            args.domain, cutoffs, args.rank_rtol,
            args.singular_rtol, args.block_rows
        )
        default_output = (
            f"verification_results/PETER_WEYL_CONSTRAINT_STABILITY_"
            f"{args.domain}_J_SWEEP.json"
        )
    else:
        result = run(args.domain, args.jmax2, args.rank_rtol,
                     args.singular_rtol, args.block_rows)
        default_output = (
            f"verification_results/PETER_WEYL_CONSTRAINT_STABILITY_"
            f"{args.domain}_J{args.jmax2}.json"
        )
    output = args.output or Path(default_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
