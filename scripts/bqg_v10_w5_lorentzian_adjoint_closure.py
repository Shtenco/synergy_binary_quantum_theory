#!/usr/bin/env python3
"""Build a W5-seeded, adjoint-closed finite basis for the raw Lorentzian KKV map.

The script computes BOTH the forward epsilon K-K-V map and its independently
ordered adjoint. It closes the reachable Gauss basis under both maps before
comparing the direct adjoint columns with the conjugate transpose of the
assembled forward matrix.

This is a structural adjoint-consistency gate only. It does not choose the
physical Hermitian ordering, restore beta/kappa/hbar factors, or claim a
physical Hamiltonian/kernel result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import deque
from pathlib import Path
import sys


W5_K = (0, 2, 4, 6, 8)
NODES = (0, 1, 2, 3, 4)
from bqg_sparse_linear_accumulator import SparseLinearAccumulator

def build_pruning_coverage(static_inventory):
    """Separate static site classification from runtime zero-loss evidence."""
    inventory_complete = bool(static_inventory.get("complete", False))
    fixed_mask_classes = {
        "post_assembly_fixed_operator_mask_requires_error_bound",
        "post_assembly_fixed_operator_mask_with_frobenius_spectral_error_bound",
    }
    fixed_masks = [
        site for site in static_inventory.get("threshold_sites", [])
        if site.get("classification") in fixed_mask_classes
    ]
    certified_bound_api_available = any(
        site.get("classification") ==
        "post_assembly_fixed_operator_mask_with_frobenius_spectral_error_bound"
        for site in fixed_masks
    )
    return {
        # Overall coverage intentionally remains false until runtime loss,
        # cache-equivalence, support-cap and convergence controls are closed.
        "complete": False,
        "status": (
            "STATIC_THRESHOLD_INVENTORY_PASS_RUNTIME_LOSS_CONTROLS_PENDING"
            if inventory_complete else
            "STATIC_THRESHOLD_INVENTORY_INCOMPLETE"
        ),
        "static_threshold_inventory_complete": inventory_complete,
        "static_inventory": {
            "module_count": int(static_inventory.get("module_count", 0)),
            "threshold_site_count": int(static_inventory.get("threshold_site_count", 0)),
            "unclassified_sites": list(static_inventory.get("unclassified_sites", [])),
            "plain_lru_caches": list(static_inventory.get("plain_lru_caches", [])),
            "audited_lru_caches": list(static_inventory.get("audited_lru_caches", [])),
            "post_assembly_fixed_operator_masks": fixed_masks,
        },
        "runtime_loss_controls_complete": False,
        "fixed_operator_mask_bound_api_available": certified_bound_api_available,
        "fixed_operator_mask_bound_applied_in_production": False,
        "fixed_operator_mask_error_bound_complete": False,
        "remaining_checks": [
            "wire the certified matrix-mask API into any production call site that uses this mask and persist its error ledger",
            "demonstrate that the post-assembly mask is absent from or explicitly bounded in the actual W5 execution path",
            "cache-on/cache-off pruning-ledger equality for the full composition chain",
            "complete spin-wall and state/depth-cap loss separation",
            "full complex support convergence across Jmax2 and amplitude-threshold ladders",
        ],
        "excluded_from_amplitude_bypass": [
            "exact SU(2) representation and parity selection rules",
            "spin-wall truncation and state/depth caps",
            "zero-aware spectral volume/nullspace convention",
            "post-assembly fixed operator matrix masks",
        ],
        "shadow_no_prune_interpretation": (
            "shadow_no_prune bypasses only wired nonzero amplitude thresholds. "
            "A passing static inventory proves site classification, not zero runtime loss. "
            "Post-assembly fixed matrix masks and support caps are not removed by this bypass. "
            "The shadow result remains partial until runtime loss, cache-equivalence, "
            "cap-accounting, fixed-mask error-bound and support-convergence controls pass."
        ),
    }


def epsilon_triples(neighbors):
    import itertools
    out = []
    for r, omitted in enumerate(neighbors):
        tri = tuple(x for x in neighbors if x != omitted)
        for p in itertools.permutations(tri):
            positions = [tri.index(x) for x in p]
            inversions = sum(
                positions[i] > positions[j]
                for i in range(3) for j in range(i + 1, 3)
            )
            sign = (-1) ** r * (-1 if inversions % 2 else 1)
            out.append((sign, p))
    if len(out) != 24 or len({p for _, p in out}) != 24:
        raise AssertionError("the K5 node epsilon orbit must contain 24 unique triples")
    return out


def w5_input_key(K):
    if K not in W5_K:
        raise ValueError(f"W5 K must be one of {W5_K}")
    return ((4,) * 10, (K, 0, 0, 0, 0))


def sum_sparse_maps(states):
    """Coherently sum complete node maps before any state-level thresholding."""
    accumulator = SparseLinearAccumulator()
    for state in states:
        accumulator.add(state)
    return accumulator.finish()


def sparse_norm(state):
    """Stable Euclidean norm for sparse complex amplitudes."""
    norm = 0.0
    for amplitude in state.values():
        norm = math.hypot(norm, abs(amplitude))
    return norm


def prune_scalar_channel(state):
    # The final rank-(0+1) covariant leg can reach total source J=0 only
    # from J=0 or J=1. Apply this exact selection rule after the middle leg.
    return {key: amp for key, amp in state.items() if key[2] in (0, 2)}


def triple_forward(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
    """Tr_aux[C_a(K) C_b(K) C_c(V)] acting right-to-left."""
    out = SparseLinearAccumulator()
    max_leak = 0.0
    max_diag = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                q, leak = cv_component(psi_cov, source_v, c, k, i, Jmax2)
                max_leak = max(max_leak, float(leak))
                if not q:
                    continue
                q, d1 = ck_component(q, source_v, b, j, k, Jmax2)
                _merge_diag(max_diag, d1)
                if not q:
                    continue
                q = prune_scalar_channel(q)
                if not q:
                    continue
                q, d2 = ck_component(q, source_v, a, i, j, Jmax2)
                _merge_diag(max_diag, d2)
                out.add(q)
    return out.finish(), max_leak, max_diag


def triple_adjoint(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
    """Directly apply T_abc^dagger = Tr_aux[C_c(V) C_b(K) C_a(K)].

    K_sine is anti-Hermitian and V is Hermitian. The two K-adjoint minus
    signs cancel. The auxiliary indices are transposed explicitly:
      K_a[j,i] -> K_b[k,j] -> V_c[i,k]
    in right-to-left action order.
    """
    out = SparseLinearAccumulator()
    max_leak = 0.0
    max_diag = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                q, d1 = ck_component(psi_cov, source_v, a, j, i, Jmax2)
                _merge_diag(max_diag, d1)
                if not q:
                    continue
                q, d2 = ck_component(q, source_v, b, k, j, Jmax2)
                _merge_diag(max_diag, d2)
                if not q:
                    continue
                q = prune_scalar_channel(q)
                if not q:
                    continue
                q, leak = cv_component(q, source_v, c, i, k, Jmax2)
                max_leak = max(max_leak, float(leak))
                out.add(q)
    return out.finish(), max_leak, max_diag


def node_kkv_maps(gauss_state, source_v, Jmax2, PW, CV, RAW):
    """Compute forward B_v and direct reversed-order B_v^dagger on one input."""
    psi_cov = CV.gauss_to_covariant(gauss_state, source_v)
    triples = epsilon_triples(PW.NEIG[source_v])
    B_acc, Bdag_acc = SparseLinearAccumulator(), SparseLinearAccumulator()
    diagnostics = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    for epsilon, (a, b, c) in triples:
        forward, leak_f, diag_f = triple_forward(
            psi_cov, source_v, a, b, c, Jmax2,
            RAW.COMP.C_volume_component, RAW.KCOMP.C_K_component
        )
        adjoint, leak_a, diag_a = triple_adjoint(
            psi_cov, source_v, a, b, c, Jmax2,
            RAW.COMP.C_volume_component, RAW.KCOMP.C_K_component
        )
        B_acc.add(forward, epsilon)
        Bdag_acc.add(adjoint, epsilon)
        diagnostics["CV_complete_basis_leakage"] = max(
            diagnostics["CV_complete_basis_leakage"], leak_f, leak_a
        )
        _merge_diag(diagnostics, diag_f)
        _merge_diag(diagnostics, diag_a)
    diagnostics["CK_outer_complete_basis_leakage"] = max(
        float(diagnostics.get("outer_complete_basis_leakage", 0.0)),
        float(diagnostics.get("CK_outer_complete_basis_leakage", 0.0)),
    )
    diagnostics["CK_internal_volume_sector_leakage"] = max(
        float(diagnostics.get("internal_volume_sector_leakage", 0.0)),
        float(diagnostics.get("CK_internal_volume_sector_leakage", 0.0)),
    )
    diagnostics["CK_complete_charge_basis_leakage"] = max(
        float(diagnostics.get("complete_charge_basis_leakage", 0.0)),
        float(diagnostics.get("CK_complete_charge_basis_leakage", 0.0)),
    )
    return B_acc.finish(), Bdag_acc.finish(), diagnostics


def _merge_diag(dst, src):
    for name, value in src.items():
        if isinstance(value, (int, float)):
            dst[name] = max(dst.get(name, 0.0), float(value))


def covariant_scalar_to_gauss(key, source_v):
    spins, Kother, J2, M2, K12, K34 = key
    if J2 != 0 or M2 != 0 or K12 != K34:
        raise ValueError("epsilon-traced output is not a Gauss scalar")
    labels = list(Kother)
    if len(labels) != 5 or labels[source_v] != -1:
        raise ValueError("malformed covariant source sentinel")
    labels[source_v] = K12
    if any(x < 0 for x in labels):
        raise ValueError("unclosed non-source Gauss label")
    return (tuple(spins), tuple(labels))


def covariant_to_gauss_map(state, source_v):
    out, rejected, rejected_norm2 = {}, 0, 0.0
    for key, amp in state.items():
        try:
            gkey = covariant_scalar_to_gauss(key, source_v)
        except ValueError:
            rejected += 1
            rejected_norm2 += abs(amp) ** 2
            continue
        out[gkey] = out.get(gkey, 0j) + amp
    return out, rejected, math.sqrt(rejected_norm2)


def _matrix_adjoint_expected_columns(forward_columns, basis_keys, key_index):
    n = len(basis_keys)
    expected = [dict() for _ in range(n)]
    for col, state in enumerate(forward_columns):
        for row_key, amp in state.items():
            row = key_index.get(row_key)
            if row is None:
                continue
            expected[row][basis_keys[col]] = amp.conjugate()
    return expected


def sparse_matrix_residual(actual_columns, expected_columns):
    err_norm = 0.0
    scale_norm = 0.0
    for actual, expected in zip(actual_columns, expected_columns):
        for key in set(actual) | set(expected):
            a = actual.get(key, 0j)
            b = expected.get(key, 0j)
            err_norm = math.hypot(err_norm, abs(a - b))
            scale_norm = math.hypot(scale_norm, abs(b))
    return err_norm, scale_norm


def sparse_frobenius(columns):
    norm = 0.0
    for column in columns:
        for amplitude in column.values():
            norm = math.hypot(norm, abs(amplitude))
    return norm


def sparse_sum_difference_norm(columns_a, columns_b, sign=1.0):
    norm = 0.0
    for a, b in zip(columns_a, columns_b):
        for key in set(a) | set(b):
            norm = math.hypot(norm, abs(a.get(key, 0j) + sign * b.get(key, 0j)))
    return norm

def closure_acceptance(
    forward_norm,
    adjoint_norm,
    adjoint_relative_residual,
    scalar_projection_rejected_terms,
    max_component_leakage,
    max_spin_reached,
    sine_engine_result,
    spin_wall_omitted_branches=0,
    max_intermediate_spin_reached=0.0,
):
    """Conjunctive acceptance: algebraic closure AND the underlying LSINE gate.

    A small direct-adjoint residual cannot override a failed representation,
    leakage, spin-wall, or sine-Hermitian operator gate.
    """
    return bool(
        forward_norm > 1e-10
        and adjoint_norm > 1e-10
        and adjoint_relative_residual < 1e-8
        and scalar_projection_rejected_terms == 0
        and max_component_leakage < 1e-9
        and max_spin_reached <= 3.5 + 1e-12
        and spin_wall_omitted_branches == 0
        and max_intermediate_spin_reached <= 3.5 + 1e-12
        and isinstance(sine_engine_result, dict)
        and sine_engine_result.get("passed") is True
    )


def sparse_ordering_parts(forward_columns, adjoint_columns):
    """Return B_sym=(B+B†)/2, B_anti=(B-B†)/2 and B_H,sym at hbar=1.

    The final object is only the algebraic two-Hermitian-K conversion
    B_H,sym=-B_sym/hbar^2 evaluated at hbar=1. It is not the full physical
    Lorentzian constraint; its external regularization coefficient is open.
    """
    if len(forward_columns) != len(adjoint_columns):
        raise ValueError("forward and adjoint column counts differ")
    symmetric, antihermitian, hermitian_k_hbar1 = [], [], []
    for forward, adjoint in zip(forward_columns, adjoint_columns):
        sym, anti, kh = {}, {}, {}
        for key in set(forward) | set(adjoint):
            f = forward.get(key, 0j)
            d = adjoint.get(key, 0j)
            s = 0.5 * (f + d)
            a = 0.5 * (f - d)
            h = -s
            if s != 0:
                sym[key] = s
                kh[key] = h
            if a != 0:
                anti[key] = a
        symmetric.append(sym)
        antihermitian.append(anti)
        hermitian_k_hbar1.append(kh)
    return symmetric, antihermitian, hermitian_k_hbar1




def _factorized_branch_metrics(tensors, scalar_amp, include_payload=False):
    """Measure a factorized spin-network branch, not merely its scalar prefactor."""
    import hashlib
    import numpy as np

    digest = hashlib.sha256()
    branch_l1 = float(abs(scalar_amp))
    branch_l2 = float(abs(scalar_amp))
    branch_max = float(abs(scalar_amp))
    shapes = []
    payload = []
    for tensor in tensors:
        arr = np.ascontiguousarray(np.asarray(tensor))
        digest.update(repr(tuple(arr.shape)).encode("ascii"))
        digest.update(arr.dtype.str.encode("ascii"))
        digest.update(arr.tobytes(order="C"))
        absolute = np.abs(arr)
        local_l1 = float(absolute.sum())
        local_l2 = float(np.linalg.norm(absolute.reshape(-1)))
        local_max = float(absolute.max()) if absolute.size else 0.0
        branch_l1 *= local_l1
        branch_l2 *= local_l2
        branch_max *= local_max
        shapes.append(list(arr.shape))
        if include_payload:
            payload.append({
                "shape": list(arr.shape),
                "dtype": arr.dtype.str,
                "real": np.asarray(arr.real).reshape(-1).astype(float).tolist(),
                "imag": np.asarray(arr.imag).reshape(-1).astype(float).tolist(),
            })
    return {
        "branch_l1": branch_l1,
        "branch_l2": branch_l2,
        "branch_max_abs": branch_max,
        "tensor_sha256": digest.hexdigest(),
        "tensor_shapes": shapes,
        "tensor_payload": payload if include_payload else None,
    }


def _make_spin_wall_probe(original_apply_hit_branch, edge_index_map, wall, stats):
    """Wrap primitive hits and retain the exact nonzero branches rejected at the upper spin wall.

    The widened call is diagnostic only; the returned production branches are
    still generated with the original wall. Since one hit changes doubled spin
    by exactly +/-1, wall+2 is sufficient to expose every upper-wall candidate.
    """
    def probed(branch, x, y, i, j, Jmax2):
        edge = tuple(sorted((x, y)))
        edge_index = edge_index_map[edge]
        input_spins, _tensors, _input_amp = branch
        if Jmax2 == wall and input_spins[edge_index] == wall:
            widened = original_apply_hit_branch(branch, x, y, i, j, wall + 2)
            for candidate in widened:
                out_spins, _out_tensors, out_amp = candidate
                if out_spins[edge_index] > wall and abs(out_amp) > 0.0:
                    metrics = _factorized_branch_metrics(
                        _out_tensors, out_amp, include_payload=len(stats["samples"]) < 5
                    )
                    magnitude = metrics["branch_l2"]
                    stats["count"] += 1
                    stats["l1"] += metrics["branch_l1"]
                    stats["l2_sq"] += magnitude * magnitude
                    stats["max_abs"] = max(stats["max_abs"], metrics["branch_max_abs"])
                    edge_label = f"{edge[0]}-{edge[1]}"
                    stats["by_edge"][edge_label] = stats["by_edge"].get(edge_label, 0) + 1
                    if len(stats["samples"]) < 25:
                        stats["samples"].append({
                            "edge": list(edge),
                            "auxiliary_indices": [int(i), int(j)],
                            "input_spins": list(input_spins),
                            "omitted_output_spins": list(out_spins),
                            "scalar_factor_real": float(out_amp.real),
                            "scalar_factor_imag": float(out_amp.imag),
                            "branch_l1": metrics["branch_l1"],
                            "branch_l2": metrics["branch_l2"],
                            "branch_max_abs": metrics["branch_max_abs"],
                            "tensor_sha256": metrics["tensor_sha256"],
                            "tensor_shapes": metrics["tensor_shapes"],
                            "tensor_payload": metrics["tensor_payload"],
                        })
        actual = original_apply_hit_branch(branch, x, y, i, j, Jmax2)
        for candidate in actual:
            out_spins, _out_tensors, out_amp = candidate
            reached = max(out_spins, default=0) / 2
            stats["max_intermediate_spin_reached"] = max(
                stats["max_intermediate_spin_reached"], float(reached)
            )
            if max(out_spins, default=0) > 7 and abs(out_amp) > 0.0:
                metrics = _factorized_branch_metrics(
                    _out_tensors,
                    out_amp,
                    include_payload=len(stats["intermediate_samples"]) < 5,
                )
                magnitude = metrics["branch_l2"]
                stats["intermediate_above_3_5_count"] += 1
                stats["intermediate_above_3_5_l1"] += metrics["branch_l1"]
                stats["intermediate_above_3_5_l2_sq"] += magnitude * magnitude
                if len(stats["intermediate_samples"]) < 25:
                    stats["intermediate_samples"].append({
                        "edge": list(edge),
                        "auxiliary_indices": [int(i), int(j)],
                        "output_spins": list(out_spins),
                        "scalar_factor_real": float(out_amp.real),
                        "scalar_factor_imag": float(out_amp.imag),
                        "branch_l1": metrics["branch_l1"],
                        "branch_l2": metrics["branch_l2"],
                        "branch_max_abs": metrics["branch_max_abs"],
                        "tensor_sha256": metrics["tensor_sha256"],
                        "tensor_shapes": metrics["tensor_shapes"],
                        "tensor_payload": metrics["tensor_payload"],
                    })
        return actual
    return probed


def _spin_wall_diagnostics(stats):
    return {
        "max_intermediate_spin_reached": float(stats["max_intermediate_spin_reached"]),
        "intermediate_above_3_5_branch_count": int(stats["intermediate_above_3_5_count"]),
        "intermediate_above_3_5_pathwise_l1": float(stats["intermediate_above_3_5_l1"]),
        "intermediate_above_3_5_pathwise_l2": math.sqrt(float(stats["intermediate_above_3_5_l2_sq"])),
        "intermediate_above_3_5_samples": list(stats["intermediate_samples"]),
        "spin_wall_omitted_branch_count": int(stats["count"]),
        "spin_wall_omitted_pathwise_l1": float(stats["l1"]),
        "spin_wall_omitted_pathwise_l2": math.sqrt(float(stats["l2_sq"])),
        "spin_wall_omitted_max_abs": float(stats["max_abs"]),
        "spin_wall_omitted_by_edge": dict(stats["by_edge"]),
        "spin_wall_omitted_samples": list(stats["samples"]),
        "spin_wall_probe_note": (
            "Counts nonzero primitive path branches immediately above the configured "
            "upper doubled-spin wall before production rejection; L1/L2 are pathwise "
            "sums and are not coherently aggregated output norms."
        ),
    }

def _closure_cap_reasons(closure_size, input_depth, max_states, max_depth):
    """Return every active support cap independently; None means unbounded."""
    reasons = {}
    if max_states is not None and closure_size >= max_states:
        reasons["max_states"] = int(max_states)
    if max_depth is not None and input_depth >= max_depth:
        reasons["max_depth"] = int(max_depth)
    return reasons


def run_closure(source_nodes=NODES, Jmax2=10, max_states=500, max_depth=8, amp_tol=1e-12, pruning_mode="production", spectral_zero_factor=1000.0, cache_mode="on", seed_Ks=None):
    source_nodes = tuple(sorted(set(source_nodes)))
    if not source_nodes or any(v not in NODES for v in source_nodes):
        raise ValueError("source_nodes must be a nonempty subset of 0..4")
    selected_seed_Ks = tuple(W5_K if seed_Ks is None else sorted(set(seed_Ks)))
    if not selected_seed_Ks or any(K not in W5_K for K in selected_seed_Ks):
        raise ValueError(f"seed_Ks must be a nonempty subset of {W5_K}")
    if Jmax2 < 10:
        raise ValueError("W5 j=2 inputs require Jmax2>=10 as an initial support wall")
    if ((max_states is not None and max_states < len(selected_seed_Ks))
            or (max_depth is not None and max_depth < 0) or amp_tol < 0):
        raise ValueError("invalid closure limits")
    if pruning_mode not in ("production", "shadow_no_prune"):
        raise ValueError("pruning_mode must be production or shadow_no_prune")
    if not math.isfinite(float(spectral_zero_factor)) or float(spectral_zero_factor) <= 0:
        raise ValueError("spectral_zero_factor must be finite and positive")
    if cache_mode not in ("on", "off"):
        raise ValueError("cache_mode must be on or off")

    HERE = Path(__file__).resolve().parent
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import k5_peter_weyl_safe_hda_column as PW
    import peter_weyl_covariant_volume_leg_gate as CV
    import peter_weyl_lorentzian_sine_ordered_triple_gate as LSINE
    import peter_weyl_lorentzian_ordered_triple_gate as RAW
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM
    import bqg_pruning_audit as PRUNE
    import bqg_pruning_coverage_inventory as PRUNE_INVENTORY
    static_pruning_inventory = PRUNE_INVENTORY.run_inventory()
    ZVM.ZERO_FACTOR = float(spectral_zero_factor)

    # A lower-level pruning ledger must not depend on whether a prior caller
    # warmed an lru_cache before capture. Clear every imported operator cache
    # in this dependency family before the audited execution.
    cleared_cache_functions = 0
    for module_name, module in list(sys.modules.items()):
        module_file = getattr(module, "__file__", None)
        if module_file and Path(module_file).resolve().parent == HERE:
            for candidate in vars(module).values():
                clear = getattr(candidate, "cache_clear", None)
                if callable(clear):
                    clear()
                    cleared_cache_functions += 1

    seeds = [w5_input_key(K) for K in selected_seed_Ks]
    closure_keys = list(seeds)
    seen = set(closure_keys)
    depth = {key: 0 for key in seeds}
    queue = deque(seeds)
    forward_columns, adjoint_columns = [], []
    leakage_max = {}
    aggregate_leakage = {}
    rejected_count = 0
    rejected_norm2 = 0.0
    truncated = False
    omitted_support = 0

    # Full coefficient-level loss ledger.  Records are intentionally not
    # sample-capped: every coefficient removed by the V10 output threshold or
    # by a reachable-support cap is retained with its complex value and keys.
    coefficient_ledger = {
        "amp_tol": {"forward": [], "adjoint": []},
        "closure_caps": {"forward": [], "adjoint": []},
        "summary": {},
    }
    coherent_sum_ledger = []

    def node_contribution_rows(direction, source_v, input_key, state):
        rows = []
        for output_key, amplitude in state.items():
            z = complex(amplitude)
            rows.append({
                "source_node": int(source_v),
                "direction": direction,
                "input_key": repr(input_key),
                "output_key": repr(output_key),
                "real": float(z.real),
                "imag": float(z.imag),
                "abs": float(abs(z)),
            })
        return rows

    def coherent_metrics(rows, state):
        path_l1 = sum(row["abs"] for row in rows)
        path_l2 = math.sqrt(sum(row["abs"] ** 2 for row in rows))
        coherent_l1 = sum(abs(z) for z in state.values())
        coherent_l2 = sparse_norm(state)
        return {
            "pre_sum_branch_coefficient_count": len(rows),
            "pre_sum_pathwise_l1": float(path_l1),
            "pre_sum_pathwise_l2": float(path_l2),
            "post_sum_coherent_support": len(state),
            "post_sum_coherent_l1": float(coherent_l1),
            "post_sum_coherent_l2": float(coherent_l2),
            "l1_cancellation_gap": float(path_l1 - coherent_l1),
        }

    def record_dropped(stage, direction, input_key, output_key, amplitude,
                       threshold=None, depth_value=None, cap=None, cap_reasons=None):
        z = complex(amplitude)
        row = {
            "stage": stage,
            "direction": direction,
            "input_key": repr(input_key),
            "output_key": repr(output_key),
            "real": float(z.real),
            "imag": float(z.imag),
            "abs": float(abs(z)),
        }
        if threshold is not None:
            row["threshold"] = float(threshold)
        if depth_value is not None:
            row["input_depth"] = int(depth_value)
        if cap is not None:
            if isinstance(cap, dict):
                row["cap_limits"] = dict(cap)
            else:
                row["cap"] = int(cap)
        if cap_reasons:
            row["cap_reasons"] = dict(cap_reasons)
        row["actual_discarded"] = bool(not PRUNE.keep(z, threshold or 0.0)) if stage == "amp_tol" else True
        row["pruning_mode"] = pruning_mode
        coefficient_ledger[stage][direction].append(row)

    def summarize_ledger():
        summary = {}
        for stage in ("amp_tol", "closure_caps"):
            for direction in ("forward", "adjoint"):
                rows = coefficient_ledger[stage][direction]
                actual = [row for row in rows if row.get("actual_discarded", True)]
                summary[f"{stage}.{direction}"] = {
                    "candidate_coefficient_count": len(rows),
                    "candidate_pathwise_l1": float(sum(row["abs"] for row in rows)),
                    "candidate_pathwise_l2": math.sqrt(sum(row["abs"] ** 2 for row in rows)),
                    "actual_discarded_coefficient_count": len(actual),
                    "actual_discarded_pathwise_l1": float(sum(row["abs"] for row in actual)),
                    "actual_discarded_pathwise_l2": math.sqrt(sum(row["abs"] ** 2 for row in actual)),
                    "max_candidate_abs": max((row["abs"] for row in rows), default=0.0),
                    "cap_reason_counts": {
                        reason: sum(reason in row.get("cap_reasons", {}) for row in rows)
                        for reason in ("max_states", "max_depth")
                    } if stage == "closure_caps" else {},
                }
        coefficient_ledger["summary"] = summary
        return coefficient_ledger

    old_raw_run = RAW.run
    old_apply_hit_branch = PW.apply_hit_branch
    spin_wall_stats = {
        "count": 0, "l1": 0.0, "l2_sq": 0.0, "max_abs": 0.0,
        "by_edge": {}, "samples": [],
        "max_intermediate_spin_reached": 0.0,
        "intermediate_above_3_5_count": 0,
        "intermediate_above_3_5_l1": 0.0,
        "intermediate_above_3_5_l2_sq": 0.0,
        "intermediate_samples": [],
    }
    PW.apply_hit_branch = _make_spin_wall_probe(
        old_apply_hit_branch, PW.EIDX, Jmax2, spin_wall_stats
    )
    sine_engine_result = {}

    def closure_executor(_engine_node=0, a=None, b=None, c=None):
        nonlocal rejected_count, rejected_norm2, truncated, omitted_support
        while queue:
            key = queue.popleft()
            f_node_maps, ad_node_maps = [], []
            f_node_rows, ad_node_rows = [], []
            for source_v in source_nodes:
                f_cov, a_cov, diag = node_kkv_maps(
                    {key: 1.0 + 0j}, source_v, Jmax2, PW, CV, RAW
                )
                f, nf, nfnorm = covariant_to_gauss_map(f_cov, source_v)
                ad, na, nanorm = covariant_to_gauss_map(a_cov, source_v)
                rejected_count += nf + na
                rejected_norm2 += nfnorm**2 + nanorm**2
                f_node_rows.extend(node_contribution_rows("forward", source_v, key, f))
                ad_node_rows.extend(node_contribution_rows("adjoint", source_v, key, ad))
                f_node_maps.append(f)
                ad_node_maps.append(ad)
                for name, value in diag.items():
                    leakage_max[f"node_{source_v}_{name}"] = max(
                        leakage_max.get(f"node_{source_v}_{name}", 0.0),
                        float(value),
                    )
                    aggregate_leakage[name] = max(
                        aggregate_leakage.get(name, 0.0), float(value)
                    )

            # Sum every complete node contribution with complex fsum before
            # any downstream support or amplitude threshold is considered.
            f_total = sum_sparse_maps(f_node_maps)
            ad_total = sum_sparse_maps(ad_node_maps)

            # Preserve the pre-sum node coefficients and post-sum coherent
            # vectors as distinct evidence. Pathwise and coherent norms are
            # intentionally never merged into one loss number.
            coherent_sum_ledger.append({
                "input_key": repr(key),
                "forward": {
                    "node_contributions": f_node_rows,
                    **coherent_metrics(f_node_rows, f_total),
                },
                "adjoint": {
                    "node_contributions": ad_node_rows,
                    **coherent_metrics(ad_node_rows, ad_total),
                },
            })

            # Record every coefficient removed after coherent node summation.
            # Forward and direct-adjoint ledgers are kept strictly separate.
            for output, amplitude in f_total.items():
                if abs(amplitude) <= amp_tol:
                    record_dropped("amp_tol", "forward", key, output, amplitude,
                                   threshold=amp_tol)
            for output, amplitude in ad_total.items():
                if abs(amplitude) <= amp_tol:
                    record_dropped("amp_tol", "adjoint", key, output, amplitude,
                                   threshold=amp_tol)
            f_total = {k: z for k, z in f_total.items() if PRUNE.keep(z, amp_tol)}
            ad_total = {k: z for k, z in ad_total.items() if PRUNE.keep(z, amp_tol)}
            forward_columns.append(f_total)
            adjoint_columns.append(ad_total)

            # A state-cap/depth-cap rejection is recorded coefficient by
            # coefficient, not merely as a count of omitted support keys.
            for output in set(f_total) | set(ad_total):
                if output not in seen:
                    cap_reasons = _closure_cap_reasons(
                        len(closure_keys), depth[key], max_states, max_depth
                    )
                    if cap_reasons:
                        truncated = True
                        omitted_support += 1
                        cap_payload = {
                            "max_states": max_states,
                            "max_depth": max_depth,
                        }
                        if output in f_total:
                            record_dropped("closure_caps", "forward", key, output,
                                           f_total[output], depth_value=depth[key],
                                           cap=cap_payload, cap_reasons=cap_reasons)
                        if output in ad_total:
                            record_dropped("closure_caps", "adjoint", key, output,
                                           ad_total[output], depth_value=depth[key],
                                           cap=cap_payload, cap_reasons=cap_reasons)
                        continue
                    seen.add(output)
                    depth[output] = depth[key] + 1
                    closure_keys.append(output)
                    queue.append(output)

        # The intercepted result is adapted to the return schema expected by
        # LSINE.run. That wrapper supplies the physical sine K=[V,H_E] patch
        # and restores all monkey-patched functions in a finally block.
        norm_f = sparse_frobenius(forward_columns)
        return {
            "passed": bool(not truncated and rejected_count == 0),
            "output_support": sum(len(x) for x in forward_columns),
            "output_norm": norm_f,
            "final_nonscalar_J_weight_fraction": 0.0 if rejected_norm2 == 0 else 1.0,
            "max_spin_reached": max(
                max((max(key[0]) for key in closure_keys), default=0) / 2,
                spin_wall_stats["max_intermediate_spin_reached"],
            ),
            "max_diagnostics": {
                "CV_complete_basis_leakage": aggregate_leakage.get("CV_complete_basis_leakage", 0.0),
                "CK_outer_complete_basis_leakage": aggregate_leakage.get("CK_outer_complete_basis_leakage", 0.0),
                "CK_internal_volume_sector_leakage": aggregate_leakage.get("CK_internal_volume_sector_leakage", 0.0),
                "CK_complete_charge_basis_leakage": aggregate_leakage.get("CK_complete_charge_basis_leakage", 0.0),
            },
        }

    RAW.run = closure_executor
    from contextlib import nullcontext
    with PRUNE.capture() as lower_level_pruning:
        pruning_context = PRUNE.bypass() if pruning_mode == "shadow_no_prune" else nullcontext()
        cache_context = PRUNE.cache_bypass() if cache_mode == "off" else nullcontext()
        with cache_context:
            with pruning_context:
                try:
                    # LSINE.run installs the exact sine-Hermitian K=[V,H_E] engine.
                    # The context-local cache policy applies through the complete
                    # nested operator dependency chain.
                    sine_engine_result = LSINE.run(source_nodes[0])
                    spectral_sensitivity = ZVM.spectral_audit_snapshot()
                finally:
                    RAW.run = old_raw_run
                    PW.apply_hit_branch = old_apply_hit_branch
    lower_level_pruning_ledger = lower_level_pruning.snapshot()
    lower_level_pruning_ledger["mode"] = pruning_mode
    lower_level_pruning_ledger["coverage"] = build_pruning_coverage(
        static_pruning_inventory
    )
    lower_level_pruning_ledger["coverage"]["instrumented_modules"] = [
        "k5_peter_weyl_safe_hda_column.py",
        "peter_weyl_covariant_volume_leg_gate.py",
        "peter_weyl_covariant_K_leg_gate.py",
        "peter_weyl_covariant_composition_gate.py",
        "peter_weyl_lorentzian_K_block_gate.py",
        "peter_weyl_lorentzian_ordered_triple_gate.py",
        "peter_weyl_euclidean_sine_ordering_gate.py",
        "peter_weyl_covariant_K_projection_audit_gate.py",
        "bqg_v10_w5_lorentzian_adjoint_closure.py",
    ]
    lower_level_pruning_ledger["cache_policy"] = {
        "cleared_lru_cache_function_count": int(cleared_cache_functions),
        "cleared_before_capture": True,
        "mode": cache_mode,
        "policy": (
            "all W5 dependency LRU caches are audit-aware; cache_mode=off bypasses "
            "every audit-aware cache context-locally while preserving operator bodies"
            if cache_mode == "off" else
            "all imported W5 dependency LRU caches are cleared before the audited closure and then enabled"
        ),
    }

    complete = not truncated and not queue and len(forward_columns) == len(closure_keys)
    if not complete:
        # Preserve the real amplitudes computed before a closure cap is hit.
        # This makes bounded pilots scientifically useful without pretending
        # that a truncated basis supports a full matrix-adjoint verdict.
        def top_amplitudes(columns, limit=10):
            rows = []
            for col_index, state in enumerate(columns):
                ranked = sorted(state.items(), key=lambda kv: abs(kv[1]), reverse=True)
                rows.append({
                    "input_key": repr(closure_keys[col_index]),
                    "support_size": len(state),
                    "column_norm": sparse_norm(state),
                    "largest_amplitudes": [
                        {"output_key": repr(key), "real": float(z.real),
                         "imag": float(z.imag), "abs": float(abs(z))}
                        for key, z in ranked[:limit]
                    ],
                })
            return rows

        return {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "source_nodes": list(source_nodes),
            "pruning_mode": pruning_mode,
            "cache_mode": cache_mode,
            "Jmax2": Jmax2,
            "Jmax": Jmax2 / 2,
            "seed_count": len(seeds),
            "seed_Ks": list(selected_seed_Ks),
            "closure_size_discovered": len(closure_keys),
            "columns_computed": len(forward_columns),
            "max_states": max_states,
            "max_depth": max_depth,
            "omitted_new_support_keys": omitted_support,
            "coefficient_loss_ledger": summarize_ledger(),
            "lower_level_pruning_ledger": lower_level_pruning_ledger,
            "coherent_node_sum_ledger": coherent_sum_ledger,
            "spectral_sensitivity": spectral_sensitivity,
            "partial_diagnostics": {
                "forward_frobenius_norm": sparse_frobenius(forward_columns),
                "direct_adjoint_frobenius_norm": sparse_frobenius(adjoint_columns),
                "forward_column_support_sizes": [len(x) for x in forward_columns],
                "direct_adjoint_column_support_sizes": [len(x) for x in adjoint_columns],
                "forward_column_norms": [sparse_norm(x) for x in forward_columns],
                "direct_adjoint_column_norms": [sparse_norm(x) for x in adjoint_columns],
                "scalar_projection_rejected_terms": rejected_count,
                "scalar_projection_rejected_norm": math.sqrt(rejected_norm2),
                **_spin_wall_diagnostics(spin_wall_stats),
                "component_leakage_by_node_and_type": leakage_max,
                "component_leakage_max_by_type": aggregate_leakage,
            },
            "partial_forward_columns_top_amplitudes": top_amplitudes(forward_columns),
            "partial_direct_adjoint_columns_top_amplitudes": top_amplitudes(adjoint_columns),
            "adjoint_residual": None,
            "claim_scope": "real computed columns and diagnostics on a truncated reachable subspace; no full adjoint-consistency verdict",
        }

    basis_keys = sorted(closure_keys, key=repr)
    key_index = {key: i for i, key in enumerate(basis_keys)}
    # Reorder computed columns from BFS order to the deterministic basis order.
    bfs_index = {key: i for i, key in enumerate(closure_keys)}
    forward_sorted = [forward_columns[bfs_index[key]] for key in basis_keys]
    adjoint_sorted = [adjoint_columns[bfs_index[key]] for key in basis_keys]

    # Verify every output is represented in the same closed basis.
    outside = [
        key for cols in (forward_sorted, adjoint_sorted)
        for col in cols for key in col if key not in key_index
    ]
    if outside:
        raise AssertionError("closure algorithm terminated with outputs outside basis")

    expected_adj = _matrix_adjoint_expected_columns(forward_sorted, basis_keys, key_index)
    err, scale = sparse_matrix_residual(adjoint_sorted, expected_adj)
    relative = err / max(scale, 1e-30)
    bnorm = sparse_frobenius(forward_sorted)
    bdag_norm = sparse_frobenius(adjoint_sorted)
    sym_norm = sparse_sum_difference_norm(forward_sorted, expected_adj, sign=+1.0)
    anti_norm = sparse_sum_difference_norm(forward_sorted, expected_adj, sign=-1.0)
    basis_payload = "\n".join(repr(k) for k in basis_keys).encode("utf-8")
    basis_id = hashlib.sha256(basis_payload).hexdigest()

    def encode_columns(cols):
        encoded = []
        for col in cols:
            encoded.append([
                [key_index[key], float(z.real), float(z.imag)]
                for key, z in sorted(col.items(), key=lambda kv: key_index[kv[0]])
            ])
        return encoded

    Bsym_columns, Banti_columns, BHsym_hbar1_columns = sparse_ordering_parts(
        forward_sorted, adjoint_sorted
    )
    hard_leak_names = (
        "CV_complete_basis_leakage",
        "CK_outer_complete_basis_leakage",
        "CK_internal_volume_sector_leakage",
    )
    max_leak = max(
        (value for name, value in leakage_max.items()
         if any(name.endswith(hard_name) for hard_name in hard_leak_names)),
        default=0.0,
    )
    historical_charge_diag = sine_engine_result.get(
        "historical_primitive_charge_basis_diagnostic",
        {"value": aggregate_leakage.get("CK_complete_charge_basis_leakage", 0.0),
         "hard_acceptance": False},
    )
    max_spin_reached = max(
        max((max(key[0]) for key in basis_keys), default=0) / 2,
        spin_wall_stats["max_intermediate_spin_reached"],
    )
    passed = closure_acceptance(
        forward_norm=bnorm,
        adjoint_norm=bdag_norm,
        adjoint_relative_residual=relative,
        scalar_projection_rejected_terms=rejected_count,
        max_component_leakage=max_leak,
        max_spin_reached=max_spin_reached,
        sine_engine_result=sine_engine_result,
        spin_wall_omitted_branches=spin_wall_stats["count"],
        max_intermediate_spin_reached=spin_wall_stats["max_intermediate_spin_reached"],
    )
    algebraic_status = "PASS_ADJOINT_CONSISTENCY" if passed else "FAIL_ADJOINT_OR_LEAKAGE_GATE"
    reported_status = (
        "INCOMPLETE_NO_PRUNING_COVERAGE"
        if pruning_mode == "shadow_no_prune"
        else algebraic_status
    )
    return {
        "status": reported_status,
        "algebraic_status": algebraic_status,
        "lower_level_pruning_coverage_complete": False,
        "source_nodes": list(source_nodes),
        "pruning_mode": pruning_mode,
        "cache_mode": cache_mode,
        "Jmax2": Jmax2,
        "Jmax": Jmax2 / 2,
        "seed_count": len(seeds),
        "seed_Ks": list(selected_seed_Ks),
        "closure_size": len(basis_keys),
        "closure_depth_max": max(depth.values(), default=0),
        "basis_id_sha256": basis_id,
        "basis_keys_repr": [repr(k) for k in basis_keys],
        "forward_B_columns_sparse": encode_columns(forward_sorted),
        "direct_Bdag_columns_sparse": encode_columns(adjoint_sorted),
        "hermitian_ordering_candidate": {
            "B_comm_symmetric_definition": "0.5*(B_comm + B_comm_dagger)",
            "B_comm_antihermitian_definition": "0.5*(B_comm - B_comm_dagger)",
            "B_H_symmetric_definition": "-B_comm_symmetric/hbar^2",
            "B_H_symmetric_numeric_columns_at_hbar_1": encode_columns(BHsym_hbar1_columns),
            "B_comm_symmetric_columns": encode_columns(Bsym_columns),
            "B_comm_antihermitian_columns": encode_columns(Banti_columns),
            "symmetric_part_is_hermitian_if_direct_adjoint_is_exact": True,
            "antihermitian_part_is_antihermitian_if_direct_adjoint_is_exact": True,
            "finite_precision_qualification": "The direct-adjoint relative residual above controls the numerical departure from exact Hermiticity/anti-Hermiticity.",
            "physical_ordering_selected": False,
            "external_beta_kappa_planck_trace_regulator_prefactor_fixed": False,
            "warning": "This is the minimal symmetric completion of the computed raw KKV ordering, not a proof of the unique physical BQG Thiemann ordering."
        },
        "coefficient_loss_ledger": summarize_ledger(),
        "lower_level_pruning_ledger": lower_level_pruning_ledger,
        "coherent_node_sum_ledger": coherent_sum_ledger,
        "spectral_sensitivity": spectral_sensitivity,
        "diagnostics": {
            "direct_adjoint_relative_residual": relative,
            "direct_adjoint_absolute_residual": err,
            "forward_frobenius_norm": bnorm,
            "direct_adjoint_frobenius_norm": bdag_norm,
            "hermitian_part_frobenius_norm": 0.5 * sym_norm,
            "antihermitian_part_frobenius_norm": 0.5 * anti_norm,
            "scalar_projection_rejected_terms": rejected_count,
            "scalar_projection_rejected_norm": math.sqrt(rejected_norm2),
            "max_component_leakage": max_leak,
            "max_spin_reached": max_spin_reached,
            **_spin_wall_diagnostics(spin_wall_stats),
            "sine_engine_gate_passed": sine_engine_result.get("passed") is True,
            "sine_engine_gate_diagnostics": sine_engine_result.get("max_diagnostics", {}),
            "component_leakage_by_type": leakage_max,
            "historical_primitive_charge_basis_diagnostic": historical_charge_diag,
            "hard_leakage_gate_names": list(hard_leak_names),
        },
        "physical_status": {
            "physical_hermitian_ordering_selected": False,
            "beta_kappa_hbar_prefactor_fixed": False,
            "full_physical_HL": False,
            "physical_joint_kernel_verdict": False,
        },
        "node_sum": "unweighted coherent sum of selected source-node KKV maps and their direct adjoints",
        "claim_scope": "finite W5-seeded reachable component at declared Jmax2; raw structural node-summed KKV and independently ordered adjoint only",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nodes", type=int, nargs="+", default=list(NODES))
    ap.add_argument("--Jmax2", type=int, default=10)
    ap.add_argument("--max-states", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument(
        "--unbounded-closure", action="store_true",
        help="disable both reachable-state and path-depth caps (finite Jmax2 still applies)",
    )
    ap.add_argument("--amp-tol", type=float, default=1e-12)
    ap.add_argument("--spectral-zero-factor", type=float, default=1000.0)
    ap.add_argument("--cache-mode", choices=("on", "off"), default="on")
    ap.add_argument("--seed-K", type=int, nargs="+", choices=W5_K, default=list(W5_K),
                    help="W5 seed labels to include; default runs all five")
    ap.add_argument("--pruning-mode", choices=("production", "shadow_no_prune"), default="production")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    max_states = None if args.unbounded_closure else args.max_states
    max_depth = None if args.unbounded_closure else args.max_depth
    result = run_closure(
        args.nodes, args.Jmax2, max_states, max_depth,
        args.amp_tol, args.pruning_mode, args.spectral_zero_factor,
        args.cache_mode, args.seed_K,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "source_nodes": result["source_nodes"],
        "seed_Ks": result.get("seed_Ks"),
        "Jmax2": result["Jmax2"],
        "closure_size": result.get("closure_size", result.get("closure_size_discovered")),
        "direct_adjoint_relative_residual": result.get("diagnostics", {}).get("direct_adjoint_relative_residual"),
        "max_component_leakage": result.get("diagnostics", {}).get("max_component_leakage"),
        "output": str(args.output),
        "claim_scope": result["claim_scope"],
    }, indent=2))
    return 0 if result["status"] == "PASS_ADJOINT_CONSISTENCY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
