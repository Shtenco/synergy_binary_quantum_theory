#!/usr/bin/env python3
"""Direct all-M SU(2) commutator audit of the actual W5 sine-H_E operator.

This deliberately bypasses the cached highest-weight wrapper. It evaluates the
same production H_E map on every magnetic component of a fixed charged J=1
four-leg state at each of the five K5 source nodes and checks
[J_+, H_E]=[J_-, H_E]=0 on the complete sparse output.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bqg_pruning_audit as PRUNE
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_sine_composition_gate as SINE


def apply_local_ladder(state, node, direction):
    if direction not in (-1, 1):
        raise ValueError("direction must be -1 or +1")
    out = {}
    for (spins, labels), amp in state.items():
        lab = labels[node]
        if lab[0] != "C":
            raise ValueError(f"node {node} is not in the charged basis")
        _, j2, m2, k12, k34 = lab
        j = j2 / 2.0
        m = m2 / 2.0
        if direction > 0:
            if m2 >= j2:
                continue
            factor = math.sqrt(max(0.0, j * (j + 1.0) - m * (m + 1.0)))
            next_m2 = m2 + 2
        else:
            if m2 <= -j2:
                continue
            factor = math.sqrt(max(0.0, j * (j + 1.0) - m * (m - 1.0)))
            next_m2 = m2 - 2
        labs = list(labels)
        labs[node] = CK.label_charge(j2, next_m2, k12, k34)
        key = (spins, tuple(labs))
        out[key] = out.get(key, 0j) + amp * factor
    return out


def add_sparse(dst, src, scale=1.0):
    for key, amp in src.items():
        dst[key] = dst.get(key, 0j) + scale * amp


def apply_local_jz(state, node):
    out = {}
    for (spins, labels), amp in state.items():
        lab = labels[node]
        if lab[0] != "C":
            raise ValueError(f"node {node} is not in the charged basis")
        m = int(lab[2]) / 2.0
        out[(spins, labels)] = complex(amp) * m
    return out


def norm(state):
    if not state:
        return 0.0
    magnitudes = []
    for amp in state.values():
        z = complex(amp)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError("operator output contains a non-finite amplitude")
        magnitudes.append(abs(z))
    value = math.hypot(*magnitudes)
    if not math.isfinite(value):
        raise ValueError("operator output norm is non-finite")
    return value


def current_rss_bytes():
    """Return current resident memory on Linux runners, or None if unavailable."""
    try:
        with open("/proc/self/statm", "r", encoding="ascii") as stream:
            resident_pages = int(stream.read().split()[1])
        return resident_pages * int(os.sysconf("SC_PAGE_SIZE"))
    except (OSError, ValueError, IndexError, AttributeError):
        return None


def _cache_modules():
    return (
        CK, CK.PW, CK.CH, CK.CV, CK.KG,
    )


def cache_snapshot(modules=None):
    """Capture bounded, JSON-safe statistics for every exposed LRU cache."""
    modules = _cache_modules() if modules is None else tuple(modules)
    rows = []
    for module in modules:
        for name, value in vars(module).items():
            info = getattr(value, "cache_info", None)
            if callable(info):
                stats = info()
                rows.append({
                    "module": getattr(module, "__name__", type(module).__name__),
                    "function": name,
                    "hits": int(stats.hits),
                    "misses": int(stats.misses),
                    "maxsize": stats.maxsize,
                    "currsize": int(stats.currsize),
                })
    return rows


def clear_operator_caches(modules=None):
    """Clear all audit-aware LRU caches in the active W5 operator modules."""
    modules = _cache_modules() if modules is None else tuple(modules)
    cleared = []
    for module in modules:
        for name, value in vars(module).items():
            clear = getattr(value, "cache_clear", None)
            if callable(clear):
                clear()
                cleared.append(f"{getattr(module, '__name__', type(module).__name__)}.{name}")
    return cleared


def require_nonzero_output(state, *, label):
    """Reject empty, numerically zero, or non-finite operator results."""
    if not state:
        raise AssertionError(f"{label} returned an empty sparse state")
    value = norm(state)
    if not math.isfinite(value) or value == 0.0:
        raise AssertionError(f"{label} has no resolvable nonzero norm")
    return value


def relative_error(a, b):
    if not a or not b:
        raise AssertionError("commutator comparison cannot use an empty state")
    error = math.hypot(*(
        abs(a.get(key, 0j) - b.get(key, 0j))
        for key in set(a) | set(b)
    ))
    scale = max(norm(a), norm(b))
    if not math.isfinite(error) or not math.isfinite(scale) or scale == 0.0:
        raise AssertionError("commutator comparison has a zero or invalid norm scale")
    return error / scale, error, scale


def wrong_sector_fraction(state, *, node, expected_j2, expected_m2):
    """Measure forbidden (J,M) weight against the full nonzero output norm."""
    magnitudes = [abs(complex(amp)) for amp in state.values()]
    scale = max(magnitudes, default=0.0)
    if not math.isfinite(scale) or scale == 0.0:
        raise AssertionError("charged-sector output has zero or invalid total weight")
    total = math.fsum((magnitude / scale) ** 2 for magnitude in magnitudes)
    wrong = 0.0
    for ((_, labels), amp), magnitude in zip(state.items(), magnitudes):
        label = labels[node]
        if label[0] != "C":
            raise AssertionError(f"operator output lost charged basis at node {node}")
        if int(label[1]) != expected_j2 or int(label[2]) != expected_m2:
            wrong += (magnitude / scale) ** 2
    fraction = wrong / total
    if not math.isfinite(fraction):
        raise AssertionError("wrong-sector fraction is non-finite")
    return fraction


def spectator_non_singlet_weight_fractions(state, *, source_node):
    """Measure non-singlet weight at every spectator vertex of the full output.

    The charged complete-JM basis is needed while each graph-changing term is
    assembled, but the fully summed scalar H_E must leave every non-source
    vertex in its Gauss-singlet (J=0) sector. A charged label with J2=0 is also
    a singlet and is accepted; only J2>0 contributes to leakage.
    """
    if not isinstance(source_node, int) or isinstance(source_node, bool) or not 0 <= source_node < 5:
        raise ValueError("source_node must be an integer in [0, 4]")
    if not state:
        raise AssertionError("spectator Gauss check cannot use an empty output")
    magnitudes = []
    for amp in state.values():
        z = complex(amp)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError("spectator Gauss check contains a non-finite amplitude")
        magnitudes.append(abs(z))
    scale = max(magnitudes, default=0.0)
    if scale == 0.0 or not math.isfinite(scale):
        raise AssertionError("spectator Gauss check has zero or invalid output norm")
    total = math.fsum((magnitude / scale) ** 2 for magnitude in magnitudes)
    fractions = {}
    for vertex in range(5):
        if vertex == source_node:
            continue
        wrong = 0.0
        for ((_, labels), _amp), magnitude in zip(state.items(), magnitudes):
            label = labels[vertex]
            if label[0] == "G":
                continue
            if label[0] != "C":
                raise AssertionError(f"operator output has an unknown basis label at spectator node {vertex}")
            if int(label[1]) != 0:
                wrong += (magnitude / scale) ** 2
        fraction = wrong / total
        if not math.isfinite(fraction):
            raise AssertionError(f"spectator Gauss leakage is non-finite at node {vertex}")
        fractions[str(vertex)] = fraction
    return fractions


def make_input(M2, node=0, K12=0, K34=2):
    if not isinstance(node, int) or isinstance(node, bool) or not 0 <= node < 5:
        raise ValueError("node must be an integer in [0, 4]")
    if not isinstance(M2, int) or isinstance(M2, bool) or M2 not in (-2, 0, 2):
        raise ValueError("J=1 doubled magnetic label M2 must be one of -2, 0, 2")
    if not isinstance(K12, int) or isinstance(K12, bool) or not isinstance(K34, int) or isinstance(K34, bool):
        raise ValueError("recoupling labels must be integers")
    spins = (4,) * 10
    labels = [("G", 0) for _ in range(5)]
    labels[node] = ("C", 2, int(M2), int(K12), int(K34))
    return spins, tuple(labels)


def json_default(value):
    """Losslessly encode complex diagnostics as explicit [real, imag] pairs."""
    if isinstance(value, (complex, np.complexfloating)):
        z = complex(value)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError("non-finite complex value cannot enter evidence JSON")
        return [float(z.real), float(z.imag)]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError("non-finite floating value cannot enter evidence JSON")
        return result
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def run(Jmax2=10, pruning_mode="production", nodes=(0,), channels=((0, 2), (2, 0)), m_pair=None):
    if pruning_mode not in ("production", "shadow_no_prune"):
        raise ValueError("invalid pruning_mode")
    if Jmax2 < 4:
        raise ValueError("Jmax2 must be at least 4 for this fixed W5 seed")
    nodes = tuple(nodes)
    if not nodes or any(not isinstance(node, int) or isinstance(node, bool) or not 0 <= node < 5 for node in nodes):
        raise ValueError("nodes must be a non-empty subset of [0, 1, 2, 3, 4]")
    if len(set(nodes)) != len(nodes):
        raise ValueError("nodes must not contain duplicates")
    channels = tuple((int(a), int(b)) for a, b in channels)
    if not channels or len(set(channels)) != len(channels):
        raise ValueError("channels must be a non-empty unique set of (K12,K34) pairs")
    if any(channel not in ((0, 2), (2, 0)) for channel in channels):
        raise ValueError("W5 charged channels must be selected from ((0,2),(2,0))")
    m_values = (-2, 0, 2) if m_pair is None else tuple(m_pair)
    if m_values not in ((-2, 0, 2), (-2, 0), (0, 2)):
        raise ValueError("m_pair must be None, (-2,0), or (0,2)")
    directions = (1, -1)
    cases = []
    all_outputs = {}
    superposition_rows = []
    cache_memory_rows = []
    context = PRUNE.bypass() if pruning_mode == "shadow_no_prune" else None

    def evaluate():
        for node in nodes:
            he_by_channel_m = {}
            for K12, K34 in channels:
                he_by_m = {}
                for M2 in m_values:
                    projection_diagnostics = {}
                    with PRUNE.capture() as collector:
                        result = SINE.complete_HE_sine(
                            make_input(M2, node, K12, K34), node, Jmax2,
                            charged_nodes=(node,),
                            projection_diagnostics=projection_diagnostics,
                        )
                    require_nonzero_output(
                        result[0], label=f"H_E(node={node}, K12={K12}, K34={K34}, M2={M2})"
                    )
                    all_outputs[(node, K12, K34, M2)] = result[0]
                    rss_before_clear = current_rss_bytes()
                    caches_before_clear = cache_snapshot()
                    # The direct sine-H_E path does not use the cached
                    # highest-weight wrapper. It can still populate lower-level
                    # recoupling/volume caches. Clear them between independent
                    # magnetic columns so memory cannot grow across M sectors.
                    clear_operator_caches()
                    cache_memory_rows.append({
                        "source_node": node,
                        "K12": K12,
                        "K34": K34,
                        "M2": M2,
                        "current_rss_bytes_before_clear": rss_before_clear,
                        "current_rss_bytes_after_clear": current_rss_bytes(),
                        "caches_before_clear": caches_before_clear,
                        "caches_after_clear": cache_snapshot(),
                    })
                    he_by_m[M2] = {
                        "state": result[0],
                        "norm": norm(result[0]),
                        "volume_leakage": float(result[1]),
                        "basis_leakage": float(result[2]),
                        "basis_leakage_by_vertex": projection_diagnostics,
                        "pruning_events": collector.snapshot()["event_count"],
                    }
                he_by_channel_m[(K12, K34)] = he_by_m

            # Independent implementation-level linearity check: place complex
            # coefficients on the input branches, then run the real sequence,
            # recoupling, projection and coherent-accumulation path.
            # Keep the real operator-composition check inside each shard.
            # With both recoupling channels selected, mix their M=0 inputs.
            # With one channel selected, mix its M=-1 and M=+1 inputs instead.
            if len(channels) >= 2:
                (ka, kb) = channels[:2]
                key_a = make_input(0, node, *ka)
                key_b = make_input(0, node, *kb)
            else:
                k12, k34 = channels[0]
                key_a = make_input(m_values[0], node, k12, k34)
                key_b = make_input(m_values[-1], node, k12, k34)
            alpha = 0.7 + 0.2j
            beta = -0.3 + 0.9j
            mixed_projection_diagnostics = {}
            mixed, mixed_vleak, mixed_bleak = SINE.complete_HE_sine_state(
                {key_a: alpha, key_b: beta},
                node, Jmax2, charged_nodes=(node,),
                projection_diagnostics=mixed_projection_diagnostics,
            )
            require_nonzero_output(
                mixed, label=f"H_E complex superposition at node {node}"
            )
            expected_terms = {}
            for key, coefficient in ((key_a, alpha), (key_b, beta)):
                label = key[1][node]
                unit_output = all_outputs[(node, label[3], label[4], label[2])]
                for output_key, amplitude in unit_output.items():
                    expected_terms.setdefault(output_key, []).append(coefficient * amplitude)
            expected = {
                output_key: complex(
                    math.fsum(value.real for value in terms),
                    math.fsum(value.imag for value in terms),
                )
                for output_key, terms in expected_terms.items()
                if any(value != 0j for value in terms)
            }
            super_rel, super_abs, super_scale = relative_error(mixed, expected)
            superposition_rows.append({
                "source_node": node,
                "input_basis_states": 2,
                "complex_coefficients": [
                    [alpha.real, alpha.imag], [beta.real, beta.imag]
                ],
                "output_support_dimension": len(mixed),
                "relative_residual": super_rel,
                "absolute_residual": super_abs,
                "comparison_scale": super_scale,
                "volume_leakage": float(mixed_vleak),
                "basis_leakage": float(mixed_bleak),
                "basis_leakage_by_vertex": mixed_projection_diagnostics,
                "scope": (
                    "Actual complex coefficients propagated through branch construction, "
                    "graph-changing sequences, recoupling projection and coherent accumulation; "
                    "compared with independently evaluated unit-column outputs"
                ),
            })

            for K12, K34 in channels:
              he_by_m = he_by_channel_m[(K12, K34)]
              for direction in directions:
                source_M2_values = tuple(m for m in m_values if m + 2 * direction in m_values)
                for M2 in source_M2_values:
                    input_key = make_input(M2, node, K12, K34)
                    raised_input = apply_local_ladder({input_key: 1.0 + 0j}, node, direction)
                    if len(raised_input) != 1:
                        raise AssertionError("ladder action should produce one adjacent M sector")
                    _, ladder_factor = next(iter(raised_input.items()))
                    H_input = he_by_m[M2]["state"]
                    H_raised = he_by_m[M2 + 2 * direction]["state"]
                    lhs = apply_local_ladder(H_input, node, direction)
                    rhs = {key: amp * ladder_factor for key, amp in H_raised.items()}
                    rel, abs_err, scale = relative_error(lhs, rhs)
                    jz_output = apply_local_jz(H_input, node)
                    jz_expected = {key: amp * (M2 / 2.0) for key, amp in H_input.items()}
                    jz_difference = {
                        key: jz_output.get(key, 0j) - jz_expected.get(key, 0j)
                        for key in set(jz_output) | set(jz_expected)
                    }
                    jz_abs_err = norm(jz_difference)
                    jz_scale = norm(H_input)
                    jz_rel = jz_abs_err / jz_scale
                    wrong_fraction = wrong_sector_fraction(
                        H_input, node=node, expected_j2=2, expected_m2=M2
                    )
                    spectator_fractions = spectator_non_singlet_weight_fractions(
                        H_input, source_node=node
                    )
                    raised_spectator_fractions = spectator_non_singlet_weight_fractions(
                        H_raised, source_node=node
                    )
                    max_spectator_fraction = max(
                        max(spectator_fractions.values(), default=0.0),
                        max(raised_spectator_fractions.values(), default=0.0),
                    )
                    cases.append({
                        "source_node": node,
                        "K12": K12,
                        "K34": K34,
                        "direction": "J_plus" if direction > 0 else "J_minus",
                        "input_M2": M2,
                        "output_M2": M2 + 2 * direction,
                        "ladder_factor": ladder_factor,
                        "input_HE_support": len(H_input),
                        "raised_HE_support": len(H_raised),
                        "input_HE_norm": he_by_m[M2]["norm"],
                        "raised_HE_norm": he_by_m[M2 + 2 * direction]["norm"],
                        "commutator_relative_residual": rel,
                        "commutator_absolute_residual": abs_err,
                        "comparison_scale": scale,
                        "jz_commutator_relative_residual": jz_rel,
                        "jz_commutator_absolute_residual": jz_abs_err,
                        "jz_comparison_scale": jz_scale,
                        "wrong_JM_weight_fraction": wrong_fraction,
                        "wrong_JM_relative_amplitude": math.sqrt(wrong_fraction),
                        "spectator_non_singlet_weight_fractions": spectator_fractions,
                        "raised_spectator_non_singlet_weight_fractions": raised_spectator_fractions,
                        "max_spectator_non_singlet_weight_fraction": max_spectator_fraction,
                        "volume_leakage": max(
                            he_by_m[M2]["volume_leakage"],
                            he_by_m[M2 + 2 * direction]["volume_leakage"],
                        ),
                        "basis_leakage": max(
                            he_by_m[M2]["basis_leakage"],
                            he_by_m[M2 + 2 * direction]["basis_leakage"],
                        ),
                        # Preserve the local projection evidence for both
                        # unit columns participating in this ladder comparison.
                        # Aggregate leakage alone cannot identify the vertex or
                        # distinguish charged-basis from Gauss-basis loss.
                        "input_basis_leakage_by_vertex": he_by_m[M2]["basis_leakage_by_vertex"],
                        "raised_basis_leakage_by_vertex": he_by_m[M2 + 2 * direction]["basis_leakage_by_vertex"],
                        "pruning_events": (
                            he_by_m[M2]["pruning_events"]
                            + he_by_m[M2 + 2 * direction]["pruning_events"]
                        ),
                    })

    # Keep the direct, uncached highest-weight-wrapper path: every M column
    # enters complete_HE_sine and the actual graph-changing branch engine.
    # Exact lower-level audit-aware memoization remains enabled. Its cache key
    # includes pruning mode and it replays the original ordered pruning trace,
    # so this does not replace the operator with a cached matrix or suppress
    # pruning evidence. The former global cache_bypass disabled every audited
    # recoupling/branch cache and made this direct control needlessly expensive.
    if context is None:
        evaluate()
    else:
        with context:
            evaluate()

    adjoint_rows = []
    linearity_rows = []
    operator_columns_sparse = []
    for node in nodes:
        basis_keys = [
            make_input(M2, node, K12, K34)
            for K12, K34 in channels for M2 in m_values
        ]
        matrix = np.zeros((len(basis_keys), len(basis_keys)), dtype=np.complex128)
        # Reuse the exact same unit columns that were evaluated under the
        # selected production or shadow-no-prune context.
        for column, (K12, K34, M2) in enumerate(
            ( (k1, k2, m) for k1, k2 in channels for m in m_values )
        ):
            output = all_outputs[(node, K12, K34, M2)]
            for row, key in enumerate(basis_keys):
                matrix[row, column] = output.get(key, 0j)
        output_keys = sorted({
            key for K12, K34 in channels for M2 in m_values
            for key in all_outputs[(node, K12, K34, M2)]
        }, key=repr)
        output_index = {key: row for row, key in enumerate(output_keys)}
        # Persist the complete sparse outputs, including support outside the
        # declared seed subspace. The downstream CI aggregator compares each
        # production column with its shadow_no_prune counterpart.
        for K12, K34 in channels:
            for M2 in m_values:
                sparse_output = all_outputs[(node, K12, K34, M2)]
                operator_columns_sparse.append({
                    "source_node": node,
                    "K12": K12,
                    "K34": K34,
                    "input_M2": M2,
                    "entries": [
                        [repr(key), float(complex(amplitude).real), float(complex(amplitude).imag)]
                        for key, amplitude in sorted(sparse_output.items(), key=lambda item: repr(item[0]))
                    ],
                })
        operator_matrix = np.zeros((len(output_keys), len(basis_keys)), dtype=np.complex128)
        for column, (K12, K34, M2) in enumerate(
            ((k1, k2, m) for k1, k2 in channels for m in m_values)
        ):
            for key, amplitude in all_outputs[(node, K12, K34, M2)].items():
                operator_matrix[output_index[key], column] = amplitude
        matrix_norm = float(np.linalg.norm(matrix, "fro"))
        # The graph-changing map is rectangular: its full codomain is
        # output_keys, not the original seed basis_keys. Keep the seed
        # compression as a diagnostic only, and measure nonzero output on
        # the complete codomain matrix.
        full_output_matrix_norm = float(np.linalg.norm(operator_matrix, "fro"))
        adjoint_abs = float(np.linalg.norm(matrix - matrix.conj().T, "fro"))
        linear_rng = np.random.default_rng(1729 + node)
        x = linear_rng.normal(size=len(basis_keys)) + 1j * linear_rng.normal(size=len(basis_keys))
        y = linear_rng.normal(size=len(basis_keys)) + 1j * linear_rng.normal(size=len(basis_keys))
        alpha = 0.7 + 0.2j
        beta = -0.3 + 0.9j
        lhs = operator_matrix @ (alpha * x + beta * y)
        rhs = alpha * (operator_matrix @ x) + beta * (operator_matrix @ y)
        linearity_abs = float(np.linalg.norm(lhs - rhs))
        linearity_scale = max(float(np.linalg.norm(lhs)), float(np.linalg.norm(rhs)), 1e-30)
        linearity_rel = linearity_abs / linearity_scale
        linearity_rows.append({
            "source_node": node,
            "input_basis_dimension": len(basis_keys),
            "output_support_dimension": len(output_keys),
            "operator_matrix_shape": [len(output_keys), len(basis_keys)],
            "linearity_absolute_residual": linearity_abs,
            "linearity_relative_residual": linearity_rel,
            "scope": "matrix algebra identity only; not an independent test of the branch implementation",
        })
        adjoint_rel = adjoint_abs / matrix_norm if matrix_norm > 1e-14 else None
        adjoint_rows.append({
            "source_node": node,
            "basis_dimension": len(basis_keys),
            "matrix_frobenius_norm": matrix_norm,
            "adjoint_absolute_residual": adjoint_abs,
            "adjoint_relative_residual": adjoint_rel,
            "adjointness_compression_informative": matrix_norm > 1e-14,
            "full_output_matrix_frobenius_norm": full_output_matrix_norm,
            "full_output_codomain_dimension": len(output_keys),
            "scope": (
                "seed-subspace compression only; its adjointness residual is null when "
                "the graph-changing map has no support back into the seed subspace. "
                "The full rectangular map norm is reported separately; this is not "
                "a global self-adjointness proof."
            ),
        })

    informative_adjoint = [
        x["adjoint_relative_residual"] for x in adjoint_rows
        if x["adjointness_compression_informative"]
    ]
    max_adjoint = max(informative_adjoint, default=None)
    adjointness_status = (
        "NOT_FULL_SEED_SUBSPACE_SHARD_DIAGNOSTIC_ONLY" if m_pair is not None
        else ("SEED_SUBSPACE_COMPRESSION_CHECKED"
              if len(informative_adjoint) == len(adjoint_rows)
              else "INCOMPLETE_SEED_SUBSPACE_COMPRESSION")
    )
    full_output_nonzero = all(
        x["full_output_matrix_frobenius_norm"] > 1e-14 for x in adjoint_rows
    )
    max_linearity = max((x["linearity_relative_residual"] for x in linearity_rows), default=0.0)
    max_commutator = max((x["commutator_relative_residual"] for x in cases), default=0.0)
    max_jz_commutator = max((x["jz_commutator_relative_residual"] for x in cases), default=0.0)
    max_wrong_fraction = max((x["wrong_JM_weight_fraction"] for x in cases), default=0.0)
    max_spectator_non_singlet_fraction = max(
        (x["max_spectator_non_singlet_weight_fraction"] for x in cases),
        default=0.0,
    )
    max_actual_superposition = max((x["relative_residual"] for x in superposition_rows), default=float("inf"))
    passed = (
        len(cases) == 2 * (len(m_values) - 1) * len(channels) * len(nodes)
        and max_commutator < 1e-8
        and full_output_nonzero
        and (m_pair is not None or max_adjoint is None or max_adjoint < 1e-8)
        and max_linearity < 1e-12
        and max_actual_superposition < 1e-10
        and len(superposition_rows) == len(nodes)
        # Seed-subspace adjointness compression is diagnostic only for a
        # graph-changing rectangular map. A zero compression is explicitly
        # reported as INCOMPLETE_SEED_SUBSPACE_COMPRESSION; it must not veto
        # the independent full-codomain SU(2) covariance result.
        and max_jz_commutator < 1e-12
        and max_wrong_fraction < 1e-18
        and max_spectator_non_singlet_fraction < 1e-18
        and max((x["volume_leakage"] for x in cases), default=0.0) < 1e-10
        and max((x["basis_leakage"] for x in cases), default=0.0) < 1e-10
    )
    return {
        "schema": "bqg-w5-direct-all-M-su2-commutator-v1",
        "status": "PASS_DIRECT_ALL_M_SU2_COMMUTATOR" if passed else "FAIL_DIRECT_ALL_M_SU2_COMMUTATOR",
        "passed": bool(passed),
        "Jmax2": Jmax2,
        "source_nodes": list(nodes),
        "input_spins_doubled": [4] * 10,
        "input_charged_sector": {"J2": 2, "channels_K12_K34": [list(x) for x in channels]},
        "magnetic_components_tested_M2": list(m_values),
        "adjointness_compression_rows": adjoint_rows,
        "adjointness_compression_status": adjointness_status,
        "max_adjointness_compression_relative_residual": max_adjoint,
        "all_full_output_matrices_nonzero": full_output_nonzero,
        "matrix_superposition_identity_rows": linearity_rows,
        "max_matrix_superposition_identity_relative_residual": max_linearity,
        "actual_complex_superposition_rows": superposition_rows,
        "max_actual_complex_superposition_relative_residual": max_actual_superposition,
        "pruning_mode": pruning_mode,
        "top_level_cached_operator_wrapper_used": False,
        "lower_level_audit_aware_cache_enabled": True,
        "commutator_target": "[J_+, H_E]=[J_-, H_E]=[J_z, H_E]=0 on the actual complete W5 sine-H_E map",
        "max_commutator_relative_residual": max_commutator,
        "max_jz_commutator_relative_residual": max_jz_commutator,
        "max_wrong_JM_weight_fraction": max_wrong_fraction,
        "max_spectator_non_singlet_weight_fraction": max_spectator_non_singlet_fraction,
        "generators_checked": ["J_plus", "J_minus", "J_z"],
        "cases": cases,
        "operator_columns_sparse": operator_columns_sparse,
        "cache_memory_diagnostics": cache_memory_rows,
        "cache_memory_policy": "snapshot and clear all exposed operator LRU caches after each independent M column",
        "operator_columns_encoding": "[repr(output_key), real_amplitude, imaginary_amplitude]; complete output support, no seed-subspace projection",
        "scope": (
            "Direct complete declared-M covariance check on the finite charged J=1 recoupling channels at each requested K5 source node; a resumable shard checks one adjacent-M pair. "
            "The adjointness result is only the compression to the declared seed subspace; it is not a global self-adjointness proof. "
            "It is not a proof for every W5 seed, charged channel, cutoff, graph refinement, or the physical continuum."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--pruning-mode", choices=("production", "shadow_no_prune"), default="production")
    parser.add_argument("--nodes", type=int, nargs="+", default=[0])
    parser.add_argument(
        "--channel", type=str, default=None,
        help="Optional single W5 recoupling channel as K12_K34 (e.g. 0_2); "
             "omitting it runs both declared channels.",
    )
    parser.add_argument("--m-pair", type=str, default=None,
                        help="Optional adjacent doubled-M pair: -2_0 or 0_2; omit for all M.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    channels = ((0, 2), (2, 0))
    m_pair = None
    if args.m_pair is not None:
        if args.m_pair == "neg2_0":
            m_pair = (-2, 0)
        elif args.m_pair == "0_2":
            m_pair = (0, 2)
        else:
            parser.error("--m-pair must be neg2_0 or 0_2")
    if args.channel is not None:
        parts = args.channel.replace("_", ",").split(",")
        if len(parts) != 2:
            parser.error("--channel must have the form K12,K34")
        try:
            channels = ((int(parts[0]), int(parts[1])),)
        except ValueError:
            parser.error("--channel values must be integers")
    result = run(args.Jmax2, args.pruning_mode, args.nodes, channels=channels, m_pair=m_pair)
    # Complex amplitudes are part of the scientific evidence, not an error.
    # Encode them explicitly as [real, imag], rejecting non-finite values.
    text = json.dumps(result, indent=2, default=json_default, allow_nan=False)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
