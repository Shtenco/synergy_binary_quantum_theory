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
from pathlib import Path
import sys

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


def run(Jmax2=10, pruning_mode="production", nodes=(0,), channels=((0, 2), (2, 0))):
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
    directions = (1, -1)
    cases = []
    context = PRUNE.bypass() if pruning_mode == "shadow_no_prune" else None

    def evaluate():
        for node in nodes:
            he_by_channel_m = {}
            for K12, K34 in channels:
                he_by_m = {}
                for M2 in (-2, 0, 2):
                    with PRUNE.capture() as collector:
                        result = SINE.complete_HE_sine(
                            make_input(M2, node, K12, K34), node, Jmax2, charged_nodes=(node,)
                        )
                    require_nonzero_output(
                        result[0], label=f"H_E(node={node}, K12={K12}, K34={K34}, M2={M2})"
                    )
                    he_by_m[M2] = {
                        "state": result[0],
                        "norm": norm(result[0]),
                        "volume_leakage": float(result[1]),
                        "basis_leakage": float(result[2]),
                        "pruning_events": collector.snapshot()["event_count"],
                    }
                he_by_channel_m[(K12, K34)] = he_by_m

            for K12, K34 in channels:
              he_by_m = he_by_channel_m[(K12, K34)]
              for direction in directions:
                source_M2_values = (-2, 0) if direction > 0 else (0, 2)
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
                        "volume_leakage": max(
                            he_by_m[M2]["volume_leakage"],
                            he_by_m[M2 + 2 * direction]["volume_leakage"],
                        ),
                        "basis_leakage": max(
                            he_by_m[M2]["basis_leakage"],
                            he_by_m[M2 + 2 * direction]["basis_leakage"],
                        ),
                        "pruning_events": (
                            he_by_m[M2]["pruning_events"]
                            + he_by_m[M2 + 2 * direction]["pruning_events"]
                        ),
                    })

    if context is None:
        with PRUNE.cache_bypass():
            evaluate()
    else:
        with PRUNE.cache_bypass(), context:
            evaluate()

    adjoint_rows = []
    for node in nodes:
        basis_keys = [
            make_input(M2, node, K12, K34)
            for K12, K34 in channels for M2 in (-2, 0, 2)
        ]
        matrix = np.zeros((len(basis_keys), len(basis_keys)), dtype=np.complex128)
        node_outputs = {
            (K12, K34, M2): he_by_channel_m[(K12, K34)][M2]["state"]
            for K12, K34 in channels for M2 in (-2, 0, 2)
        } if False else None
        # Re-evaluate only the finite matrix elements from the already computed
        # unit columns. The output may leave this small seed subspace; the
        # compression is a necessary (not sufficient) adjointness check.
        for column, (K12, K34, M2) in enumerate(
            ( (k1, k2, m) for k1, k2 in channels for m in (-2, 0, 2) )
        ):
            with PRUNE.capture():
                output = SINE.complete_HE_sine(
                    make_input(M2, node, K12, K34), node, Jmax2, charged_nodes=(node,)
                )[0]
            for row, key in enumerate(basis_keys):
                matrix[row, column] = output.get(key, 0j)
        matrix_norm = float(np.linalg.norm(matrix, "fro"))
        adjoint_abs = float(np.linalg.norm(matrix - matrix.conj().T, "fro"))
        adjoint_rel = adjoint_abs / max(matrix_norm, 1e-30)
        adjoint_rows.append({
            "source_node": node,
            "basis_dimension": len(basis_keys),
            "matrix_frobenius_norm": matrix_norm,
            "adjoint_absolute_residual": adjoint_abs,
            "adjoint_relative_residual": adjoint_rel,
            "scope": "compression to the declared seed subspace; not a global self-adjointness proof",
        })

    max_adjoint = max((x["adjoint_relative_residual"] for x in adjoint_rows), default=0.0)
    max_commutator = max((x["commutator_relative_residual"] for x in cases), default=0.0)
    max_jz_commutator = max((x["jz_commutator_relative_residual"] for x in cases), default=0.0)
    max_wrong_fraction = max((x["wrong_JM_weight_fraction"] for x in cases), default=0.0)
    passed = (
        len(cases) == 4 * len(channels) * len(nodes)
        and max_commutator < 1e-8
        and max_adjoint < 1e-8
        and all(x["matrix_frobenius_norm"] > 1e-14 for x in adjoint_rows)
        and max_jz_commutator < 1e-12
        and max_wrong_fraction < 1e-18
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
        "magnetic_components_tested_M2": [-2, 0, 2],
        "adjointness_compression_rows": adjoint_rows,
        "max_adjointness_compression_relative_residual": max_adjoint,
        "pruning_mode": pruning_mode,
        "commutator_target": "[J_+, H_E]=[J_-, H_E]=[J_z, H_E]=0 on the actual complete W5 sine-H_E map",
        "max_commutator_relative_residual": max_commutator,
        "max_jz_commutator_relative_residual": max_jz_commutator,
        "max_wrong_JM_weight_fraction": max_wrong_fraction,
        "generators_checked": ["J_plus", "J_minus", "J_z"],
        "cases": cases,
        "scope": (
            "Direct all-M covariance check on the declared finite charged J=1 recoupling channels at each requested K5 source node. "
            "The adjointness result is only the compression to the declared seed subspace; it is not a global self-adjointness proof. "
            "It is not a proof for every W5 seed, charged channel, cutoff, graph refinement, or the physical continuum."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--pruning-mode", choices=("production", "shadow_no_prune"), default="production")
    parser.add_argument("--nodes", type=int, nargs="+", default=[0])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.Jmax2, args.pruning_mode, args.nodes)
    text = json.dumps(result, indent=2)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
