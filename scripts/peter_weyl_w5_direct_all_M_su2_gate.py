#!/usr/bin/env python3
"""Direct all-M SU(2) commutator audit of the actual W5 sine-H_E operator.

This deliberately bypasses the cached highest-weight wrapper. It evaluates the
same production H_E map on every magnetic component of a fixed charged J=1
four-leg state and checks [J_+, H_E]=[J_-, H_E]=0 on the complete sparse output.
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


def norm(state):
    return math.sqrt(sum(abs(amp) ** 2 for amp in state.values()))


def relative_error(a, b):
    err2 = 0.0
    for key in set(a) | set(b):
        err2 += abs(a.get(key, 0j) - b.get(key, 0j)) ** 2
    scale = max(norm(a), norm(b), 1e-30)
    return math.sqrt(err2) / scale, math.sqrt(err2), scale


def make_input(M2):
    spins = (4,) * 10
    labels = (
        ("C", 2, int(M2), 0, 2),
        ("G", 0), ("G", 0), ("G", 0), ("G", 0),
    )
    return spins, labels


def run(Jmax2=10, pruning_mode="production"):
    if pruning_mode not in ("production", "shadow_no_prune"):
        raise ValueError("invalid pruning_mode")
    if Jmax2 < 4:
        raise ValueError("Jmax2 must be at least 4 for this fixed W5 seed")
    directions = (1, -1)
    cases = []
    context = PRUNE.bypass() if pruning_mode == "shadow_no_prune" else None

    def evaluate():
        for direction in directions:
            source_M2_values = (-2, 0) if direction > 0 else (0, 2)
            for M2 in source_M2_values:
                input_key = make_input(M2)
                raised_input = apply_local_ladder({input_key: 1.0 + 0j}, 0, direction)
                if len(raised_input) != 1:
                    raise AssertionError("ladder action should produce one adjacent M sector")
                raised_key, ladder_factor = next(iter(raised_input.items()))
                with PRUNE.capture() as collector:
                    H_input, v_leak, basis_leak = SINE.complete_HE_sine(
                        input_key, 0, Jmax2, charged_nodes=(0,)
                    )
                    H_raised, v_leak_2, basis_leak_2 = SINE.complete_HE_sine(
                        raised_key, 0, Jmax2, charged_nodes=(0,)
                    )
                lhs = apply_local_ladder(H_input, 0, direction)
                rhs = {key: amp * ladder_factor for key, amp in H_raised.items()}
                rel, abs_err, scale = relative_error(lhs, rhs)
                total_weight = sum(abs(amp) ** 2 for (_, labels), amp in H_input.items()
                                    if labels[0][0] == "C")
                wrong_weight = sum(
                    abs(amp) ** 2 for (_, labels), amp in H_input.items()
                    if labels[0][0] != "C" or labels[0][1] != 2 or labels[0][2] != M2
                )
                wrong_fraction = wrong_weight / total_weight if total_weight > 0.0 else 0.0
                cases.append({
                    "direction": "J_plus" if direction > 0 else "J_minus",
                    "input_M2": M2,
                    "output_M2": M2 + 2 * direction,
                    "ladder_factor": ladder_factor,
                    "input_HE_support": len(H_input),
                    "raised_HE_support": len(H_raised),
                    "commutator_relative_residual": rel,
                    "commutator_absolute_residual": abs_err,
                    "comparison_scale": scale,
                    "wrong_JM_weight_fraction": wrong_fraction,
                    "wrong_JM_relative_amplitude": math.sqrt(wrong_fraction),
                    "volume_leakage": max(float(v_leak), float(v_leak_2)),
                    "basis_leakage": max(float(basis_leak), float(basis_leak_2)),
                    "pruning_events": collector.snapshot()["event_count"],
                })

    if context is None:
        with PRUNE.cache_bypass():
            evaluate()
    else:
        with PRUNE.cache_bypass(), context:
            evaluate()

    max_commutator = max((x["commutator_relative_residual"] for x in cases), default=0.0)
    max_wrong_fraction = max((x["wrong_JM_weight_fraction"] for x in cases), default=0.0)
    passed = (
        len(cases) == 4
        and max_commutator < 1e-8
        and max_wrong_fraction < 1e-18
        and max((x["volume_leakage"] for x in cases), default=0.0) < 1e-10
        and max((x["basis_leakage"] for x in cases), default=0.0) < 1e-10
    )
    return {
        "schema": "bqg-w5-direct-all-M-su2-commutator-v1",
        "status": "PASS_DIRECT_ALL_M_SU2_COMMUTATOR" if passed else "FAIL_DIRECT_ALL_M_SU2_COMMUTATOR",
        "passed": bool(passed),
        "Jmax2": Jmax2,
        "source_node": 0,
        "input_spins_doubled": [4] * 10,
        "input_charged_sector": {"J2": 2, "K12": 0, "K34": 2},
        "pruning_mode": pruning_mode,
        "commutator_target": "[J_+, H_E]=[J_-, H_E]=0 on the actual complete W5 sine-H_E map",
        "max_commutator_relative_residual": max_commutator,
        "max_wrong_JM_weight_fraction": max_wrong_fraction,
        "cases": cases,
        "scope": (
            "Direct all-M covariance check on one finite charged J=1 seed and one K5 node. "
            "It is not a proof for every W5 seed, cutoff, graph refinement, or the physical continuum."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--pruning-mode", choices=("production", "shadow_no_prune"), default="production")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.Jmax2, args.pruning_mode)
    text = json.dumps(result, indent=2)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
