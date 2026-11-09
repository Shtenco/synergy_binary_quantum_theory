#!/usr/bin/env python3
"""Falsification-first W5 scalar-covariance diagnostic across Peter-Weyl cutoffs.

This runner deliberately stops before KKV, five-node closure, physical ordering,
Gram assembly, coefficient fitting, or any joint physical-kernel calculation.
It evaluates one W5 seed at one source node, decomposes each oriented Euclidean
node contribution into T and T^dagger, and compares the legacy plus and sine
Hermitian combinations on identical inputs and cutoffs.
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

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_composition_gate as COMP


def add(dst, src, scale=1.0, tol=0.0):
    CK.add_state(dst, src, scale, tol)


def norm(state):
    return math.sqrt(sum(abs(z) ** 2 for z in state.values()))


def weights_jm(state, node):
    out = {}
    for (_, labels), amp in state.items():
        lab = labels[node]
        if lab[0] != "C":
            raise RuntimeError(f"node {node} lost charged label: {lab!r}")
        _, j2, m2, _, _ = lab
        k = (int(j2), int(m2))
        out[k] = out.get(k, 0.0) + abs(amp) ** 2
    return {
        f"J={j2/2:g},M={m2/2:g}": float(value)
        for (j2, m2), value in sorted(out.items())
    }


def sector_metrics(state, node, expected_j2=0, expected_m2=0):
    w = weights_jm(state, node)
    total = sum(w.values())
    good = w.get(f"J={expected_j2/2:g},M={expected_m2/2:g}", 0.0)
    bad = max(0.0, total - good)
    return {
        "weight_by_JM": w,
        "total_weight": float(total),
        "expected_sector_weight": float(good),
        "wrong_sector_weight": float(bad),
        "wrong_sector_fraction": float(bad / max(total, 1e-300)),
    }


def support_counts(state):
    return {
        f"amp_gt_{tol:g}": sum(abs(z) > tol for z in state.values())
        for tol in (1e-13, 1e-11, 1e-9)
    }


def compare_maps(a, b):
    keys = set(a) | set(b)
    diff2 = sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys)
    b2 = sum(abs(z) ** 2 for z in b.values())
    sa, sb = set(a), set(b)
    return {
        "relative_l2_difference": math.sqrt(diff2) / max(math.sqrt(b2), 1e-300),
        "support_a": len(sa),
        "support_b": len(sb),
        "common_support": len(sa & sb),
        "a_only_support": len(sa - sb),
        "b_only_support": len(sb - sa),
    }


def evaluate(cutoff, node, input_k):
    import peter_weyl_zeroaware_volume_migration_experiment as ZVM
    ZVM.patch_and_clear()

    gauss_key = ((4,) * 10, tuple(input_k if v == node else 0 for v in range(5)))
    cov = CV.gauss_to_covariant({gauss_key: 1.0 + 0j}, node)
    complete = COMP.covariant_to_complete(cov, node)
    if len(complete) != 1:
        raise RuntimeError(f"expected one covariant seed, got {len(complete)}")
    input_key, input_amp = next(iter(complete.items()))
    if abs(input_amp - 1.0) > 1e-12:
        raise RuntimeError(f"Gauss-to-covariant map is not unit-normalized: {input_amp}")

    base = CK.branch_from_key(input_key)
    total_t = {}
    total_tdag = {}
    total_plus = {}
    total_sine = {}
    per_spec = []
    max_volume_leak = 0.0
    max_basis_leak = 0.0

    for sign, spec in PW.oriented_specs(node):
        spec_t = {}
        spec_tdag = {}
        spec_volume_leak = 0.0
        spec_basis_leak = 0.0

        for coef, seq0 in PW.T_sequences(*spec):
            branches, vleak = CK.apply_sequence_to_branch(
                base, seq0, node, cutoff
            )
            spec_volume_leak = max(spec_volume_leak, float(vleak))
            for branch in branches:
                projected, bleak = CK.project_branch_complete_charges(
                    branch, (node,), tol=0.0
                )
                spec_basis_leak = max(spec_basis_leak, float(bleak))
                add(spec_t, projected, coef)

            adj_seq = PW.adjoint_sequence(seq0)
            adj_branches, adj_vleak = CK.apply_sequence_to_branch(
                base, adj_seq, node, cutoff
            )
            spec_volume_leak = max(spec_volume_leak, float(adj_vleak))
            for branch in adj_branches:
                projected, bleak = CK.project_branch_complete_charges(
                    branch, (node,), tol=0.0
                )
                spec_basis_leak = max(spec_basis_leak, float(bleak))
                add(spec_tdag, projected, coef)

        add(total_t, spec_t, sign)
        add(total_tdag, spec_tdag, sign)
        add(total_plus, spec_t, 0.5 * sign)
        add(total_plus, spec_tdag, 0.5 * sign)
        add(total_sine, spec_t, -0.5j * sign)
        add(total_sine, spec_tdag, +0.5j * sign)

        spec_plus = {}
        spec_sine = {}
        add(spec_plus, spec_t, 0.5)
        add(spec_plus, spec_tdag, 0.5)
        add(spec_sine, spec_t, -0.5j)
        add(spec_sine, spec_tdag, +0.5j)
        per_spec.append({
            "oriented_spec": list(spec),
            "epsilon_sign": int(sign),
            "T_support": len(spec_t),
            "T_norm": norm(spec_t),
            "T_sector": sector_metrics(spec_t, node),
            "Tdag_support": len(spec_tdag),
            "Tdag_norm": norm(spec_tdag),
            "Tdag_sector": sector_metrics(spec_tdag, node),
            "plus_support": len(spec_plus),
            "plus_norm": norm(spec_plus),
            "plus_sector": sector_metrics(spec_plus, node),
            "sine_support": len(spec_sine),
            "sine_norm": norm(spec_sine),
            "sine_sector": sector_metrics(spec_sine, node),
            "max_volume_leakage": spec_volume_leak,
            "max_complete_basis_leakage": spec_basis_leak,
        })
        max_volume_leak = max(max_volume_leak, spec_volume_leak)
        max_basis_leak = max(max_basis_leak, spec_basis_leak)

    return {
        "cutoff_doubled_spin": cutoff,
        "Jmax": cutoff / 2,
        "source_node": node,
        "input_K_labels": [input_k if v == node else 0 for v in range(5)],
        "input_spins": [4] * 10,
        "input_covariant_support": len(cov),
        "input_covariant_norm": norm(cov),
        "oriented_specs_completed": len(per_spec),
        "T_total": {
            "support": len(total_t), "norm": norm(total_t),
            "sector": sector_metrics(total_t, node),
            "support_thresholds": support_counts(total_t),
        },
        "Tdag_total": {
            "support": len(total_tdag), "norm": norm(total_tdag),
            "sector": sector_metrics(total_tdag, node),
            "support_thresholds": support_counts(total_tdag),
        },
        "H_plus": {
            "definition": "(T+T^dagger)/2",
            "support": len(total_plus), "norm": norm(total_plus),
            "sector": sector_metrics(total_plus, node),
            "support_thresholds": support_counts(total_plus),
        },
        "H_sine": {
            "definition": "(T-T^dagger)/(2i)",
            "support": len(total_sine), "norm": norm(total_sine),
            "sector": sector_metrics(total_sine, node),
            "support_thresholds": support_counts(total_sine),
        },
        "plus_vs_sine": compare_maps(total_plus, total_sine),
        "max_volume_leakage": max_volume_leak,
        "max_complete_basis_leakage": max_basis_leak,
        "per_oriented_spec": per_spec,
        "_maps": {
            "T_total": total_t,
            "Tdag_total": total_tdag,
            "H_plus": total_plus,
            "H_sine": total_sine,
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--node", type=int, default=0, choices=range(5))
    ap.add_argument("--input-K", type=int, default=0, choices=(0, 2, 4, 6, 8))
    ap.add_argument("--cutoffs", type=int, nargs="+", default=(10, 12, 14))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if any(c < 10 or c % 2 for c in args.cutoffs):
        raise ValueError("W5 cutoffs must be even doubled-spin walls >= 10")

    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def checkpoint(status):
        # Persist completed cutoffs immediately. A long W5 run must leave
        # recoverable evidence even if a later cutoff times out.
        partial = {
            "status": status,
            "physical_joint_kernel_computed": False,
            "five_node_closure_computed": False,
            "physical_coefficients_fitted": False,
            "node": args.node,
            "input_K": args.input_K,
            "cutoffs_requested": list(args.cutoffs),
            "cutoffs_completed": [row["cutoff_doubled_spin"] for row in rows],
            "cutoff_runs": [
                {k: v for k, v in row.items() if k != "_maps"}
                for row in rows
            ],
            "cutoff_trends": [],
            "interpretation_rule": (
                "A nonzero wrong-sector weight is diagnostic evidence, not "
                "grounds for projection or threshold tuning. Compare H_plus "
                "and H_sine and inspect oriented-spec and T/Tdag contributions."
            ),
        }
        args.output.write_text(json.dumps(partial, indent=2) + "\n", encoding="utf-8")

    checkpoint("RUNNING_W5_SCALAR_COVARIANCE_DIAGNOSTIC")
    for cutoff in args.cutoffs:
        row = evaluate(cutoff, args.node, args.input_K)
        rows.append(row)
        checkpoint("PARTIAL_W5_SCALAR_COVARIANCE_DIAGNOSTIC")
        print(json.dumps({
            "completed_cutoff": cutoff,
            "H_plus": row["H_plus"],
            "H_sine": row["H_sine"],
            "max_volume_leakage": row["max_volume_leakage"],
            "max_complete_basis_leakage": row["max_complete_basis_leakage"],
            "output": str(args.output),
        }, indent=2), flush=True)

    convergence = []
    for prev, curr in zip(rows, rows[1:]):
        for op in ("H_plus", "H_sine", "T_total", "Tdag_total"):
            convergence.append({
                "from_cutoff": prev["cutoff_doubled_spin"],
                "to_cutoff": curr["cutoff_doubled_spin"],
                "map": op,
                **compare_maps(prev["_maps"][op], curr["_maps"][op]),
                "norm_ratio": (
                    curr[op]["norm"] / max(prev[op]["norm"], 1e-300)
                ),
                "support_delta": curr[op]["support"] - prev[op]["support"],
                "wrong_sector_fraction_change": (
                    curr[op]["sector"]["wrong_sector_fraction"]
                    - prev[op]["sector"]["wrong_sector_fraction"]
                ),
            })

    result = {
        "status": "W5 scalar-covariance diagnostic; no KKV or physical kernel",
        "physical_joint_kernel_computed": False,
        "five_node_closure_computed": False,
        "physical_coefficients_fitted": False,
        "node": args.node,
        "input_K": args.input_K,
        "cutoff_runs": [
            {k: v for k, v in row.items() if k != "_maps"}
            for row in rows
        ],
        "cutoff_trends": convergence,
        "interpretation_rule": (
            "A nonzero wrong-sector weight is a hard scalar-covariance blocker. "
            "Do not remove it by projection or threshold tuning. Compare H_plus "
            "and H_sine on identical input/cutoff; inspect per-oriented-spec and "
            "T/Tdag contributions before assigning a unique cause."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "node": args.node,
        "input_K": args.input_K,
        "cutoff_summary": [
            {
                "Jmax2": row["cutoff_doubled_spin"],
                "H_plus_norm": row["H_plus"]["norm"],
                "H_plus_support": row["H_plus"]["support"],
                "H_plus_wrong_sector_fraction": row["H_plus"]["sector"]["wrong_sector_fraction"],
                "H_sine_norm": row["H_sine"]["norm"],
                "H_sine_support": row["H_sine"]["support"],
                "H_sine_wrong_sector_fraction": row["H_sine"]["sector"]["wrong_sector_fraction"],
                "max_volume_leakage": row["max_volume_leakage"],
                "max_basis_leakage": row["max_complete_basis_leakage"],
            }
            for row in rows
        ],
        "output": str(args.output),
        "physical_joint_kernel_computed": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
