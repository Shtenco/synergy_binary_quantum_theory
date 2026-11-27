#!/usr/bin/env python3
"""Sine-ordered [V,H_E] commutator column and complete reverse-row audit.

This is an audit-only finite-regulator runner. It computes one seed column
K|seed> with K=[V,H_sine], then computes K|b> for every b in the support of
that seed column and checks the complete reverse row <seed|K|b> against
-<b|K|seed>* . It is a full reverse-row check for one seed, not a full matrix
anti-Hermiticity proof on the entire reachable space.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_composition_gate as COMP
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def add(dst, src, scale=1.0):
    CK.add_state(dst, src, scale=scale, tol=0.0)


def norm(state):
    return math.sqrt(math.fsum(abs(z) ** 2 for z in state.values()))


def apply_h_sine_key(key, node, jmax2):
    """Apply (T-T†)/(2i) with no outer projection/accumulation pruning."""
    base = CK.branch_from_key(key)
    out = {}
    max_volume_leak = 0.0
    max_basis_leak = 0.0
    for sign, spec in PW.oriented_specs(node):
        v, a, b, c = spec
        for adj in (False, True):
            # (T-T†)/(2i) = -i*T/2 + i*T†/2
            prefactor = (-0.5j if not adj else +0.5j) * sign
            for coefficient, seq0 in PW.T_sequences(v, a, b, c):
                seq = PW.adjoint_sequence(seq0) if adj else seq0
                branches, volume_leak = CK.apply_sequence_to_branch(
                    base, seq, node, jmax2
                )
                max_volume_leak = max(max_volume_leak, float(volume_leak))
                for branch in branches:
                    projected, basis_leak = CK.project_branch_complete_charges(
                        branch, (node,), tol=0.0
                    )
                    max_basis_leak = max(max_basis_leak, float(basis_leak))
                    add(out, projected, prefactor * coefficient)
    return out, max_volume_leak, max_basis_leak


def apply_h_sine_state(state, node, jmax2):
    out = {}
    max_volume_leak = 0.0
    max_basis_leak = 0.0
    for key, amplitude in state.items():
        column, volume_leak, basis_leak = apply_h_sine_key(key, node, jmax2)
        max_volume_leak = max(max_volume_leak, volume_leak)
        max_basis_leak = max(max_basis_leak, basis_leak)
        add(out, column, amplitude)
    return out, max_volume_leak, max_basis_leak


def apply_v_state(state, node):
    # Uses the established complete-charge volume block; its internal
    # coefficient/output cutoffs are recorded explicitly in the report.
    return CK.apply_volume_complete_state(state, node)


def apply_k_column(key, node, jmax2):
    seed = {key: 1.0 + 0.0j}
    h_seed, hv1, hb1 = apply_h_sine_state(seed, node, jmax2)
    vh = apply_v_state(h_seed, node)
    v_seed = apply_v_state(seed, node)
    hv, hv2, hb2 = apply_h_sine_state(v_seed, node, jmax2)
    out = {}
    add(out, vh, +1.0)
    add(out, hv, -1.0)
    return out, {
        "H_sine_volume_leakage": max(hv1, hv2),
        "H_sine_basis_leakage": max(hb1, hb2),
    }


def make_seed(node, input_k):
    gauss_key = (
        (4,) * 10,
        tuple(input_k if v == node else 0 for v in range(5)),
    )
    covariant = CV.gauss_to_covariant({gauss_key: 1.0 + 0.0j}, node)
    complete = COMP.covariant_to_complete(covariant, node)
    if len(complete) != 1:
        raise RuntimeError(f"expected one complete seed, got {len(complete)}")
    key, amplitude = next(iter(complete.items()))
    if abs(amplitude - 1.0) > 1e-12:
        raise RuntimeError(f"seed normalization defect: {amplitude!r}")
    return key


def run(node, input_k, jmax2, max_reverse_columns):
    ZVM.patch_and_clear()
    seed_key = make_seed(node, input_k)
    k_seed, leakage = apply_k_column(seed_key, node, jmax2)
    reverse_basis = dict(k_seed)
    reverse_basis.setdefault(seed_key, 0.0j)  # include the diagonal matrix element
    ranked = sorted(reverse_basis.items(), key=lambda item: abs(item[1]), reverse=True)
    selected = ranked if max_reverse_columns <= 0 else ranked[:max_reverse_columns]

    rows = []
    max_abs_defect = 0.0
    max_relative_defect = 0.0
    for b_key, k_ba in selected:
        k_on_b, leak_b = apply_k_column(b_key, node, jmax2)
        k_ab = complex(k_on_b.get(seed_key, 0.0j))
        defect = abs(k_ba + k_ab.conjugate())
        scale = max(abs(k_ba), abs(k_ab), 1e-30)
        relative = defect / scale
        max_abs_defect = max(max_abs_defect, defect)
        max_relative_defect = max(max_relative_defect, relative)
        rows.append({
            "abs_K_ba": abs(k_ba),
            "abs_K_ab": abs(k_ab),
            "antihermitian_absolute_defect": defect,
            "antihermitian_relative_defect": relative,
            "reverse_column_support": len(k_on_b),
            "reverse_column_leakage": leak_b,
        })

    return {
        "status": "sine-ordered finite-regulator K commutator seed-column/reverse-row audit",
        "passed": bool(
            len(k_seed) > 0
            and len(rows) == len(ranked)
            and max_relative_defect < 5e-7
        ),
        "operator": "K_comm_sine=[V,H_E_sine], H_E_sine=(T-Tdagger)/(2i)",
        "node": node,
        "input_K": input_k,
        "Jmax2": jmax2,
        "seed_support": 1,
        "K_seed_support": len(k_seed),
        "reverse_basis_support_including_seed": len(reverse_basis),
        "K_seed_norm": norm(k_seed),
        "reverse_columns_requested": len(ranked),
        "reverse_columns_evaluated": len(rows),
        "reverse_row_coverage": len(rows) / max(len(ranked), 1),
        "max_antihermitian_absolute_defect": max_abs_defect,
        "max_antihermitian_relative_defect": max_relative_defect,
        "seed_column_leakage": leakage,
        "reverse_rows": rows,
        "cutoff_ledger": {
            "outer_complete_charge_projection": 0.0,
            "outer_state_accumulation": 0.0,
            "charged/Gauss recoupling coefficient_cutoff": 1e-13,
            "complete-volume matrix coefficient cutoff": 1e-13,
            "complete-volume output cutoff": 1e-11,
            "path-branch amplitude cutoff": 0.0,
            "path-branch regulator": "finite doubled-spin ceiling Jmax2; no amplitude pruning in apply_path_branch/apply_hit_branch",
            "floating_dtype": "complex128 / float64 numerical backend",
        },
        "scope": (
            "Checks the entire reverse row associated with one seed column, "
            "not every pair in the full reachable-space matrix. Internal "
            "operator cutoffs remain active. No continuum or physical-HDA claim."
        ),
        "two_leg_conversion": (
            "If this finite commutator is established anti-Hermitian on the "
            "same domain, K_H=(i*hbar)^-1 K_comm; two such legs contribute "
            "-hbar^-2. This does not fix the full Lorentzian Hamiltonian prefactor."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node", type=int, default=0, choices=range(5))
    parser.add_argument("--input-K", type=int, default=0, choices=(0, 2, 4, 6, 8))
    parser.add_argument("--Jmax2", type=int, default=10)
    parser.add_argument("--max-reverse-columns", type=int, default=0,
                        help="0 evaluates every state in the seed-column support")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.Jmax2 < 10 or args.Jmax2 % 2:
        raise ValueError("Jmax2 must be an even doubled-spin cutoff >= 10")
    result = run(args.node, args.input_K, args.Jmax2, args.max_reverse_columns)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
