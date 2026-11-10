#!/usr/bin/env python3
"""Isolate W5 scalar leakage from projection/accumulation tolerances.

Uses the exact captured K=0, node-0, Jmax2=10 canonical input. The operator,
branch generation, sequences, coefficients and charge basis are unchanged.
Only the two numerical pruning tolerances and accumulation order are varied.
The scalar-restoration guard is not disabled.
"""
from __future__ import annotations
import argparse, json, math, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_zeroaware_volume_migration_experiment as ZVM
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_sine_composition_gate as SINEK
import k5_peter_weyl_safe_hda_column as PW


def make_key():
    return ((4,) * 10,
            (CK.label_charge(0, 0, 0, 0),) + (CK.label_gauss(0),) * 4)


def summary(state):
    by_jm = defaultdict(list)
    for (_, labels), amp in state.items():
        lab = labels[0]
        if lab[0] != "C":
            raise RuntimeError(f"source node lost charged label: {lab!r}")
        by_jm[(int(lab[1]), int(lab[2]))].append(abs(amp) ** 2)
    weights = {f"J2={j},M2={m}": math.fsum(v)
               for (j, m), v in by_jm.items()}
    good = math.fsum(by_jm.get((0, 0), []))
    bad = math.fsum(v for (j, m), vals in by_jm.items()
                    if (j, m) != (0, 0) for v in vals)
    total = math.fsum(abs(z) ** 2 for z in state.values())
    return {
        "support": len(state), "norm": math.sqrt(total),
        "allowed_J0_M0_weight": good, "wrong_sector_weight": bad,
        "wrong_sector_amplitude_norm": math.sqrt(bad),
        "wrong_sector_by_JM": {k: v for k, v in weights.items()
                               if k != "J2=0,M2=0"},
    }


def map_difference(a, b):
    keys = set(a) | set(b)
    d2 = math.fsum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys)
    b2 = math.fsum(abs(z) ** 2 for z in b.values())
    return {"relative_L2": math.sqrt(d2) / max(math.sqrt(b2), 1e-300),
            "support_intersection": len(set(a) & set(b)),
            "support_symmetric_difference": len(set(a) ^ set(b))}


def patched_tolerances(project_tol, add_tol):
    old_project, old_add = CK.project_branch_complete_charges, CK.add_state
    def project(branch, charged_nodes=(0, 1), tol=1e-11):
        return old_project(branch, charged_nodes, tol=project_tol)
    def add(dst, src, scale=1.0, tol=1e-11):
        return old_add(dst, src, scale, tol=add_tol)
    CK.project_branch_complete_charges, CK.add_state = project, add
    return old_project, old_add


def run_flat(key, project_tol, add_tol):
    old_project, old_add = patched_tolerances(project_tol, add_tol)
    try:
        return SINEK.complete_HE_sine(key, 0, 10, charged_nodes=(0,))
    finally:
        CK.project_branch_complete_charges, CK.add_state = old_project, old_add


def run_grouped(key, project_tol, add_tol, reverse=False):
    old_project, old_add = patched_tolerances(project_tol, add_tol)
    try:
        base, out = CK.branch_from_key(key), {}
        specs = list(PW.oriented_specs(0))
        if reverse:
            specs.reverse()
        max_v = max_b = 0.0
        for sign, spec in specs:
            spec_out = {}
            for adj in (False, True):
                pref = (+0.5j if adj else -0.5j) * sign
                for coef, seq0 in PW.T_sequences(*spec):
                    seq = PW.adjoint_sequence(seq0) if adj else seq0
                    branches, vl = CK.apply_sequence_to_branch(base, seq, 0, 10)
                    max_v = max(max_v, float(vl))
                    for branch in branches:
                        part, bl = CK.project_branch_complete_charges(branch, (0,))
                        max_b = max(max_b, float(bl))
                        CK.add_state(spec_out, part, pref * coef)
            CK.add_state(out, spec_out)
        return out, max_v, max_b
    finally:
        CK.project_branch_complete_charges, CK.add_state = old_project, old_add


def run_fsum(key, project_tol=0.0):
    old_project = CK.project_branch_complete_charges
    def project(branch, charged_nodes=(0, 1), tol=1e-11):
        return old_project(branch, charged_nodes, tol=project_tol)
    CK.project_branch_complete_charges = project
    try:
        base, terms = CK.branch_from_key(key), defaultdict(list)
        max_v = max_b = 0.0
        for sign, spec in PW.oriented_specs(0):
            for adj in (False, True):
                pref = (+0.5j if adj else -0.5j) * sign
                for coef, seq0 in PW.T_sequences(*spec):
                    seq = PW.adjoint_sequence(seq0) if adj else seq0
                    branches, vl = CK.apply_sequence_to_branch(base, seq, 0, 10)
                    max_v = max(max_v, float(vl))
                    for branch in branches:
                        part, bl = CK.project_branch_complete_charges(branch, (0,))
                        max_b = max(max_b, float(bl))
                        scale = pref * coef
                        for k, amp in part.items():
                            terms[k].append(scale * amp)
        out = {k: complex(math.fsum(z.real for z in vals),
                           math.fsum(z.imag for z in vals))
               for k, vals in terms.items()}
        return {k: z for k, z in out.items() if z != 0j}, max_v, max_b
    finally:
        CK.project_branch_complete_charges = old_project


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    ZVM.patch_and_clear()
    key = make_key()
    expected = ((4,) * 10, (("C", 0, 0, 0, 0),) + (("G", 0),) * 4)
    if key != expected:
        raise RuntimeError(f"input key mismatch: {key!r}")

    configs = [
        ("flat_p1e-11_a1e-11", 1e-11, 1e-11, "flat", False),
        ("flat_p0_a1e-11", 0.0, 1e-11, "flat", False),
        ("flat_p1e-11_a0", 1e-11, 0.0, "flat", False),
        ("flat_p0_a0", 0.0, 0.0, "flat", False),
        ("grouped_p1e-11_a1e-11", 1e-11, 1e-11, "grouped", False),
        ("grouped_p0_a0", 0.0, 0.0, "grouped", False),
        ("grouped_reverse_p0_a0", 0.0, 0.0, "grouped", True),
    ]
    states = {}
    diagnostics = {}
    for name, pt, at, mode, reverse in configs:
        if mode == "flat":
            state, vl, bl = run_flat(key, pt, at)
        else:
            state, vl, bl = run_grouped(key, pt, at, reverse)
        states[name] = state
        diagnostics[name] = {"summary": summary(state),
                             "max_volume_leakage": vl,
                             "max_basis_leakage": bl}
    ref, vl, bl = run_fsum(key, 0.0)
    states["flat_p0_fsum"] = ref
    diagnostics["flat_p0_fsum"] = {
        "summary": summary(ref), "max_volume_leakage": vl,
        "max_basis_leakage": bl,
    }
    baseline = states["flat_p1e-11_a1e-11"]
    for name in diagnostics:
        diagnostics[name]["relative_L2_to_unpruned_fsum"] = map_difference(states[name], ref)
        diagnostics[name]["relative_L2_to_production_baseline"] = map_difference(states[name], baseline)

    result = {
        "status": "W5 scalar tolerance/aggregation ablation; diagnostic only",
        "input_key": repr(key), "source_node": 0, "Jmax2": 10,
        "operator_changed": False, "guard_disabled": False,
        "operator_path": "unchanged complete_HE_sine; identical oriented specs, sequences, adjoints and branch coefficients",
        "projection_tolerance": "tol passed to project_branch_complete_charges (production 1e-11, unpruned 0)",
        "accumulator_tolerance": "tol passed to add_state (production 1e-11, unpruned 0)",
        "remaining_internal_basis_coefficient_cutoff": 1e-13,
        "aggregation_modes": ["flat interleaved branches", "one accumulator per oriented spec", "reversed oriented-spec order", "math.fsum real/imag reference"],
        "variants": diagnostics,
        "interpretation": "A lower forbidden-sector residual is insufficient by itself. Compare allowed-sector weight, full-map L2 residual, support, and volume/basis leakage. The production scalar-preservation guard remains untouched.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
