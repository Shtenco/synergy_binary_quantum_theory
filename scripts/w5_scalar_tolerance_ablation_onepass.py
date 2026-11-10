#!/usr/bin/env python3
"""One-pass W5 pruning/aggregation ablation on the exact captured seed.

Every T-sequence branch is generated once. Its complete-charge projection is
evaluated at the production and zero projection tolerances, then fed to the
production/zero accumulator variants in parallel. This preserves the same
operator, sequence order, coefficients, source key and internal basis cutoff
while avoiding seven repeated expensive branch-generation passes.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_zeroaware_volume_migration_experiment as ZVM
import peter_weyl_covariant_K_leg_gate as CK
import k5_peter_weyl_safe_hda_column as PW


PROJECTION_PRODUCTION = 1e-11
ACCUMULATOR_PRODUCTION = 1e-11


def make_key():
    return (
        (4,) * 10,
        (CK.label_charge(0, 0, 0, 0),) + (CK.label_gauss(0),) * 4,
    )


def summary(state):
    by_jm = defaultdict(list)
    for (_, labels), amp in state.items():
        lab = labels[0]
        if lab[0] != "C":
            raise RuntimeError(f"source node lost charged label: {lab!r}")
        by_jm[(int(lab[1]), int(lab[2]))].append(abs(amp) ** 2)
    weights = {
        f"J2={j},M2={m}": math.fsum(values)
        for (j, m), values in by_jm.items()
    }
    good = math.fsum(by_jm.get((0, 0), []))
    bad = math.fsum(
        value
        for (j, m), values in by_jm.items()
        if (j, m) != (0, 0)
        for value in values
    )
    total = math.fsum(abs(value) ** 2 for value in state.values())
    return {
        "support": len(state),
        "norm": math.sqrt(total),
        "allowed_J0_M0_weight": good,
        "wrong_sector_weight": bad,
        "wrong_sector_amplitude_norm": math.sqrt(bad),
        "wrong_sector_by_JM": {
            key: value for key, value in weights.items()
            if key != "J2=0,M2=0"
        },
    }


def map_difference(a, b):
    keys = set(a) | set(b)
    d2 = math.fsum(abs(a.get(key, 0j) - b.get(key, 0j)) ** 2 for key in keys)
    b2 = math.fsum(abs(value) ** 2 for value in b.values())
    return {
        "relative_L2": math.sqrt(d2) / max(math.sqrt(b2), 1e-300),
        "support_intersection": len(set(a) & set(b)),
        "support_symmetric_difference": len(set(a) ^ set(b)),
    }


def add(dst, src, scale, tol):
    CK.add_state(dst, src, scale, tol=tol)


def run_one_pass(key):
    base = CK.branch_from_key(key)
    specs = list(PW.oriented_specs(0))

    flat_p1_a1 = {}
    flat_p0_a1 = {}
    flat_p1_a0 = {}
    flat_p0_a0 = {}

    # Preserve per-spec maps so grouped and reversed-grouped accumulation can
    # be evaluated after the same single sequence-generation pass.
    grouped_spec_p1_a1 = []
    grouped_spec_p0_a0 = []

    # Exact same global term reference as the prior run_fsum implementation.
    fsum_terms = defaultdict(list)

    max_volume_leakage = 0.0
    max_basis_leakage_p1 = 0.0
    max_basis_leakage_p0 = 0.0
    branch_count = 0
    sequence_count = 0

    for sign, spec in specs:
        spec_p1_a1 = {}
        spec_p0_a0 = {}

        for adj in (False, True):
            prefactor = (+0.5j if adj else -0.5j) * sign
            for coefficient, sequence0 in PW.T_sequences(*spec):
                sequence = PW.adjoint_sequence(sequence0) if adj else sequence0
                branches, volume_leakage = CK.apply_sequence_to_branch(
                    base, sequence, 0, 10
                )
                max_volume_leakage = max(
                    max_volume_leakage, float(volume_leakage)
                )
                sequence_count += 1

                for branch in branches:
                    branch_count += 1

                    # Evaluate both projection tolerances on the identical
                    # generated branch. The internal 1e-13 basis cutoff is
                    # untouched in both calls.
                    part_p1, leakage_p1 = CK.project_branch_complete_charges(
                        branch, (0,), tol=PROJECTION_PRODUCTION
                    )
                    part_p0, leakage_p0 = CK.project_branch_complete_charges(
                        branch, (0,), tol=0.0
                    )
                    max_basis_leakage_p1 = max(
                        max_basis_leakage_p1, float(leakage_p1)
                    )
                    max_basis_leakage_p0 = max(
                        max_basis_leakage_p0, float(leakage_p0)
                    )
                    scale = prefactor * coefficient

                    # Four flat combinations isolate projection and
                    # accumulator pruning independently.
                    add(flat_p1_a1, part_p1, scale, ACCUMULATOR_PRODUCTION)
                    add(flat_p0_a1, part_p0, scale, ACCUMULATOR_PRODUCTION)
                    add(flat_p1_a0, part_p1, scale, 0.0)
                    add(flat_p0_a0, part_p0, scale, 0.0)

                    # The grouped variants use exactly the same within-spec
                    # order and are combined only after the spec is complete.
                    add(spec_p1_a1, part_p1, scale, ACCUMULATOR_PRODUCTION)
                    add(spec_p0_a0, part_p0, scale, 0.0)

                    for state_key, amplitude in part_p0.items():
                        fsum_terms[state_key].append(scale * amplitude)

        grouped_spec_p1_a1.append(spec_p1_a1)
        grouped_spec_p0_a0.append(spec_p0_a0)

    grouped_p1_a1 = {}
    for spec_state in grouped_spec_p1_a1:
        add(grouped_p1_a1, spec_state, 1.0, ACCUMULATOR_PRODUCTION)

    grouped_p0_a0 = {}
    for spec_state in grouped_spec_p0_a0:
        add(grouped_p0_a0, spec_state, 1.0, 0.0)

    grouped_reverse_p0_a0 = {}
    for spec_state in reversed(grouped_spec_p0_a0):
        add(grouped_reverse_p0_a0, spec_state, 1.0, 0.0)

    fsum_reference = {
        state_key: complex(
            math.fsum(value.real for value in values),
            math.fsum(value.imag for value in values),
        )
        for state_key, values in fsum_terms.items()
    }
    fsum_reference = {
        key: value for key, value in fsum_reference.items() if value != 0j
    }

    states = {
        "flat_p1e-11_a1e-11": flat_p1_a1,
        "flat_p0_a1e-11": flat_p0_a1,
        "flat_p1e-11_a0": flat_p1_a0,
        "flat_p0_a0": flat_p0_a0,
        "grouped_p1e-11_a1e-11": grouped_p1_a1,
        "grouped_p0_a0": grouped_p0_a0,
        "grouped_reverse_p0_a0": grouped_reverse_p0_a0,
        "flat_p0_fsum": fsum_reference,
    }

    diagnostics = {
        name: {
            "summary": summary(state),
            "max_volume_leakage": max_volume_leakage,
            "max_basis_leakage": (
                max_basis_leakage_p1
                if name in {"flat_p1e-11_a1e-11", "flat_p1e-11_a0",
                            "grouped_p1e-11_a1e-11"}
                else max_basis_leakage_p0
            ),
        }
        for name, state in states.items()
    }

    baseline = states["flat_p1e-11_a1e-11"]
    for name, state in states.items():
        diagnostics[name]["relative_L2_to_unpruned_fsum"] = map_difference(
            state, fsum_reference
        )
        diagnostics[name]["relative_L2_to_production_baseline"] = map_difference(
            state, baseline
        )

    return {
        "status": "W5 scalar tolerance/aggregation ablation; diagnostic only",
        "input_key": repr(key),
        "source_node": 0,
        "Jmax2": 10,
        "operator_changed": False,
        "guard_disabled": False,
        "branch_generation_passes": 1,
        "generated_sequence_count": sequence_count,
        "generated_branch_count": branch_count,
        "operator_path": (
            "unchanged complete T/T-dagger sequence path; identical oriented "
            "specs, adjoints, coefficients, branches and source key"
        ),
        "projection_tolerance": (
            "production 1e-11 versus zero, applied to each identical generated branch"
        ),
        "accumulator_tolerance": "production 1e-11 versus zero",
        "remaining_internal_basis_coefficient_cutoff": 1e-13,
        "aggregation_modes": [
            "flat interleaved branches",
            "one accumulator per oriented spec",
            "reversed oriented-spec order",
            "math.fsum real/imaginary reference",
        ],
        "variants": diagnostics,
        "interpretation": (
            "A lower forbidden-sector residual is insufficient by itself. "
            "Compare allowed-sector weight, full-map L2 residual, support, "
            "and volume/basis leakage. The production scalar-preservation "
            "guard remains untouched."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    ZVM.patch_and_clear()
    key = make_key()
    expected = (
        (4,) * 10,
        (("C", 0, 0, 0, 0),) + (("G", 0),) * 4,
    )
    if key != expected:
        raise RuntimeError(f"input key mismatch: {key!r}")

    result = run_one_pass(key)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
