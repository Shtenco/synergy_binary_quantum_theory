#!/usr/bin/env python3
"""Deterministically aggregate execution shards of the preregistered H_L witness.

In addition to the scientific logical-return diagnostics, the aggregate now
preserves the complete final scalar Gauss-basis outgoing state.  This is not a
new observable: it prevents loss of the expensive Peter-Weyl column so that the
same calculation can feed the preregistered Lorentzian Gram/master assembly
without recomputation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL


def decode_state(rows):
    out = {}
    for row in rows:
        key = (
            tuple(int(x) for x in row["spins"]),
            tuple(int(x) for x in row["Kother"]),
            int(row["J2"]),
            int(row["M2"]),
            int(row["K12"]),
            int(row["K34"]),
        )
        amp = complex(float(row["amp"][0]), float(row["amp"][1]))
        out[key] = out.get(key, 0j) + amp
    return out


def encode_gauss_state(state):
    rows = []
    for key, amp in sorted(state.items(), key=lambda kv: repr(kv[0])):
        spins, Ks = key
        rows.append({
            "spins": [int(x) for x in spins],
            "K_labels": [int(x) for x in Ks],
            "amp": [float(complex(amp).real), float(complex(amp).imag)],
        })
    return rows


def run(input_dir: Path, shards=8, source_v=0, input_index=0):
    paths = sorted(input_dir.glob("shard_*.json"))
    by_id = {}
    for path in paths:
        row = json.loads(path.read_text(encoding="utf-8"))
        by_id[int(row["shard"])] = row
    missing = [i for i in range(shards) if i not in by_id]
    if missing:
        raise RuntimeError(f"missing shards: {missing}")

    total = {}
    term_rows = []
    term_ids = []
    global_diag = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    max_spin = 0.0
    shard_pass = True

    for sid in range(shards):
        row = by_id[sid]
        if int(row["shards"]) != shards:
            raise RuntimeError(f"shard-count mismatch in shard {sid}")
        if int(row["source_node"]) != source_v or int(row["input_logical_basis_index"]) != input_index:
            raise RuntimeError(f"input mismatch in shard {sid}")
        shard_pass = shard_pass and bool(row["passed"])
        state = decode_state(row["signed_partial_state"])
        FULL.add(total, state)
        term_rows.extend(row["terms"])
        term_ids.extend(int(x) for x in row["assigned_term_ids"])
        FULL.update_diag(global_diag, row["max_diagnostics"])
        max_spin = max(max_spin, float(row["max_spin_reached"]))

    sorted_ids = sorted(term_ids)
    exact_terms = sorted_ids == list(range(24)) and len(term_ids) == 24 and len(set(term_ids)) == 24

    full_scalar = FULL.scalar_diagnostics(total)
    gauss, reverse = FULL.project_covariant_J0_to_gauss(total, source_v)
    logical = FULL.logical_projection(gauss)
    logical_norm = math.sqrt(FULL.norm2(logical))
    full_norm = math.sqrt(FULL.norm2(total))
    gauss_norm = math.sqrt(FULL.norm2(gauss))
    logical_nonzero = logical_norm > FULL.NONZERO_TOL

    basis = FULL.RAW.PW.basis_full_jhalf()
    initial = basis[input_index]
    initial_return = complex(logical.get(initial, 0j))
    ranked = sorted(logical.items(), key=lambda kv: abs(kv[1]), reverse=True)
    logical_rows = [
        {
            "K_labels": list(key[1]),
            "abs_amp": float(abs(amp)),
            "amp": [float(amp.real), float(amp.imag)],
        }
        for key, amp in ranked
    ]

    max_term_frac = 0.0
    max_term_nearzero_nonscalar = 0.0
    all_terms_pass = True
    for t in term_rows:
        all_terms_pass = all_terms_pass and bool(t["passed"])
        sd = t["scalar_diagnostics"]
        if float(sd["norm"]) > FULL.NONZERO_TOL:
            max_term_frac = max(max_term_frac, float(sd["nonscalar_weight_fraction"]))
        else:
            max_term_nearzero_nonscalar = max(max_term_nearzero_nonscalar, float(sd["nonscalar_norm"]))

    hard = {
        "all_shards_passed": bool(shard_pass),
        "exactly_the_24_preregistered_terms": bool(exact_terms),
        "all_term_integrity_checks_passed": bool(all_terms_pass),
        "all_ordered_terms_scalar_within_frozen_threshold": bool(
            max_term_frac < 1e-8 and max_term_nearzero_nonscalar < FULL.NONZERO_TOL
        ),
        "full_epsilon_sum_scalar_within_frozen_threshold": bool(FULL.scalar_ok(full_scalar)),
        "CV_complete_basis_leakage_below_1e-9": bool(global_diag["CV_complete_basis_leakage"] < 1e-9),
        "CK_outer_complete_basis_leakage_below_1e-9": bool(global_diag["CK_outer_complete_basis_leakage"] < 1e-9),
        "CK_internal_volume_sector_leakage_below_1e-9": bool(global_diag["CK_internal_volume_sector_leakage"] < 1e-9),
        "spin_cutoff_respected": bool(max_spin <= FULL.JMAX2 / 2 + 1e-12),
        "J0_covariant_to_gauss_mapping_valid": bool(
            not reverse["invalid_J0_covariant_keys"]
            and reverse["mapping_collisions"] == 0
            and reverse["mapped_covariant_J0_basis_states"] == reverse["distinct_gauss_basis_states"]
        ),
        "scalar_covariant_norm_equals_gauss_norm": bool(abs(full_norm - gauss_norm) < 1e-8 * max(1.0, full_norm)),
    }
    passed = bool(all(hard.values()))

    if logical_nonzero:
        science_status = "NONZERO_FINITE_LOGICAL_RETURN_WITNESS"
        implication = (
            "The preregistered finite sine-Lorentzian epsilon sum has nonzero direct logical return on the first tested column. "
            "This establishes P_logical H_L^raw P_logical != 0 on the declared finite habitat, but not a physical scalar pole."
        )
    else:
        science_status = "FIRST_LOGICAL_COLUMN_ZERO_NO_GLOBAL_CONCLUSION"
        implication = (
            "The preregistered first logical column has zero direct return within the frozen tolerance. "
            "No operator-level zero is inferred without the remaining 31 columns."
        )

    return {
        "status": "sharded deterministic aggregate of full epsilon-oriented sine-Lorentzian logical-return witness",
        "passed": passed,
        "science_status": science_status,
        "execution": {
            "shards": shards,
            "term_ids": sorted_ids,
            "summation_schedule_only": True,
            "scientific_protocol_changed": False,
        },
        "source_node": source_v,
        "input_logical_basis_index": input_index,
        "input_K_labels": list(initial[1]),
        "logical_basis_dimension": len(basis),
        "Jmax": FULL.JMAX2 / 2,
        "euclidean_ordering": "H_E^sine=(T-T^dagger)/(2i)",
        "K_definition": "K=[V,H_E^sine]",
        "epsilon_assembly": "sum_r (-1)^r sum_perm sgn(perm) Tr_aux[C_a(K) C_b(K) C_c(V)]",
        "ordered_term_count": 24,
        "full_output_support": len(total),
        "full_output_norm": full_norm,
        "full_scalar_diagnostics": full_scalar,
        "max_spin_reached": max_spin,
        "max_diagnostics": global_diag,
        "historical_primitive_charge_basis_diagnostic": {
            "value": global_diag["CK_complete_charge_basis_leakage"],
            "hard_acceptance": False,
        },
        "reverse_projection": reverse,
        "complete_gauss_outgoing_column": {
            "support": len(gauss),
            "norm": gauss_norm,
            "basis": "Peter-Weyl Gauss basis (spins,K_labels)",
            "state": encode_gauss_state(gauss),
            "reuse": "direct input to Lorentzian Gram/master assembly; no logical projection required",
        },
        "logical_projection": {
            "support": len(logical),
            "norm": logical_norm,
            "fraction_of_full_norm": logical_norm / max(full_norm, 1e-300),
            "nonzero_detection_threshold": FULL.NONZERO_TOL,
            "nonzero": logical_nonzero,
            "initial_state_return_amplitude": [float(initial_return.real), float(initial_return.imag)],
            "amplitudes": logical_rows,
        },
        "hard_integrity_checks": hard,
        "ordered_terms": sorted(term_rows, key=lambda t: int(t["term_id"])),
        "scientific_implication": implication,
        "claim_boundary": (
            "Finite raw Lorentzian constraint-amplitude witness only. It is not the full physical projector, "
            "not Gamma_scalar(omega,k), and not evidence for dark matter or dark energy."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--shards", type=int, default=8)
    ap.add_argument("--source-node", type=int, default=0)
    ap.add_argument("--input-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = run(args.input_dir, args.shards, args.source_node, args.input_index)
    text = json.dumps(out, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "complete_gauss_outgoing_column"}, indent=2))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
