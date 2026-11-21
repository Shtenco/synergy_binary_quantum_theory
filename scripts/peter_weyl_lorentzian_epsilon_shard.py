#!/usr/bin/env python3
"""Execution-only sharding of the preregistered 24-term Lorentzian epsilon sum.

Scientific protocol is unchanged from
PETER_WEYL_LORENTZIAN_LOGICAL_RETURN_PREREGISTRATION.md and
peter_weyl_lorentzian_epsilon_logical_return_gate.py.

Eight default shards each evaluate three of the exact 24 ordered triples with a
shared in-process cache.  Every shard stores its signed sparse partial state and
term diagnostics for deterministic aggregation.  Sharding changes only the
summation schedule, never the frozen operator, cutoff, signs or thresholds.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_zeroaware_volume_migration_experiment as ZVM


def encode_state(state):
    rows = []
    for key, amp in sorted(state.items(), key=lambda kv: repr(kv[0])):
        spins, Kother, J2, M2, K12, K34 = key
        rows.append({
            "spins": list(spins),
            "Kother": list(Kother),
            "J2": int(J2),
            "M2": int(M2),
            "K12": int(K12),
            "K34": int(K34),
            "amp": [float(amp.real), float(amp.imag)],
        })
    return rows


def term_spec(source_v, term_id):
    neighbors = tuple(FULL.RAW.PW.NEIG[source_v])
    if not (0 <= term_id < 24):
        raise ValueError("term_id must be in [0,23]")
    r = term_id // 6
    pidx = term_id % 6
    omitted = neighbors[r]
    base = tuple(x for x in neighbors if x != omitted)
    perms = tuple(itertools.permutations(base))
    perm = perms[pidx]
    face_sign = -1 if r % 2 else 1
    sign = face_sign * FULL.parity(base, perm)
    return {
        "term_id": term_id,
        "face_index": r,
        "permutation_index": pidx,
        "omitted_neighbor": omitted,
        "base_face": base,
        "ordered_edges": perm,
        "sign": sign,
    }


def run(shard=0, shards=8, source_v=0, input_index=0):
    if shards <= 0 or not (0 <= shard < shards):
        raise ValueError("invalid shard/shards")

    ZVM.patch_and_clear()
    basis = FULL.RAW.PW.basis_full_jhalf()
    if len(basis) != 32:
        raise RuntimeError(f"expected 32 logical inputs, found {len(basis)}")
    initial = basis[input_index]
    psi = FULL.RAW.CV.gauss_to_covariant({initial: 1 + 0j}, source_v)

    assigned = [t for t in range(24) if t % shards == shard]
    if not assigned:
        raise RuntimeError("empty shard")

    partial = {}
    terms = []
    global_diag = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    max_spin = 0.0

    old, caches = FULL.install_sine_ordering()
    try:
        for term_id in assigned:
            spec = term_spec(source_v, term_id)
            a, b, c = spec["ordered_edges"]
            state, diag, paths = FULL.ordered_triple_state(
                psi, source_v, a, b, c
            )
            FULL.add(partial, state, scale=spec["sign"])
            FULL.update_diag(global_diag, diag)
            sd = FULL.scalar_diagnostics(state)
            mspin = FULL.max_spin(state)
            max_spin = max(max_spin, mspin)

            if sd["norm"] > FULL.NONZERO_TOL:
                scalar_ok = sd["nonscalar_weight_fraction"] < 1e-8
            else:
                scalar_ok = sd["nonscalar_norm"] < FULL.NONZERO_TOL

            hard = {
                "term_scalar": bool(scalar_ok),
                "CV_leakage": bool(diag["CV_complete_basis_leakage"] < 1e-9),
                "CK_outer_leakage": bool(diag["CK_outer_complete_basis_leakage"] < 1e-9),
                "CK_internal_volume_leakage": bool(diag["CK_internal_volume_sector_leakage"] < 1e-9),
                "spin_cutoff": bool(mspin <= FULL.JMAX2 / 2 + 1e-12),
            }
            terms.append({
                **{k: (list(v) if isinstance(v, tuple) else v) for k, v in spec.items()},
                "support": len(state),
                "scalar_diagnostics": sd,
                "max_spin_reached": mspin,
                "hard_checks": hard,
                "passed": bool(all(hard.values())),
                "path_supports": paths,
            })

        cache_info = {
            name: {
                "hits": fun.cache_info().hits,
                "misses": fun.cache_info().misses,
                "currsize": fun.cache_info().currsize,
            }
            for name, fun in caches.items()
        }
    finally:
        FULL.restore_ordering(old)

    partial_diag = FULL.scalar_diagnostics(partial)
    passed = bool(all(t["passed"] for t in terms))
    return {
        "status": "execution shard of preregistered full epsilon sine-Lorentzian sum",
        "passed": passed,
        "science_status": "EXECUTION_SHARD_ONLY",
        "source_node": source_v,
        "input_logical_basis_index": input_index,
        "input_K_labels": list(initial[1]),
        "Jmax": FULL.JMAX2 / 2,
        "shard": shard,
        "shards": shards,
        "assigned_term_ids": assigned,
        "term_count": len(assigned),
        "terms": terms,
        "signed_partial_support": len(partial),
        "signed_partial_scalar_diagnostics": partial_diag,
        "max_spin_reached": max_spin,
        "max_diagnostics": global_diag,
        "signed_partial_state": encode_state(partial),
        "runtime_exact_cache": cache_info,
        "claim_boundary": (
            "Execution decomposition only. No shard has independent physical meaning; "
            "only the deterministic aggregate of all 24 preregistered terms is interpreted."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--shards", type=int, default=8)
    ap.add_argument("--source-node", type=int, default=0)
    ap.add_argument("--input-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = run(args.shard, args.shards, args.source_node, args.input_index)
    text = json.dumps(out, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "signed_partial_state"}, indent=2))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
