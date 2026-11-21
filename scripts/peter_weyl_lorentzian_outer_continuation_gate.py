#!/usr/bin/env python3
"""Continue a certified Lorentzian middle-prefix through only the outer C(K).

This gate consumes one serialized nonzero middle-prefix packet produced by
peter_weyl_lorentzian_middle_prefix_gate.py and evaluates the two frozen
ordered triples sharing that (b,c) prefix.  It does not recompute C(V) or the
middle C(K), and it derives epsilon signs from the same face_sign*parity
enumeration as the preregistered 24-term evaluator.

Scientific scope: one source node and one logical boundary input only.  A
zero/nonzero pair result is not a global H_L statement and is not a physical
projector result.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import peter_weyl_lorentzian_epsilon_logical_return_gate as FULL
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

SURVIVING_PAIR_INDICES = (0, 1, 3, 4, 6, 7)


def encode_state(state):
    rows = []
    for key, amp in sorted(state.items(), key=lambda kv: repr(kv[0])):
        spins, Kother, J2, M2, K12, K34 = key
        rows.append({
            "spins": list(spins), "Kother": list(Kother),
            "J2": int(J2), "M2": int(M2), "K12": int(K12), "K34": int(K34),
            "amp": [float(amp.real), float(amp.imag)],
        })
    return rows


def decode_state(rows):
    out = {}
    for row in rows:
        key = (
            tuple(int(x) for x in row["spins"]),
            tuple(int(x) for x in row["Kother"]),
            int(row["J2"]), int(row["M2"]), int(row["K12"]), int(row["K34"]),
        )
        amp = complex(float(row["amp"][0]), float(row["amp"][1]))
        if abs(amp) > FULL.TOL:
            out[key] = amp
    return out


def state_hash(state):
    payload = json.dumps(encode_state(state), separators=(",", ":"), sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


def ordered_pairs(source_v):
    n = tuple(FULL.RAW.PW.NEIG[source_v])
    return tuple((b, c) for b in n for c in n if b != c)


def frozen_triples(source_v):
    neighbors = tuple(FULL.RAW.PW.NEIG[source_v])
    rows = []
    term_index = 0
    for r, omitted in enumerate(neighbors):
        base = tuple(x for x in neighbors if x != omitted)
        face_sign = -1 if r % 2 else 1
        for perm in itertools.permutations(base):
            sign = face_sign * FULL.parity(base, perm)
            a, b, c = perm
            rows.append({
                "term_index": term_index, "a": a, "b": b, "c": c,
                "omitted": omitted, "epsilon_sign": int(sign),
            })
            term_index += 1
    if len(rows) != 24:
        raise RuntimeError("frozen epsilon enumeration must contain 24 terms")
    return rows


def continue_prefix(packet, prefix_run_id=None, prefix_head_sha=None):
    errors = []
    if packet.get("passed") is not True:
        errors.append("upstream prefix packet did not pass integrity checks")
    if packet.get("science_status") != "MIDDLE_PREFIX_NONZERO":
        errors.append("upstream packet must be a measured nonzero middle prefix")
    source_v = int(packet.get("source_node", -1))
    input_index = int(packet.get("input_logical_basis_index", -1))
    pair = packet.get("ordered_pair", {})
    b, c = int(pair.get("b", -1)), int(pair.get("c", -1))
    pair_index = int(pair.get("pair_index", -1))
    if pair_index not in SURVIVING_PAIR_INDICES:
        errors.append(f"pair_index {pair_index} is not in the six certified surviving classes")
    pairs = ordered_pairs(source_v) if source_v >= 0 else ()
    if pair_index < 0 or pair_index >= len(pairs) or pairs[pair_index] != (b, c):
        errors.append("ordered pair does not match frozen pair-index enumeration")
    if abs(float(packet.get("Jmax", -1)) - FULL.JMAX2 / 2) > 1e-12:
        errors.append("prefix Jmax does not match frozen Lorentzian evaluator")
    paths = packet.get("paths", [])
    if not isinstance(paths, list) or len(paths) != 8:
        errors.append("prefix packet must contain exactly eight serialized auxiliary paths")
    if errors:
        return {"passed": False, "science_status": "INVALID_PREFIX_PACKET", "errors": errors}

    triples = [r for r in frozen_triples(source_v) if r["b"] == b and r["c"] == c]
    if len(triples) != 2:
        raise RuntimeError("every ordered middle prefix must have exactly two outer continuations")

    ZVM.patch_and_clear()
    pair_total = {}
    global_diag = {
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    max_spin = 0.0
    term_rows = []
    old, caches = FULL.install_sine_ordering()
    try:
        for triple in triples:
            term = {}
            path_rows = []
            for prow in paths:
                i, j, k = (int(x) for x in prow["indices"])
                s2 = decode_state(prow.get("middle_state", []))
                s2_scalar = {key: amp for key, amp in s2.items() if int(key[2]) in (0, 2)}
                if s2_scalar:
                    s3, d3 = FULL.RAW.KCOMP.C_K_component(
                        s2_scalar, source_v, int(triple["a"]), i, j, FULL.JMAX2
                    )
                    FULL.update_diag(global_diag, d3)
                else:
                    s3, d3 = {}, {}
                FULL.add(term, s3)
                max_spin = max(max_spin, FULL.max_spin(s2), FULL.max_spin(s3))
                path_rows.append({
                    "indices": [i, j, k],
                    "middle_support": len(s2),
                    "scalar_relevant_middle_support": len(s2_scalar),
                    "after_outer_CK_support": len(s3),
                    "after_outer_CK_norm": FULL.norm(s3),
                })
            FULL.add(pair_total, term, scale=int(triple["epsilon_sign"]))
            sd = FULL.scalar_diagnostics(term)
            term_rows.append({
                **triple,
                "unsigned_term_support": len(term),
                "unsigned_term_norm": FULL.norm(term),
                "unsigned_term_state_sha256": state_hash(term),
                "scalar_diagnostics": sd,
                "paths": path_rows,
            })
        cache_info = {
            name: {"hits": fun.cache_info().hits, "misses": fun.cache_info().misses, "currsize": fun.cache_info().currsize}
            for name, fun in caches.items()
        }
    finally:
        FULL.restore_ordering(old)

    pair_scalar = FULL.scalar_diagnostics(pair_total)
    hard = {
        "exactly_two_outer_continuations": len(term_rows) == 2,
        "all_epsilon_signs_unit": all(abs(int(r["epsilon_sign"])) == 1 for r in term_rows),
        "outer_complete_basis_leakage_below_1e-9": global_diag["CK_outer_complete_basis_leakage"] < 1e-9,
        "internal_volume_sector_leakage_below_1e-9": global_diag["CK_internal_volume_sector_leakage"] < 1e-9,
        "spin_cutoff_respected": max_spin <= FULL.JMAX2 / 2 + 1e-12,
        "prefix_paths_preserved": len(paths) == 8,
    }
    pair_zero = FULL.norm(pair_total) <= FULL.NONZERO_TOL
    return {
        "status": "Lorentzian outer continuation from certified serialized middle prefix",
        "passed": bool(all(hard.values())),
        "science_status": "PAIR_OUTER_CONTINUATION_ZERO" if pair_zero else "PAIR_OUTER_CONTINUATION_NONZERO",
        "source_node": source_v,
        "input_logical_basis_index": input_index,
        "input_K_labels": packet.get("input_K_labels"),
        "Jmax": FULL.JMAX2 / 2,
        "ordered_pair": {"b": b, "c": c, "pair_index": pair_index},
        "upstream_prefix_provenance": {
            "run_id": prefix_run_id,
            "head_sha": prefix_head_sha,
            "science_status": packet.get("science_status"),
        },
        "continued_full_terms": term_rows,
        "signed_pair_support": len(pair_total),
        "signed_pair_norm": FULL.norm(pair_total),
        "signed_pair_state_sha256": state_hash(pair_total),
        "signed_pair_scalar_diagnostics": pair_scalar,
        "signed_pair_state": encode_state(pair_total),
        "max_spin_reached": max_spin,
        "max_diagnostics": global_diag,
        "hard_integrity_checks": hard,
        "runtime_exact_cache": cache_info,
        "claim_boundary": "Finite continuation for one middle prefix and one preregistered logical input only; no global H_L, physical projector, pole, dark matter or dark energy conclusion follows.",
        "errors": [],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefix-packet", type=Path, required=True)
    ap.add_argument("--prefix-run-id", type=int)
    ap.add_argument("--prefix-head-sha")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    packet = json.loads(args.prefix_packet.read_text(encoding="utf-8"))
    out = continue_prefix(packet, args.prefix_run_id, args.prefix_head_sha)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    slim = {k: v for k, v in out.items() if k not in ("signed_pair_state", "continued_full_terms")}
    print(json.dumps(slim, indent=2))
    return 0 if out.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
