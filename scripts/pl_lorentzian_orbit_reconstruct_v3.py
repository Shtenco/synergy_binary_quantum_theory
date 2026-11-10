#!/usr/bin/env python3
"""Validate and reconstruct the exact V3 pairing-stabilizer Lorentzian orbit.

Input is exactly two corrected-V2 six-term cache-sharing shards:
  forward first-slot=3 -> [0,2,6,8,12,14]
  adjoint first-slot=1 -> [0,1,14,15,20,21]

One direct pair in every H orbit is used as a held-out implementation check.
Only if all six pairs pass is the complete 24+24 term set materialized by the
exact oriented Peter-Weyl state action U_h.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import peter_weyl_zeroaware_volume_migration_experiment as ZVM
import pl_lorentzian_48_collect as COL
import pl_lorentzian_triple_worker as W
from pl_dual_complex import DualComplex, seed_16cell_boundary
from pl_peter_weyl_euclidean import PLPeterWeylEuclidean

VERSION = "tetrahedral-charged-volume-v2"
EXEC_DIRECT = "direct-corrected-v2"
EXEC_RECON = "pairing-stabilizer-reconstructed-v3"
TOL = 1e-8

DIRECT = {
    "forward": [0, 2, 6, 8, 12, 14],
    "adjoint": [0, 1, 14, 15, 20, 21],
}

# Frozen one-representative/one-held-out pair per free orbit.
PAIRS = {
    "forward": [(0, 6), (2, 8), (12, 14)],
    "adjoint": [(0, 1), (15, 21), (14, 20)],
}

EXPECTED_ORBITS = [
    [0, 1, 6, 7, 16, 17, 22, 23],
    [2, 4, 8, 10, 13, 15, 19, 21],
    [3, 5, 9, 11, 12, 14, 18, 20],
]


def parity(p):
    return -1 if sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) % 2 else 1


def pairing_stabilizer():
    pairs = {frozenset((0, 1)), frozenset((2, 3))}
    return tuple(
        p
        for p in itertools.permutations(range(4))
        if {frozenset((p[0], p[1])), frozenset((p[2], p[3]))} == pairs
    )


def map_node(v, h):
    bits = [(v >> (3 - i)) & 1 for i in range(4)]
    nb = [0] * 4
    for i in range(4):
        nb[h[i]] = bits[i]
    out = 0
    for i, b in enumerate(nb):
        out |= b << (3 - i)
    return out


def inverse_perm(h):
    q = [0] * 4
    for i, x in enumerate(h):
        q[x] = i
    return tuple(q)


def mapped_spins(spins, edges, ei, h):
    ns = [0] * len(edges)
    for old, (a, b) in enumerate(edges):
        e = tuple(sorted((map_node(a, h), map_node(b, h))))
        ns[ei[e]] = spins[old]
    return tuple(ns)


def relerr(a, b):
    keys = set(a) | set(b)
    num = math.sqrt(sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys))
    den = math.sqrt(sum(abs(z) ** 2 for z in b.values()))
    return num / max(den, 1e-300)


def state_norm(s):
    return math.sqrt(sum(abs(a) ** 2 for a in s.values()))


def state_max_spin(s):
    return max((max(k[0]) for k in s), default=0) / 2.0


def find_unique(root: Path, name: str) -> Path:
    hits = list(root.rglob(name))
    if len(hits) != 1:
        raise RuntimeError(f"expected one {name}, found {len(hits)}: {hits[:8]}")
    return hits[0]


def load_direct(root: Path, mode: str, idx: int):
    jp = find_unique(root, f"term_{mode}_{idx}.json")
    npz = find_unique(root, f"term_{mode}_{idx}.npz")
    meta = json.loads(jp.read_text(encoding="utf-8"))
    if not meta.get("passed"):
        raise RuntimeError(f"direct worker failed: {jp}")
    if meta.get("operator_version") != VERSION:
        raise RuntimeError(f"wrong operator provenance in {jp}: {meta.get('operator_version')}")
    if meta.get("mode") != mode or int(meta.get("index")) != idx:
        raise RuntimeError(f"metadata identity mismatch in {jp}")
    return COL.load_state(npz), meta


def worker_full_tuple(D, source, idx):
    omit, base, perm, targets, coef = W.ordered_spec(D, source, idx)
    return tuple(perm) + (omit,), (omit, base, perm, targets, coef)


def build_orbits(D, source, H):
    full_to_idx = {}
    idx_to_full = {}
    for idx in range(24):
        full, _ = worker_full_tuple(D, source, idx)
        full_to_idx[full] = idx
        idx_to_full[idx] = full
    unused = set(range(24))
    orbits = []
    while unused:
        i = min(unused)
        f = idx_to_full[i]
        orb = sorted({full_to_idx[tuple(h[x] for x in f)] for h in H})
        orbits.append(orb)
        unused -= set(orb)
    return orbits, idx_to_full, full_to_idx


def h_mapping(idx_to_full, H, src_idx, dst_idx):
    fs = idx_to_full[src_idx]
    fd = idx_to_full[dst_idx]
    hits = [h for h in H if tuple(h[x] for x in fs) == fd]
    if len(hits) != 1:
        raise RuntimeError(("non-unique H transport", src_idx, dst_idx, hits))
    return hits[0]


class StateTransport:
    def __init__(self, D, G):
        self.D = D
        self.G = G
        self.edges = list(G.EDGES)
        self.ei = {e: i for i, e in enumerate(self.edges)}
        self.invs = {h: inverse_perm(h) for h in pairing_stabilizer()}
        self.local_cache = {}
        self.max_local_leak = 0.0
        self.max_phase_mod = 0.0

    def local_phase(self, v, spins, K, h, newspins):
        key = (v, spins, K, h, newspins)
        if key in self.local_cache:
            return self.local_cache[key]
        t = map_node(v, h)
        oldls = self.G.local_spins(spins, v)
        newls = self.G.local_spins(newspins, t)
        expected = [None] * 4
        for r in range(4):
            expected[h[r]] = oldls[r]
        if tuple(expected) != tuple(newls):
            raise RuntimeError(("local spin permutation mismatch", v, h, oldls, newls, expected))
        T = self.G.oriented_intertwiner(v, oldls, K)
        Tp = np.transpose(T, axes=self.invs[h])
        U = self.G.oriented_intertwiner(t, newls, K)
        z = np.vdot(U, Tp)
        leak = float(np.linalg.norm(Tp - z * U))
        mod = float(abs(abs(z) - 1.0))
        self.max_local_leak = max(self.max_local_leak, leak)
        self.max_phase_mod = max(self.max_phase_mod, mod)
        if leak > 1e-9 or mod > 1e-9:
            raise RuntimeError(("H failed to preserve K line", v, h, oldls, K, complex(z), leak, mod))
        self.local_cache[key] = z
        return z

    def map_key_amp(self, key, h):
        spins, Ks = key
        ns = mapped_spins(spins, self.edges, self.ei, h)
        nk = [None] * 16
        phase = 1 + 0j
        for v, K in enumerate(Ks):
            t = map_node(v, h)
            nk[t] = K
            phase *= self.local_phase(v, spins, K, h, ns)
        if any(x is None for x in nk):
            raise RuntimeError("node permutation incomplete")
        return (ns, tuple(nk)), phase

    def map_state(self, state, h):
        out = {}
        for key, a in state.items():
            k, z = self.map_key_amp(key, h)
            out[k] = out.get(k, 0j) + a * z
        return {k: a for k, a in out.items() if abs(a) > 1e-11}


def canonical_meta_from_rep(D, source, mode, idx, rep_idx, rep_meta, state, h):
    omit, base, perm, targets, coef = W.ordered_spec(D, source, idx)
    m = dict(rep_meta)
    m.update(
        {
            "status": "exact pairing-stabilizer reconstructed PL-S3 Lorentzian ordered term",
            "passed": True,
            "operator_version": VERSION,
            "execution_version": EXEC_RECON,
            "symmetry_reconstructed": True,
            "mode": mode,
            "index": idx,
            "source_node": source,
            "omitted_local_slot": omit,
            "base_local_slots": list(base),
            "permuted_local_slots": list(perm),
            "ordered_target_nodes": list(targets),
            "PL_epsilon_coefficient": int(coef),
            "gauss_support": len(state),
            "gauss_norm": state_norm(state),
            "gauss_max_spin": state_max_spin(state),
            "representative_index": rep_idx,
            "transport_permutation": list(h),
            "transport_parity": parity(h),
            "diagnostics_inherited_by_exact_unitary_covariance": True,
            "weighted_here": False,
            "scope_note": (
                "Term reconstructed exactly from a direct corrected-V2 representative after "
                "all six preregistered held-out orbit pairs passed. The PL epsilon coefficient "
                "is independently recomputed from ordered_spec."
            ),
        }
    )
    return m


def run(root: Path, out_dir: Path, source: int = 0):
    ZVM.patch_and_clear()
    D = DualComplex(seed_16cell_boundary())
    G = PLPeterWeylEuclidean(D)
    H = pairing_stabilizer()
    if len(H) != 8:
        raise RuntimeError(f"pairing stabilizer order changed: {len(H)}")

    orbits, idx_to_full, _ = build_orbits(D, source, H)
    if orbits != EXPECTED_ORBITS:
        raise RuntimeError(("worker orbit changed", orbits))

    direct_states = {}
    direct_meta = {}
    for mode, indices in DIRECT.items():
        for idx in indices:
            st, meta = load_direct(root, mode, idx)
            direct_states[(mode, idx)] = st
            direct_meta[(mode, idx)] = meta

    transport = StateTransport(D, G)
    validation = []
    max_pair_error = 0.0
    all_pairs_pass = True

    for mode, pairs in PAIRS.items():
        for orbit_id, (rep, held) in enumerate(pairs):
            # Hard guard: both pair members must belong to the same declared orbit.
            orb = EXPECTED_ORBITS[orbit_id]
            if rep not in orb or held not in orb:
                raise RuntimeError(("frozen pair/orbit mismatch", mode, orbit_id, rep, held, orb))
            h = h_mapping(idx_to_full, H, rep, held)
            mapped = transport.map_state(direct_states[(mode, rep)], h)
            target = direct_states[(mode, held)]
            err = relerr(mapped, target)
            support_equal = set(mapped) == set(target)
            passed = bool(
                direct_meta[(mode, rep)].get("passed")
                and direct_meta[(mode, held)].get("passed")
                and support_equal
                and err < TOL
            )
            max_pair_error = max(max_pair_error, err)
            all_pairs_pass &= passed
            validation.append(
                {
                    "mode": mode,
                    "orbit_id": orbit_id,
                    "representative_index": rep,
                    "heldout_index": held,
                    "transport_permutation": list(h),
                    "transport_parity": parity(h),
                    "support_equal": support_equal,
                    "relative_amplitude_error": err,
                    "representative_norm": state_norm(direct_states[(mode, rep)]),
                    "heldout_norm": state_norm(target),
                    "mapped_norm": state_norm(mapped),
                    "passed": passed,
                }
            )

    if not all_pairs_pass:
        raise RuntimeError(f"one or more held-out orbit pairs failed; max relative error={max_pair_error}")

    out_dir.mkdir(parents=True, exist_ok=True)
    direct_count = 0
    reconstructed_count = 0

    for mode in ("forward", "adjoint"):
        rep_by_orbit = {oid: PAIRS[mode][oid][0] for oid in range(3)}
        for idx in range(24):
            out_npz = out_dir / f"term_{mode}_{idx}.npz"
            out_json = out_dir / f"term_{mode}_{idx}.json"
            if idx in DIRECT[mode]:
                st = direct_states[(mode, idx)]
                meta = dict(direct_meta[(mode, idx)])
                meta["execution_version"] = EXEC_DIRECT
                meta["symmetry_reconstructed"] = False
                meta["direct_validation_term"] = True
                direct_count += 1
            else:
                orbit_id = next(oid for oid, orb in enumerate(EXPECTED_ORBITS) if idx in orb)
                rep = rep_by_orbit[orbit_id]
                h = h_mapping(idx_to_full, H, rep, idx)
                st = transport.map_state(direct_states[(mode, rep)], h)
                meta = canonical_meta_from_rep(
                    D, source, mode, idx, rep, direct_meta[(mode, rep)], st, h
                )
                reconstructed_count += 1
            W.save_state(out_npz, st, len(G.EDGES), D.n_tets)
            # A transported state is a new artifact; never inherit the
            # representative's NPZ digest into reconstructed-term provenance.
            meta["state_sha256"] = hashlib.sha256(out_npz.read_bytes()).hexdigest()
            meta["state_file"] = out_npz.name
            out_json.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = {
        "H_order8": len(H) == 8,
        "three_expected_free_orbits": orbits == EXPECTED_ORBITS,
        "all_six_heldout_pairs_pass": bool(all_pairs_pass and len(validation) == 6),
        "max_pair_relative_error_below_1e-8": max_pair_error < TOL,
        "local_intertwiner_transport_exact": transport.max_local_leak < 1e-9
        and transport.max_phase_mod < 1e-9,
        "direct_term_count_12": direct_count == 12,
        "reconstructed_term_count_36": reconstructed_count == 36,
        "materialized_term_count_48x2": len(list(out_dir.glob("term_*.json"))) == 48
        and len(list(out_dir.glob("term_*.npz"))) == 48,
    }

    return {
        "status": "exact V3 pairing-stabilizer Lorentzian orbit reconstruction",
        "passed": bool(all(checks.values())),
        "science_status": "ORBIT_REDUCTION_VALIDATED_BEFORE_FINAL_S_COLLECTOR",
        "operator_version": VERSION,
        "execution_version": "v3-direct-heldout-plus-exact-H-reconstruction",
        "source_node": source,
        "group_order": len(H),
        "orbits": orbits,
        "forward_direct_indices": DIRECT["forward"],
        "adjoint_direct_indices": DIRECT["adjoint"],
        "direct_term_count": direct_count,
        "reconstructed_term_count": reconstructed_count,
        "validation_pairs": validation,
        "max_pair_relative_amplitude_error": max_pair_error,
        "max_local_intertwiner_line_leakage": transport.max_local_leak,
        "max_local_phase_modulus_defect": transport.max_phase_mod,
        "checks": checks,
        "hard_guard": (
            "No reconstruction is emitted unless all six direct held-out orbit pairs pass. "
            "No GR/HDA target enters representatives, transport, signs, cutoffs or thresholds."
        ),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--certificate", type=Path, required=True)
    p.add_argument("--source", type=int, default=0)
    a = p.parse_args()
    cert = run(a.root, a.out_dir, a.source)
    a.certificate.parent.mkdir(parents=True, exist_ok=True)
    a.certificate.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(cert, indent=2, sort_keys=True))
    return 0 if cert["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
