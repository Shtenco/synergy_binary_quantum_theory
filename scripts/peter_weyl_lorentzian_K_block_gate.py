#!/usr/bin/env python3
"""Exact Gauss-basis Lorentzian prerequisite K=[V,H_E] for the safe K5 engine.

Thiemann's real-connection Lorentzian construction uses an extrinsic-curvature
operator generated from the volume and the Euclidean Hamiltonian.  Before
building the expensive triple C_e(K) Lorentzian term, this gate verifies that
the existing Peter-Weyl recoupling engine can represent the required
commutator with genuine amplitudes rather than support counting.

We use the already symmetrized regulator-safe Euclidean node Hamiltonian H_E,v
at Jmax=5/2 and the exact four-valent absolute-volume operator at the same node:

    K_v := [V_v, H_E,v].

Because both V_v and H_E,v are Hermitian, K_v must be anti-Hermitian.  The gate
checks this with reverse matrix elements on the largest output amplitudes,
checks that no spin exceeds the safe Euclidean wall, and reports the actual
sparse amplitudes.  Overall constants (including the later (1+beta^2)
Lorentzian prefactor) are deliberately not inserted here.

This is a prerequisite amplitude block, not H_L itself and not HDA closure.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import bqg_pruning_audit as PRUNE


@PRUNE.audit_lru_cache(None)
def _local_volume_column_full(key, v: int):
    """Return retained and threshold-candidate coefficients for one V column."""
    spins, Ks = key
    ls = PW.local_spins(spins, v)
    Kin = Ks[v]
    Tin = PW.oriented_intertwiner(v, ls, Kin)
    Tout = PW.apply_volume_tensor_oriented(Tin, ls, v)
    out = {}
    candidates = []
    for Kout in PW.allowed_k2_t(*ls):
        c = np.vdot(PW.oriented_intertwiner(v, ls, Kout), Tout)
        ko = (spins, tuple(Kout if u == v else Ks[u] for u in PW.VERT))
        PRUNE.record_threshold("lorentzian_K_block.local_volume.1e-12",(key,v,Kout),c,1e-12)
        if abs(c) <= 1e-12 and c != 0:
            candidates.append((ko, complex(c)))
        if PRUNE.keep(c,1e-12):
            out[ko] = complex(c)
    return tuple(out.items()), tuple(candidates)


def local_volume_column(key, v: int, incoming: complex = 1.0 + 0j):
    """Return the retained column and audit each removed coefficient after input propagation."""
    retained, candidates = _local_volume_column_full(key, v)
    for output_key, coefficient in candidates:
        PRUNE.record_threshold(
            "lorentzian_K_block.local_volume.propagated_loss",
            (key, output_key), coefficient, 1e-12, incoming=incoming,
        )
    return retained


# Preserve the cache-control/introspection surface used by the audit harness.
local_volume_column.cache_clear = _local_volume_column_full.cache_clear
local_volume_column.cache_info = _local_volume_column_full.cache_info
local_volume_column.cache_parameters = _local_volume_column_full.cache_parameters
local_volume_column.__wrapped__ = (
    lambda key, v, incoming=1.0 + 0j:
    _local_volume_column_full.__wrapped__(key, v)[0]
)


def apply_V_local(state, v: int):
    out = {}
    for key, amp in state.items():
        for ko, c in local_volume_column(key, v, incoming=amp):
            out[ko] = out.get(ko, 0j) + amp * c
    return PW.prune_state(out, 1e-10)


def apply_HE_local(state, v: int, Jmax2: int):
    return PW.prune_state(PW.apply_H_cached_state(state, v, Jmax2), 1e-9)


def apply_K_local(state, v: int, Jmax2: int):
    # K=[V,H]=V H-H V, with composition read right-to-left.
    VH = apply_V_local(apply_HE_local(state, v, Jmax2), v)
    HV = apply_HE_local(apply_V_local(state, v), v, Jmax2)
    out = {}
    PW.add_dict(out, VH, +1)
    PW.add_dict(out, HV, -1)
    return PW.prune_state(out, 1e-9)


def coeff(state, key):
    return complex(state.get(key, 0j))


def run(v=0, reverse_samples=8):
    JMAX2 = 5
    initial = PW.basis_full_jhalf()[0]
    ket = {initial: 1 + 0j}

    Vket = apply_V_local(ket, v)
    Hket = apply_HE_local(ket, v, JMAX2)
    Kket = apply_K_local(ket, v, JMAX2)

    Vnorm = math.sqrt(PW.norm2_state(Vket))
    Hnorm = math.sqrt(PW.norm2_state(Hket))
    Knorm = math.sqrt(PW.norm2_state(Kket))
    max_spin = max((max(k[0]) for k in Kket), default=0) / 2

    # Exact Hermiticity of the local V matrix on every local intertwiner block
    # reached by H|initial>.
    max_V_herm = 0.0
    local_keys = {initial, *Hket.keys()}
    for key in local_keys:
        col = dict(local_volume_column(key, v))
        for ko, a in col.items():
            rev = dict(local_volume_column(ko, v)).get(key, 0j)
            max_V_herm = max(max_V_herm, abs(a - np.conj(rev)))

    # Anti-Hermiticity K^dag=-K.  Reverse columns are expensive, so freeze the
    # largest amplitudes; this is an amplitude identity, not a norm proxy.
    ranked = sorted(Kket.items(), key=lambda kv: abs(kv[1]), reverse=True)
    reverse_rows = []
    max_K_anti = 0.0
    for b, K_ba in ranked[:reverse_samples]:
        K_on_b = apply_K_local({b: 1 + 0j}, v, JMAX2)
        K_ab = coeff(K_on_b, initial)
        defect = abs(K_ba + np.conj(K_ab))
        scale = max(abs(K_ba), abs(K_ab), 1e-30)
        rel = defect / scale
        max_K_anti = max(max_K_anti, rel)
        reverse_rows.append({
            "abs_K_ba": abs(K_ba),
            "abs_K_ab": abs(K_ab),
            "antihermitian_relative_defect": rel,
            "output_max_spin": max(b[0]) / 2,
        })

    # At j=1/2 the initial local absolute volume is scalar on the 2D
    # intertwiner sector, but H_E leaves that sector/spin assignment, making K
    # nontrivial.  Record the input expectation rather than hard-coding it.
    Vexp = coeff(Vket, initial)

    passed = (
        len(Hket) > 0
        and len(Kket) > 0
        and Hnorm > 1e-10
        and Knorm > 1e-10
        and max_spin <= 1.5 + 1e-12
        and max_V_herm < 1e-10
        and max_K_anti < 5e-7
    )
    return {
        "status": "exact safe Peter-Weyl local K=[V,H_E] amplitude gate",
        "passed": bool(passed),
        "Jmax": 2.5,
        "node": v,
        "input": "all ten links j=1/2; all five K=0",
        "V_support": len(Vket),
        "V_norm": Vnorm,
        "V_input_expectation": [Vexp.real, Vexp.imag],
        "H_E_support": len(Hket),
        "H_E_norm": Hnorm,
        "K_support": len(Kket),
        "K_norm": Knorm,
        "max_spin_reached_by_K": max_spin,
        "max_local_volume_hermiticity_error": max_V_herm,
        "reverse_matrix_element_samples": reverse_rows,
        "max_K_antihermitian_relative_defect": max_K_anti,
        "definition": "K_v=[V_v,H_E,v] with genuine four-valent V=sqrt(|J1.(J2xJ3)|)",
        "beta_note": "No beta-dependent coefficient is inserted in K. The later Lorentzian H_L structural block must factor the fixed (1+beta^2) coefficient rather than fit it.",
        "scope_note": (
            "This establishes the exact extrinsic-curvature prerequisite amplitude block in the Gauss Peter-Weyl basis. "
            "It is not yet C_i(K)C_j(K)C_k(V), not the full Lorentzian Hamiltonian, and not an HDA closure result."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--node", type=int, default=0)
    ap.add_argument("--reverse-samples", type=int, default=8)
    ap.add_argument("--output", type=Path)
    a = ap.parse_args()
    out = run(a.node, a.reverse_samples)
    text = json.dumps(out, indent=2)
    print(text)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
