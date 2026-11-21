#!/usr/bin/env python3
"""Direct logical projection of the full epsilon-oriented raw Lorentzian node operator.

This is the first finite killer test of

    P L_raw,epsilon P

on the real sine-ordered Peter-Weyl Lorentzian amplitude stack.

The operator core is the existing genuine state-to-state triple

    Tr_aux[C_a(K_sine) C_b(K_sine) C_c(V)],
    K_sine=[V,H_E^sine],

at the preregistered single-H_L wall Jmax=7/2.  For one source node the gate
sums all 24 oriented terms (four omitted-face choices times six permutations)
with the same epsilon signs checked independently by
`peter_weyl_lorentzian_triple_algebra_gate.py`.

The environment is frozen to K=0 on the other four K5 logical nodes for this
first falsifier.  The source logical input is evaluated for both K=0 and K=2.
After the traced triple, the source covariant J=0 key is mapped back to the
ordinary Gauss logical key exactly via

    (J2,M2,K12,K34)=(0,0,K,K).

The output reports:
- the 2x2 source-logical matrix with the environment held at K=0;
- all-j=1/2 logical support that changes the environment;
- nonlogical/spin-changed remainder;
- mirror Y-odd/even Pauli content of the 2x2 local block.

Important scope: L_raw,epsilon is the full oriented structural K-K-V sum but is
not yet declared to be the final Hermitian H_L normalization.  A nonzero result
means only that support and amplitudes permit a direct logical Lorentzian term.
A physical mass/anisotropy requires the correct Hermitian/prefactor completion,
route/matter coupling and RG analysis.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_K_composition_gate as KC
import peter_weyl_covariant_K_composition_cached_gate as CACHE
import peter_weyl_covariant_K_sine_composition_gate as SINEK
import peter_weyl_lorentzian_ordered_triple_gate as RAW
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

JMAX2 = 7
TOL = 1e-11

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def freeze(state):
    return tuple(sorted(state.items(), key=lambda kv: repr(kv[0])))


def add(dst, src, scale=1.0, tol=TOL):
    for key, amp in src.items():
        z = dst.get(key, 0j) + scale * amp
        if abs(z) > tol:
            dst[key] = z
        elif key in dst:
            del dst[key]


def norm2(state):
    return float(sum(abs(a) ** 2 for a in state.values()))


def parity(base, perm):
    idx = [base.index(x) for x in perm]
    inv = sum(idx[i] > idx[j] for i in range(3) for j in range(i + 1, 3))
    return -1 if inv % 2 else +1


def install_sine_cached_stack():
    """Patch the already validated Lorentzian stack once and return restore()."""
    ZVM.patch_and_clear()

    old_he = KC.CK.apply_HE_complete_key
    old_gauss = KC.KG.apply_HE_local
    old_k = KC.apply_K_complete_custom
    old_inv = KC.COMP.inverse_complete
    old_direct = KC.direct_K_covariant
    old_close = KC.COMP.close_complete
    old_ck = KC.C_K_component
    old_raw_ck = RAW.KCOMP.C_K_component
    old_raw_cv = RAW.COMP.C_volume_component

    @functools.lru_cache(maxsize=None)
    def he_reduced(canonical_key, source_v, Jmax2, charged_nodes):
        return SINEK.complete_HE_sine(
            canonical_key, source_v, Jmax2, charged_nodes=tuple(charged_nodes)
        )

    def he_sine_reduced(key, source_v, Jmax2, charged_nodes=(0, 1)):
        charged_nodes = tuple(charged_nodes)
        canonical, original = CACHE.canonicalize_scalar_charge_M(key, charged_nodes)
        state, vleak, bleak = he_reduced(canonical, source_v, Jmax2, charged_nodes)
        return CACHE.restore_scalar_charge_M(state, charged_nodes, original), vleak, bleak

    @functools.lru_cache(maxsize=None)
    def k_cached(frozen, source_v, Jmax2, charged_nodes):
        out, vleak, bleak = old_k(dict(frozen), source_v, Jmax2, tuple(charged_nodes))
        return tuple(out.items()), float(vleak), float(bleak)

    def k_wrap(state, source_v, Jmax2, charged_nodes):
        items, vleak, bleak = k_cached(freeze(state), source_v, Jmax2, tuple(charged_nodes))
        return dict(items), vleak, bleak

    @functools.lru_cache(maxsize=None)
    def inv_cached(frozen, source_v, target_v, k, j, Jmax2):
        out, leak = old_inv(dict(frozen), source_v, target_v, k, j, Jmax2)
        return tuple(out.items()), float(leak)

    def inv_wrap(state, source_v, target_v, k, j, Jmax2):
        items, leak = inv_cached(freeze(state), source_v, target_v, k, j, Jmax2)
        return dict(items), leak

    @functools.lru_cache(maxsize=None)
    def direct_cached(frozen, source_v, Jmax2):
        out, vleak, bleak = old_direct(dict(frozen), source_v, Jmax2)
        return tuple(out.items()), float(vleak), float(bleak)

    def direct_wrap(state, source_v, Jmax2):
        items, vleak, bleak = direct_cached(freeze(state), source_v, Jmax2)
        return dict(items), vleak, bleak

    @functools.lru_cache(maxsize=None)
    def close_cached(frozen, source_v, target_v, i, k, Jmax2):
        return tuple(old_close(dict(frozen), source_v, target_v, i, k, Jmax2).items())

    def close_wrap(state, source_v, target_v, i, k, Jmax2):
        return dict(close_cached(freeze(state), source_v, target_v, i, k, Jmax2))

    @functools.lru_cache(maxsize=None)
    def ck_cached(frozen, source_v, target_v, i, j, Jmax2):
        out, diag = old_ck(dict(frozen), source_v, target_v, i, j, Jmax2)
        return tuple(out.items()), tuple(sorted(diag.items()))

    def ck_wrap(state, source_v, target_v, i, j, Jmax2):
        items, diag = ck_cached(freeze(state), source_v, target_v, i, j, Jmax2)
        return dict(items), dict(diag)

    @functools.lru_cache(maxsize=None)
    def cv_cached(frozen, source_v, target_v, i, j, Jmax2):
        out, leak = old_raw_cv(dict(frozen), source_v, target_v, i, j, Jmax2)
        return tuple(out.items()), float(leak)

    def cv_wrap(state, source_v, target_v, i, j, Jmax2):
        items, leak = cv_cached(freeze(state), source_v, target_v, i, j, Jmax2)
        return dict(items), leak

    KC.CK.apply_HE_complete_key = he_sine_reduced
    KC.KG.apply_HE_local = SINEK.gauss_HE_sine_with_historical_K_cutoff
    KC.apply_K_complete_custom = k_wrap
    KC.COMP.inverse_complete = inv_wrap
    KC.direct_K_covariant = direct_wrap
    KC.COMP.close_complete = close_wrap
    KC.C_K_component = ck_wrap
    RAW.KCOMP.C_K_component = ck_wrap
    RAW.COMP.C_volume_component = cv_wrap
    if hasattr(KC.CK.HE_complete_cached, "cache_clear"):
        KC.CK.HE_complete_cached.cache_clear()

    def restore():
        KC.CK.apply_HE_complete_key = old_he
        KC.KG.apply_HE_local = old_gauss
        KC.apply_K_complete_custom = old_k
        KC.COMP.inverse_complete = old_inv
        KC.direct_K_covariant = old_direct
        KC.COMP.close_complete = old_close
        KC.C_K_component = old_ck
        RAW.KCOMP.C_K_component = old_raw_ck
        RAW.COMP.C_volume_component = old_raw_cv
        if hasattr(KC.CK.HE_complete_cached, "cache_clear"):
            KC.CK.HE_complete_cached.cache_clear()

    cache_info = {
        "HE_reduced": he_reduced,
        "K_complete": k_cached,
        "inverse": inv_cached,
        "direct_K": direct_cached,
        "close": close_cached,
        "C_K": ck_cached,
        "C_V": cv_cached,
    }
    return restore, cache_info


def update_diag(dst, name, value):
    dst[name] = max(dst.get(name, 0.0), float(value))


def ordered_triple_state(initial, source_v, a, b, c):
    """Return the actual covariant state for Tr[C_a(K)C_b(K)C_c(V)]|initial>."""
    psi = CV.gauss_to_covariant({initial: 1 + 0j}, source_v)
    total = {}
    diag = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    for i, j, k in itertools.product(range(2), repeat=3):
        s1, leakV = RAW.COMP.C_volume_component(psi, source_v, c, k, i, JMAX2)
        update_diag(diag, "CV_complete_basis_leakage", leakV)
        if not s1:
            continue
        s2, d2 = RAW.KCOMP.C_K_component(s1, source_v, b, j, k, JMAX2)
        for name, val in (
            ("CK_outer_complete_basis_leakage", d2["outer_complete_basis_leakage"]),
            ("CK_internal_volume_sector_leakage", d2["internal_volume_sector_leakage"]),
            ("CK_complete_charge_basis_leakage", d2["complete_charge_basis_leakage"]),
        ):
            update_diag(diag, name, val)
        # exact scalar-channel pruning used by the existing raw gate
        s2 = {key: amp for key, amp in s2.items() if key[2] in (0, 2)}
        if not s2:
            continue
        s3, d3 = RAW.KCOMP.C_K_component(s2, source_v, a, i, j, JMAX2)
        for name, val in (
            ("CK_outer_complete_basis_leakage", d3["outer_complete_basis_leakage"]),
            ("CK_internal_volume_sector_leakage", d3["internal_volume_sector_leakage"]),
            ("CK_complete_charge_basis_leakage", d3["complete_charge_basis_leakage"]),
        ):
            update_diag(diag, name, val)
        add(total, s3)
    return total, diag


def epsilon_sum_state(initial, source_v):
    neighbors = PW.NEIG[source_v]
    total = {}
    rows = []
    diagmax = {}
    for r, omit in enumerate(neighbors):
        base = tuple(x for x in neighbors if x != omit)
        face = (-1) ** r
        for perm in itertools.permutations(base):
            a, b, c = perm
            coef = face * parity(base, perm)
            st, diag = ordered_triple_state(initial, source_v, a, b, c)
            add(total, st, coef)
            for name, val in diag.items():
                update_diag(diagmax, name, val)
            rows.append({
                "omitted_neighbor": omit,
                "ordered_edges": [a, b, c],
                "coefficient": coef,
                "support": len(st),
                "norm": math.sqrt(norm2(st)),
            })
    return total, rows, diagmax


def covariant_to_gauss_logical(key, source_v):
    """Map a final scalar covariant key to a Gauss logical key, or return None."""
    spins, Kother, J2, M2, K12, K34 = key
    if tuple(spins) != (1,) * len(PW.EDGES):
        return None
    if J2 != 0 or M2 != 0 or K12 != K34:
        return None
    Ks = list(Kother)
    Ks[source_v] = K12
    Ks = tuple(Ks)
    if any(k not in (0, 2) for k in Ks):
        return None
    return (tuple(spins), Ks)


def project_all_logical(state, source_v):
    logical = {}
    for key, amp in state.items():
        g = covariant_to_gauss_logical(key, source_v)
        if g is not None:
            logical[g] = logical.get(g, 0j) + amp
    return {k: v for k, v in logical.items() if abs(v) > TOL}


def complex_pair(z):
    return [float(z.real), float(z.imag)]


def pauli_decompose_2(M):
    return {a: complex_pair(np.trace(A @ M) / 2.0) for a, A in PAULI.items()}


def run(source_v=0):
    restore, caches = install_sine_cached_stack()
    try:
        spins = (1,) * len(PW.EDGES)
        env_K = [0] * len(PW.VERT)
        source_basis = (0, 2)
        columns = []
        epsilon_rows = []
        diagmax = {}
        all_logical_support = []

        for Kin in source_basis:
            Ks = list(env_K)
            Ks[source_v] = Kin
            initial = (spins, tuple(Ks))
            total, rows, diag = epsilon_sum_state(initial, source_v)
            logical = project_all_logical(total, source_v)
            columns.append((initial, total, logical))
            epsilon_rows.append({"input_source_K2": Kin, "terms": rows})
            for name, val in diag.items():
                update_diag(diagmax, name, val)
            all_logical_support.append({
                "input_source_K2": Kin,
                "full_output_support": len(total),
                "full_output_norm": math.sqrt(norm2(total)),
                "all_jhalf_logical_support": len(logical),
                "all_jhalf_logical_norm": math.sqrt(norm2(logical)),
                "logical_outputs": [
                    {"Ks2": list(key[1]), "amp": complex_pair(amp), "abs_amp": abs(amp)}
                    for key, amp in sorted(logical.items(), key=lambda kv: abs(kv[1]), reverse=True)[:24]
                ],
            })

        M = np.zeros((2, 2), dtype=complex)
        env_excited_norm2 = [0.0, 0.0]
        for col, (initial, total, logical) in enumerate(columns):
            for row, Kout in enumerate(source_basis):
                target_Ks = [0] * len(PW.VERT)
                target_Ks[source_v] = Kout
                target = (spins, tuple(target_Ks))
                M[row, col] = logical.get(target, 0j)
            local_targets = {
                (spins, tuple([Kout if u == source_v else 0 for u in PW.VERT]))
                for Kout in source_basis
            }
            env_excited_norm2[col] = float(sum(abs(a) ** 2 for k, a in logical.items() if k not in local_targets))

        local_norm = float(np.linalg.norm(M))
        herm = float(np.linalg.norm(M - M.conj().T))
        antiherm = float(np.linalg.norm(M + M.conj().T))
        coeff = pauli_decompose_2(M)

        # Mirror in the logical geometry qubit is complex conjugation; a
        # Hermitian mirror-even local operator has I/X/Z and no Y.  For the raw
        # non-Hermitian structural sum report all components without assuming
        # the final Hermitian completion.
        cache_info = {
            name: {
                "hits": fn.cache_info().hits,
                "misses": fn.cache_info().misses,
                "currsize": fn.cache_info().currsize,
            }
            for name, fn in caches.items()
        }

        passed = (
            len(columns) == 2
            and all(x["full_output_support"] > 0 for x in all_logical_support)
            and max(diagmax.get("CV_complete_basis_leakage", 0.0),
                    diagmax.get("CK_outer_complete_basis_leakage", 0.0),
                    diagmax.get("CK_internal_volume_sector_leakage", 0.0),
                    diagmax.get("CK_complete_charge_basis_leakage", 0.0)) < 1e-8
        )

        return {
            "status": "direct logical projection of full epsilon-oriented raw Lorentzian K-K-V node sum",
            "passed": bool(passed),
            "source_node": source_v,
            "Jmax": JMAX2 / 2,
            "input_environment": "all other logical K=0; source K in {0,2}",
            "orientation_sum": "4 omitted faces x 6 signed permutations = 24 ordered triples per input",
            "local_2x2_raw_matrix": [[complex_pair(M[r, c]) for c in range(2)] for r in range(2)],
            "local_2x2_frobenius_norm": local_norm,
            "local_2x2_hermiticity_defect_norm": herm,
            "local_2x2_antihermiticity_defect_norm": antiherm,
            "local_2x2_pauli_coefficients": coeff,
            "environment_excited_logical_norm2_by_input": env_excited_norm2,
            "columns": all_logical_support,
            "max_diagnostics": diagmax,
            "cache_info": cache_info,
            "raw_projection_nonzero": bool(local_norm > 1e-10),
            "interpretation": (
                "If raw_projection_nonzero is false, the full oriented sine-ordered K-K-V core has no direct local logical matrix element in this frozen environment. "
                "If true, support and amplitudes permit a direct Lorentzian logical term; the final Hermitian Lorentzian completion and unbiased environment trace are then required before interpreting it as a physical anisotropy or mass."
            ),
            "scope": (
                "Finite candidate-model killer column. No (1+beta^2) prefactor, no final Hermitian H_L normalization, no route/matter sector and no physical force claim."
            ),
        }
    finally:
        restore()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--node", type=int, default=0)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.node)
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
