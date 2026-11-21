#!/usr/bin/env python3
"""Full epsilon-oriented sine-Lorentzian logical-return witness.

This is the preregistered next amplitude test from
PETER_WEYL_LORENTZIAN_LOGICAL_RETURN_PREREGISTRATION.md.

It assembles all 24 ordered terms of

    H_L^raw ~ epsilon * Tr_aux[C(K_sine) C(K_sine) C(V)]

on one all-j=1/2 logical K5 input at Jmax=7/2, preserving the exact
noncommuting geometry-operator order.  The gate's PASS means only that the
finite computation satisfies the frozen scalar/covariance/cutoff diagnostics.
A zero or nonzero logical return is a reported scientific outcome, not a
pass/fail target.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
from pathlib import Path

import peter_weyl_covariant_K_composition_gate as KC
import peter_weyl_covariant_K_composition_cached_gate as CACHE
import peter_weyl_covariant_K_sine_composition_gate as SINEK
import peter_weyl_lorentzian_ordered_triple_gate as RAW
import peter_weyl_zeroaware_volume_migration_experiment as ZVM

JMAX2 = 7
TOL = 1e-11
NONZERO_TOL = 1e-10


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


def norm(state):
    return math.sqrt(norm2(state))


def max_spin(state):
    return max((max(key[0]) for key in state), default=0) / 2


def source_J_weights(state):
    out = {}
    for key, amp in state.items():
        J2 = int(key[2])
        out[J2] = out.get(J2, 0.0) + abs(amp) ** 2
    return out


def scalar_diagnostics(state):
    weights = source_J_weights(state)
    total2 = float(sum(weights.values()))
    scalar2 = float(weights.get(0, 0.0))
    nonscalar2 = max(0.0, total2 - scalar2)
    total = math.sqrt(total2)
    nonscalar = math.sqrt(nonscalar2)
    frac = nonscalar2 / max(total2, 1e-300) if total2 else 0.0
    return {
        "norm": total,
        "nonscalar_norm": nonscalar,
        "nonscalar_weight_fraction": frac,
        "weight_by_source_J": {str(k / 2): float(v) for k, v in sorted(weights.items())},
    }


def scalar_ok(diag):
    if diag["norm"] > NONZERO_TOL:
        return diag["nonscalar_weight_fraction"] < 1e-8
    return diag["nonscalar_norm"] < NONZERO_TOL


def update_diag(dst, src):
    for name, val in src.items():
        if isinstance(val, (int, float)):
            dst[name] = max(dst.get(name, 0.0), float(val))


def parity(base, perm):
    idx = [base.index(x) for x in perm]
    inv = sum(idx[i] > idx[j] for i in range(3) for j in range(i + 1, 3))
    return -1 if inv % 2 else 1


def project_covariant_J0_to_gauss(state, source_v):
    """Invert CV.gauss_to_covariant on final source-J=0 states."""
    out = {}
    mapped_keys = []
    invalid = []
    collisions = 0
    for key, amp in state.items():
        spins, Kother, J2, M2, K12, K34 = key
        if J2 != 0:
            continue
        if M2 != 0 or K12 != K34:
            invalid.append(repr(key))
            continue
        Ks = list(Kother)
        if Ks[source_v] not in (-1, K12):
            invalid.append(repr(key))
            continue
        Ks[source_v] = K12
        if any(k < 0 for k in Ks):
            invalid.append(repr(key))
            continue
        gkey = (spins, tuple(Ks))
        if gkey in out:
            collisions += 1
        out[gkey] = out.get(gkey, 0j) + amp
        mapped_keys.append(gkey)
    return out, {
        "mapped_covariant_J0_basis_states": len(mapped_keys),
        "distinct_gauss_basis_states": len(set(mapped_keys)),
        "mapping_collisions": collisions,
        "invalid_J0_covariant_keys": invalid,
    }


def logical_projection(gauss_state):
    return {
        key: amp
        for key, amp in gauss_state.items()
        if all(int(s) == 1 for s in key[0])
    }


def ordered_triple_state(psi, source_v, a, b, c):
    total = {}
    diagmax = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    path_rows = []

    for i, j, k in itertools.product(range(2), repeat=3):
        s1, leakV = RAW.COMP.C_volume_component(
            psi, source_v, c, k, i, JMAX2
        )
        diagmax["CV_complete_basis_leakage"] = max(
            diagmax["CV_complete_basis_leakage"], float(leakV)
        )
        if not s1:
            path_rows.append({
                "indices": [i, j, k],
                "after_CV_support": 0,
                "after_middle_CK_support": 0,
                "after_middle_CK_scalar_relevant_support": 0,
                "after_final_CK_support": 0,
            })
            continue

        s2, d2 = RAW.KCOMP.C_K_component(
            s1, source_v, b, j, k, JMAX2
        )
        update_diag(diagmax, d2)

        # Exact scalar-channel pruning: the final rank-(0+1) C(K) cannot
        # couple source J=2 back to final J=0.
        s2_scalar = {key: amp for key, amp in s2.items() if key[2] in (0, 2)}
        if not s2_scalar:
            path_rows.append({
                "indices": [i, j, k],
                "after_CV_support": len(s1),
                "after_middle_CK_support": len(s2),
                "after_middle_CK_scalar_relevant_support": 0,
                "after_final_CK_support": 0,
            })
            continue

        s3, d3 = RAW.KCOMP.C_K_component(
            s2_scalar, source_v, a, i, j, JMAX2
        )
        update_diag(diagmax, d3)
        add(total, s3)
        path_rows.append({
            "indices": [i, j, k],
            "after_CV_support": len(s1),
            "after_middle_CK_support": len(s2),
            "after_middle_CK_scalar_relevant_support": len(s2_scalar),
            "after_final_CK_support": len(s3),
        })

    return total, diagmax, path_rows


def install_sine_ordering():
    """Patch the existing exact composition engine exactly as the frozen sine gate."""
    old = {
        "he_complete": KC.CK.apply_HE_complete_key,
        "he_gauss": KC.KG.apply_HE_local,
        "k_complete": KC.apply_K_complete_custom,
        "inverse": KC.COMP.inverse_complete,
        "direct_K": KC.direct_K_covariant,
        "close": KC.COMP.close_complete,
        "C_K": KC.C_K_component,
        "RAW_C_K": RAW.KCOMP.C_K_component,
        "RAW_C_V": RAW.COMP.C_volume_component,
    }

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
        out, vleak, bleak = old["k_complete"](
            dict(frozen), source_v, Jmax2, tuple(charged_nodes)
        )
        return tuple(out.items()), float(vleak), float(bleak)

    def k_wrap(state, source_v, Jmax2, charged_nodes):
        items, vleak, bleak = k_cached(
            freeze(state), source_v, Jmax2, tuple(charged_nodes)
        )
        return dict(items), vleak, bleak

    @functools.lru_cache(maxsize=None)
    def inv_cached(frozen, source_v, target_v, k, j, Jmax2):
        out, leak = old["inverse"](
            dict(frozen), source_v, target_v, k, j, Jmax2
        )
        return tuple(out.items()), float(leak)

    def inv_wrap(state, source_v, target_v, k, j, Jmax2):
        items, leak = inv_cached(
            freeze(state), source_v, target_v, k, j, Jmax2
        )
        return dict(items), leak

    @functools.lru_cache(maxsize=None)
    def direct_cached(frozen, source_v, Jmax2):
        out, vleak, bleak = old["direct_K"](dict(frozen), source_v, Jmax2)
        return tuple(out.items()), float(vleak), float(bleak)

    def direct_wrap(state, source_v, Jmax2):
        items, vleak, bleak = direct_cached(freeze(state), source_v, Jmax2)
        return dict(items), vleak, bleak

    @functools.lru_cache(maxsize=None)
    def close_cached(frozen, source_v, target_v, i, k, Jmax2):
        return tuple(old["close"](
            dict(frozen), source_v, target_v, i, k, Jmax2
        ).items())

    def close_wrap(state, source_v, target_v, i, k, Jmax2):
        return dict(close_cached(
            freeze(state), source_v, target_v, i, k, Jmax2
        ))

    @functools.lru_cache(maxsize=None)
    def ck_cached(frozen, source_v, target_v, i, j, Jmax2):
        out, diag = old["C_K"](
            dict(frozen), source_v, target_v, i, j, Jmax2
        )
        return tuple(out.items()), tuple(sorted(diag.items()))

    def ck_wrap(state, source_v, target_v, i, j, Jmax2):
        items, diag = ck_cached(
            freeze(state), source_v, target_v, i, j, Jmax2
        )
        return dict(items), dict(diag)

    @functools.lru_cache(maxsize=None)
    def cv_cached(frozen, source_v, target_v, i, j, Jmax2):
        out, leak = old["RAW_C_V"](
            dict(frozen), source_v, target_v, i, j, Jmax2
        )
        return tuple(out.items()), float(leak)

    def cv_wrap(state, source_v, target_v, i, j, Jmax2):
        items, leak = cv_cached(
            freeze(state), source_v, target_v, i, j, Jmax2
        )
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

    caches = {
        "HE_reduced": he_reduced,
        "K_complete": k_cached,
        "inverse": inv_cached,
        "direct_K": direct_cached,
        "close": close_cached,
        "C_K": ck_cached,
        "C_V": cv_cached,
    }
    return old, caches


def restore_ordering(old):
    KC.CK.apply_HE_complete_key = old["he_complete"]
    KC.KG.apply_HE_local = old["he_gauss"]
    KC.apply_K_complete_custom = old["k_complete"]
    KC.COMP.inverse_complete = old["inverse"]
    KC.direct_K_covariant = old["direct_K"]
    KC.COMP.close_complete = old["close"]
    KC.C_K_component = old["C_K"]
    RAW.KCOMP.C_K_component = old["RAW_C_K"]
    RAW.COMP.C_volume_component = old["RAW_C_V"]
    if hasattr(KC.CK.HE_complete_cached, "cache_clear"):
        KC.CK.HE_complete_cached.cache_clear()


def run(source_v=0, input_index=0):
    ZVM.patch_and_clear()
    basis = RAW.PW.basis_full_jhalf()
    if len(basis) != 32:
        raise RuntimeError(f"expected 32 all-j=1/2 logical states, found {len(basis)}")
    if not (0 <= input_index < len(basis)):
        raise ValueError("input_index outside logical basis")

    initial = basis[input_index]
    psi = RAW.CV.gauss_to_covariant({initial: 1 + 0j}, source_v)
    neighbors = tuple(RAW.PW.NEIG[source_v])
    if len(neighbors) != 4:
        raise RuntimeError("four-valent source required")

    total = {}
    ordered_rows = []
    global_diag = {
        "CV_complete_basis_leakage": 0.0,
        "CK_outer_complete_basis_leakage": 0.0,
        "CK_internal_volume_sector_leakage": 0.0,
        "CK_complete_charge_basis_leakage": 0.0,
    }
    max_term_nonscalar_fraction = 0.0
    max_term_nonscalar_norm_if_near_zero = 0.0
    max_reached_spin = 0.0

    old, caches = install_sine_ordering()
    try:
        for r, omitted in enumerate(neighbors):
            base = tuple(x for x in neighbors if x != omitted)
            face_sign = -1 if r % 2 else 1
            for perm in itertools.permutations(base):
                sign = face_sign * parity(base, perm)
                a, b, c = perm
                state, diag, paths = ordered_triple_state(
                    psi, source_v, a, b, c
                )
                add(total, state, scale=sign)
                update_diag(global_diag, diag)
                sd = scalar_diagnostics(state)
                if sd["norm"] > NONZERO_TOL:
                    max_term_nonscalar_fraction = max(
                        max_term_nonscalar_fraction,
                        sd["nonscalar_weight_fraction"],
                    )
                else:
                    max_term_nonscalar_norm_if_near_zero = max(
                        max_term_nonscalar_norm_if_near_zero,
                        sd["nonscalar_norm"],
                    )
                max_reached_spin = max(max_reached_spin, max_spin(state))
                ordered_rows.append({
                    "omitted_neighbor": omitted,
                    "face_index": r,
                    "base_face": list(base),
                    "ordered_edges": [a, b, c],
                    "sign": sign,
                    "support": len(state),
                    "scalar_diagnostics": sd,
                    "max_spin_reached": max_spin(state),
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
        restore_ordering(old)

    full_scalar = scalar_diagnostics(total)
    gauss, reverse = project_covariant_J0_to_gauss(total, source_v)
    logical = logical_projection(gauss)
    logical_norm = norm(logical)
    full_norm = norm(total)
    logical_nonzero = logical_norm > NONZERO_TOL
    initial_return = complex(logical.get(initial, 0j))

    ranked_logical = sorted(
        logical.items(), key=lambda kv: abs(kv[1]), reverse=True
    )
    logical_rows = [
        {
            "K_labels": list(key[1]),
            "abs_amp": float(abs(amp)),
            "amp": [float(amp.real), float(amp.imag)],
        }
        for key, amp in ranked_logical
    ]

    hard_checks = {
        "exactly_24_ordered_terms": len(ordered_rows) == 24,
        "all_ordered_terms_scalar_within_frozen_threshold": (
            max_term_nonscalar_fraction < 1e-8
            and max_term_nonscalar_norm_if_near_zero < NONZERO_TOL
        ),
        "full_epsilon_sum_scalar_within_frozen_threshold": scalar_ok(full_scalar),
        "CV_complete_basis_leakage_below_1e-9": (
            global_diag["CV_complete_basis_leakage"] < 1e-9
        ),
        "CK_outer_complete_basis_leakage_below_1e-9": (
            global_diag["CK_outer_complete_basis_leakage"] < 1e-9
        ),
        "CK_internal_volume_sector_leakage_below_1e-9": (
            global_diag["CK_internal_volume_sector_leakage"] < 1e-9
        ),
        "spin_cutoff_respected": max_reached_spin <= JMAX2 / 2 + 1e-12,
        "J0_covariant_to_gauss_mapping_valid": (
            not reverse["invalid_J0_covariant_keys"]
            and reverse["mapping_collisions"] == 0
            and reverse["mapped_covariant_J0_basis_states"]
            == reverse["distinct_gauss_basis_states"]
        ),
    }
    passed = bool(all(hard_checks.values()))

    if logical_nonzero:
        result_class = "NONZERO_FINITE_LOGICAL_RETURN_WITNESS"
        implication = (
            "For the declared finite habitat, sine ordering and Jmax=7/2, this one "
            "column is sufficient to establish P_logical H_L^raw P_logical != 0."
        )
    else:
        result_class = "FIRST_LOGICAL_COLUMN_ZERO_NO_GLOBAL_CONCLUSION"
        implication = (
            "The preregistered first logical column is zero within 1e-10. This does "
            "not establish P_logical H_L P_logical=0; the remaining 31 columns are required."
        )

    return {
        "status": "full epsilon-oriented sine-Lorentzian logical-return witness",
        "passed": passed,
        "science_status": result_class,
        "source_node": source_v,
        "neighbors_in_frozen_order": list(neighbors),
        "input_logical_basis_index": input_index,
        "input_K_labels": list(initial[1]),
        "logical_basis_dimension": len(basis),
        "Jmax": JMAX2 / 2,
        "euclidean_ordering": "H_E^sine=(T-T^dagger)/(2i)",
        "K_definition": "K=[V,H_E^sine]",
        "epsilon_assembly": "sum_r (-1)^r sum_perm sgn(perm) Tr_aux[C_a(K) C_b(K) C_c(V)]",
        "ordered_term_count": len(ordered_rows),
        "full_output_support": len(total),
        "full_output_norm": full_norm,
        "full_scalar_diagnostics": full_scalar,
        "max_spin_reached": max_reached_spin,
        "max_diagnostics": global_diag,
        "historical_primitive_charge_basis_diagnostic": {
            "value": global_diag["CK_complete_charge_basis_leakage"],
            "hard_acceptance": False,
            "reason": (
                "Frozen sine-C(K) convention: fixed-index primitive branches precede "
                "the complete gauge-invariant H_E sum."
            ),
        },
        "reverse_projection": reverse,
        "logical_projection": {
            "support": len(logical),
            "norm": logical_norm,
            "fraction_of_full_norm": logical_norm / max(full_norm, 1e-300),
            "nonzero_detection_threshold": NONZERO_TOL,
            "nonzero": logical_nonzero,
            "initial_state_return_amplitude": [
                float(initial_return.real), float(initial_return.imag)
            ],
            "amplitudes": logical_rows,
        },
        "hard_integrity_checks": hard_checks,
        "ordered_terms": ordered_rows,
        "runtime_exact_cache": cache_info,
        "scientific_implication": implication,
        "claim_boundary": (
            "This is a finite raw Lorentzian constraint-amplitude witness. No overall "
            "Lorentzian normalization or (1+beta^2) factor is inserted, because those "
            "cannot change structural zero/nonzero. The result is not a physical time "
            "propagator, not Gamma_scalar(omega,k), and not evidence for dark matter or dark energy."
        ),
        "next_if_nonzero": (
            "Compute the full 32x32 P_logical H_L P_logical matrix with the same frozen "
            "ordering, audit its Hermitian completion/S4 decomposition, then construct the "
            "declared finite constraint family and master projector before source dressing."
        ),
        "next_if_zero": (
            "Evaluate all remaining 31 logical columns before making any operator-level zero claim."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-node", type=int, default=0)
    ap.add_argument("--input-index", type=int, default=0)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.source_node, args.input_index)
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
