#!/usr/bin/env python3
"""Preregistered zero-aware volume migration experiment.

Purpose
-------
Test the diagnosis from peter_weyl_volume_cross_rep_audit without weakening any
existing H_E/C(K) acceptance threshold.  The only experimental change is to
map eigenvalues compatible with exact zero under a conservative backward-error
bound to zero *before* applying sqrt(abs(.)):

    tau(Q) = 1000 eps dim(Q) max(1, ||Q||_spectral),
    V      = U diag(sqrt(|lambda|) if |lambda| > tau else 0) U^dagger.

This is not an eigenvalue fit: the preregistered cross-representation audit
already requires the smallest nonzero |lambda| to be at least 1e8 tau on every
H_E-reached local quartet.  The experiment monkey-patches both independently
existing magnetic and symmetry-adapted volume implementations in one fresh
Python process, clears affected caches, then reruns the frozen invariant
covariant-K audit.  Production files are deliberately left untouched until
this experiment is green.
"""
from __future__ import annotations

import functools
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import bqg_pruning_audit as PRUNE

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import charged_intertwiner_recoupling_gate as CH
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_K_projection_audit_gate as AUD

EPS = np.finfo(float).eps
ZERO_FACTOR = 1000.0
SPECTRAL_AUDIT = {}


def zeroaware_sqrt_abs(H):
    H = 0.5 * (H + H.conj().T)
    ev, U = np.linalg.eigh(H)
    scale = float(np.max(np.abs(ev))) if len(ev) else 0.0
    tau = ZERO_FACTOR * EPS * H.shape[0] * max(1.0, scale)
    zeroed = np.abs(ev) <= tau
    retained = np.abs(ev) > tau
    vals = np.where(retained, np.sqrt(np.abs(ev)), 0.0)
    V = (U * vals) @ U.conj().T
    dim_key = str(int(H.shape[0]))
    row = SPECTRAL_AUDIT.setdefault(dim_key, {
        "dimension": int(H.shape[0]), "matrix_count": 0,
        "eigenvalue_count": 0, "zeroed_eigenvalue_count": 0,
        "retained_eigenvalue_count": 0, "max_zeroed_abs": 0.0,
        "min_retained_abs_over_tau": None, "max_tau": 0.0,
    })
    row["matrix_count"] += 1
    row["eigenvalue_count"] += int(len(ev))
    row["zeroed_eigenvalue_count"] += int(np.count_nonzero(zeroed))
    row["retained_eigenvalue_count"] += int(np.count_nonzero(retained))
    row["max_zeroed_abs"] = max(
        float(row["max_zeroed_abs"]),
        float(np.max(np.abs(ev[zeroed]))) if np.any(zeroed) else 0.0,
    )
    if np.any(retained):
        gap = float(np.min(np.abs(ev[retained])) / max(tau, np.finfo(float).tiny))
        current = row["min_retained_abs_over_tau"]
        row["min_retained_abs_over_tau"] = gap if current is None else min(float(current), gap)
    row["max_tau"] = max(float(row["max_tau"]), float(tau))
    return 0.5 * (V + V.conj().T), tau, ev


def spectral_audit_snapshot():
    return {
        "zero_factor": float(ZERO_FACTOR),
        "machine_epsilon": float(EPS),
        "by_dimension": {key: dict(value) for key, value in sorted(SPECTRAL_AUDIT.items())},
        "matrix_count": int(sum(row["matrix_count"] for row in SPECTRAL_AUDIT.values())),
        "zeroed_eigenvalue_count": int(sum(row["zeroed_eigenvalue_count"] for row in SPECTRAL_AUDIT.values())),
        "retained_eigenvalue_count": int(sum(row["retained_eigenvalue_count"] for row in SPECTRAL_AUDIT.values())),
    }


@PRUNE.audit_lru_cache(2)
def zeroaware_volume123_matrix(s1, s2, s3):
    mats = [PW.spin_mats_cached(s) for s in (s1, s2, s3)]
    d = (s1 + 1) * (s2 + 1) * (s3 + 1)
    Q = np.zeros((d, d), complex)
    for a, b, c in itertools.product(range(3), repeat=3):
        e = PW.EPS3[a, b, c]
        if e:
            Q += e * np.kron(np.kron(mats[0][a], mats[1][b]), mats[2][c])
    V, _, _ = zeroaware_sqrt_abs(Q)
    return V


def zeroaware_canonical_volume_block(Qb):
    V, _, _ = zeroaware_sqrt_abs(Qb)
    return V


def patch_and_clear(zero_factor=None):
    global ZERO_FACTOR, SPECTRAL_AUDIT
    if zero_factor is not None:
        zero_factor = float(zero_factor)
        if not math.isfinite(zero_factor) or zero_factor <= 0:
            raise ValueError("zero_factor must be finite and positive")
        ZERO_FACTOR = zero_factor
    SPECTRAL_AUDIT = {}
    # Both sides of the independent equivalence test receive the same
    # mathematically defined zero-eigenspace convention.
    PW.volume123_matrix = zeroaware_volume123_matrix
    CH.canonical_volume_block = zeroaware_canonical_volume_block
    zeroaware_volume123_matrix.cache_clear()
    # The general-J cache stores V blocks produced through CH.
    CK.canonical_volume_block_general.cache_clear()
    # Clear higher-level cached columns so no pre-patch value can survive.
    for obj in (getattr(PW, 'apply_H_cached', None),
                getattr(PW, 'apply_H_cached_state', None),
                getattr(CK, 'HE_complete_cached', None)):
        if hasattr(obj, 'cache_clear'):
            obj.cache_clear()


def run():
    patch_and_clear()
    audit = AUD.run()
    va = audit['volume_cross_representation_audit']
    he_ok = bool(audit.get('production_HE_equivalent', False))
    ck_ok = bool(audit.get('passed', False))
    min_gap = float(va.get('minimum_Q_gap_over_zero_tolerance', 0.0))
    return {
        'status': 'zero-aware Peter-Weyl volume migration experiment',
        'passed': bool(he_ok and ck_ok and min_gap >= 1e8),
        'zero_factor': ZERO_FACTOR,
        'machine_epsilon': EPS,
        'minimum_nonzero_Q_gap_over_zero_tolerance': min_gap,
        'frozen_HE_relative_error_threshold': 1e-9,
        'production_HE_equivalent': he_ok,
        'production_HE_relative_error': audit.get('production_relative_error'),
        'production_HE_support_identical': audit.get('production_support_identical'),
        'max_excluded_tail_amplitude': audit.get('max_excluded_tail_amplitude'),
        'allJ_internal_volume_sector_leakage': audit.get('allJ_internal_volume_sector_leakage'),
        'covariant_K_audit_passed': ck_ok,
        'C_matrix_Frobenius_covariant_state_norm': audit.get('C_matrix_Frobenius_covariant_state_norm'),
        'C_weight_by_source_J': audit.get('C_weight_by_source_J'),
        'C_J_greater_than_1_weight_fraction': audit.get('C_J_greater_than_1_weight_fraction'),
        'max_spin_reached': audit.get('max_spin_reached'),
        'volume_cross_representation_audit_passed': va.get('passed'),
        'max_Vmag_vs_Vblock_relative_frobenius': va.get('max_Vmag_vs_Vblock_relative_frobenius'),
        'max_Vmag_kernel_action_frobenius': va.get('max_Vmag_kernel_action_frobenius'),
        'max_Vblock_kernel_action_frobenius': va.get('max_Vblock_kernel_action_frobenius'),
        'max_Vzeroaware_kernel_action_frobenius': va.get('max_Vzeroaware_kernel_action_frobenius'),
        'q_block_gauge_audit_passed_as_diagnosis': audit['q_block_gauge_audit'].get('passed'),
        'q_block_max_fixedJ_M_difference': audit['q_block_gauge_audit'].get('max_fixedJ_M_block_difference'),
        'scope': 'single frozen H_E column plus existing finite C_e(K) invariant audit; H_L triple and joint Lorentzian HDA remain outside this experiment',
        'full_audit': audit,
    }


def main():
    out = run()
    print(json.dumps(out, indent=2))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
