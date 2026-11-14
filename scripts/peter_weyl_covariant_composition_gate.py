#!/usr/bin/env python3
"""Composition gate for matrix-covariant Peter-Weyl legs.

The existing C_e(V) and C_e(K) gates construct exact columns on a Gauss input.
A Lorentzian triple requires something stronger: after one covariant leg the
source carries J=0+1 and the next leg must act on that state without projecting
it back to J=0.

This gate implements that missing representation bridge for C_e(V):

  covariant source key
    -> complete all-J labels at source
    -> h^{-1} creates a second charged endpoint
    -> V acts in the exact source-J block
    -> h closes the target endpoint
    -> source is retained in every allowed J sector.

Acceptance is deliberately independent of H_L:
1. on the frozen Gauss input the generalized implementation must reproduce the
   independently existing C_e(V) column component-by-component;
2. sum_k h_ik h^{-1}_kj must equal delta_ij on a genuinely non-Gauss J=1
   intermediate source state;
3. a second C(V) action on that J=1 state must stay inside the SU(2) selection
   rule J_out subset {0,1,2} and be nonzero.

No Lorentzian coefficient, beta, HDA normalization or tolerance is fitted.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW
import bqg_pruning_audit as PRUNE
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_lorentzian_K_block_gate as KG

TOL = 1e-11


def add(dst, src, scale=1.0, tol=TOL):
    """Coherently merge complex amplitudes without an intermediate cutoff.

    tol is retained for audit-event classification only. Applying it to this
    incremental merge would make the operator nonlinear and order dependent:
    a sub-threshold partial sum could be deleted before later contributions
    arrive. Only an exact floating-point zero is removed here; physical and
    regulator truncations must be explicit, separately certified boundaries
    after the complete coherent sum.
    """
    factor = complex(scale)
    if not (math.isfinite(factor.real) and math.isfinite(factor.imag)):
        raise ValueError("merge scale must be finite")
    if not math.isfinite(float(tol)) or float(tol) < 0.0:
        raise ValueError("merge audit tolerance must be finite and nonnegative")
    for k, a in src.items():
        prior = complex(dst.get(k, 0j))
        incoming = factor * complex(a)
        real = math.fsum((prior.real, incoming.real))
        imag = math.fsum((prior.imag, incoming.imag))
        z = complex(real, imag)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError("merged complex amplitude must be finite")
        PRUNE.record_merge("covariant_composition.add", k, prior, incoming, z, tol)
        if z != 0j:
            dst[k] = z
        else:
            dst.pop(k, None)


def covariant_key_to_complete(key, source_v):
    spins, Kother, J2, M2, K12, K34 = key
    labels = []
    for u in PW.VERT:
        if u == source_v:
            labels.append(CK.label_charge(J2, M2, K12, K34))
        else:
            K = Kother[u]
            if K < 0:
                raise ValueError("only source_v may carry the covariant sentinel")
            labels.append(CK.label_gauss(K))
    return spins, tuple(labels)


def complete_to_covariant(state, source_v, tol=TOL):
    """Change representation labels without pruning coherent amplitudes.

    The tolerance remains a validation/diagnostic parameter only.  This map is
    a basis relabelling after the target endpoint has closed, not a license to
    erase small nonzero components.
    """
    if not math.isfinite(float(tol)) or float(tol) < 0.0:
        raise ValueError("projection audit tolerance must be finite and nonnegative")
    terms = {}
    for (spins, labels), amp in state.items():
        src = labels[source_v]
        if src[0] != "C":
            raise ValueError("source must remain in an explicit total-J sector")
        _, J2, M2, K12, K34 = src
        Kother = []
        for u, lab in enumerate(labels):
            if u == source_v:
                Kother.append(-1)
            else:
                if lab[0] != "G":
                    raise ValueError("target charge must close before covariant projection")
                Kother.append(lab[1])
        key = (spins, tuple(Kother), J2, M2, K12, K34)
        terms.setdefault(key, []).append(complex(amp))
    out = {}
    for key, values in terms.items():
        z = complex(
            math.fsum(value.real for value in values),
            math.fsum(value.imag for value in values),
        )
        if z != 0j:
            out[key] = z
        else:
            PRUNE.record_merge(
                "covariant_composition.complete_to_covariant", key, 0j, 0j, z, tol
            )
    return out


def covariant_to_complete(state, source_v):
    out = {}
    for key, amp in state.items():
        ko = covariant_key_to_complete(key, source_v)
        out[ko] = out.get(ko, 0j) + amp
    return out


def inverse_complete(state_cov, source_v, target_v, k, j, Jmax2):
    out = {}
    max_leak = 0.0
    for key, amp in state_cov.items():
        br = CK.branch_from_key(covariant_key_to_complete(key, source_v), amp)
        for hb in PW.apply_hit_branch(br, target_v, source_v, k, j, Jmax2):
            projected, leak = CK.project_branch_complete_charges(
                hb, (source_v, target_v)
            )
            max_leak = max(max_leak, leak)
            add(out, projected)
    return out, max_leak


def close_complete(state_complete, source_v, target_v, i, k, Jmax2):
    out = {}
    for key, amp in state_complete.items():
        br = CK.branch_from_key(key, amp)
        for hb in PW.apply_hit_branch(br, source_v, target_v, i, k, Jmax2):
            projected = CV.project_covariant_branches([hb], source_v)
            add(out, projected)
    return out


def direct_volume_covariant(state_cov, source_v):
    complete = covariant_to_complete(state_cov, source_v)
    return complete_to_covariant(
        CK.apply_volume_complete_state(complete, source_v), source_v
    )


def hh_inverse_component(state_cov, source_v, target_v, i, j, Jmax2):
    total = {}
    max_leak = 0.0
    for k in range(2):
        inv, leak = inverse_complete(
            state_cov, source_v, target_v, k, j, Jmax2
        )
        max_leak = max(max_leak, leak)
        add(total, close_complete(inv, source_v, target_v, i, k, Jmax2))
    return total, max_leak


def C_volume_component(state_cov, source_v, target_v, i, j, Jmax2):
    direct = direct_volume_covariant(state_cov, source_v)
    hVh = {}
    max_leak = 0.0
    for k in range(2):
        inv, leak = inverse_complete(
            state_cov, source_v, target_v, k, j, Jmax2
        )
        max_leak = max(max_leak, leak)
        Vinv = CK.apply_volume_complete_state(inv, source_v)
        add(hVh, close_complete(Vinv, source_v, target_v, i, k, Jmax2))
    out = {}
    if i == j:
        add(out, direct, +1)
    add(out, hVh, -1)
    return out, max_leak


def norm2(state):
    return float(sum(abs(a)**2 for a in state.values()))


def diff_norm(a, b):
    keys = set(a) | set(b)
    return math.sqrt(sum(abs(a.get(k,0j)-b.get(k,0j))**2 for k in keys))


def relerr(a, b):
    return diff_norm(a,b) / max(math.sqrt(norm2(b)), 1e-30)


def weight_by_J(state):
    w = {}
    for key, amp in state.items():
        J2 = key[2]
        w[J2] = w.get(J2, 0.0) + abs(amp)**2
    return {str(J2/2):float(x) for J2,x in sorted(w.items())}


def reference_CV_matrix(initial, source_v, target_v, Jmax2=3):
    Vgauss = dict(KG.local_volume_column(initial, source_v))
    Vcov = CV.gauss_to_covariant(Vgauss, source_v)
    C = [[{} for _ in range(2)] for _ in range(2)]
    for i in range(2):
        for j in range(2):
            hVh, _ = CV.inverse_then_forward(
                initial, source_v, target_v, i, j, Jmax2, True
            )
            out = {}
            if i == j:
                add(out, Vcov, +1)
            add(out, hVh, -1)
            C[i][j] = out
    return C


def choose_J1_basis_state(C):
    ranked = []
    for row in C:
        for state in row:
            for key, amp in state.items():
                if key[2] == 2:  # source J=1
                    ranked.append((abs(amp), key, amp))
    if not ranked:
        raise RuntimeError("reference C(V) has no J=1 state")
    ranked.sort(reverse=True, key=lambda x:x[0])
    _, key, amp = ranked[0]
    # Normalize a single exact basis vector.  Its original coefficient is
    # reported but does not affect the operator identity checks below.
    return {key: 1+0j}, key, amp


def run(source_v=0, target_v=1):
    initial = PW.basis_full_jhalf()[0]
    gauss_cov = CV.gauss_to_covariant({initial:1+0j}, source_v)
    ref = reference_CV_matrix(initial, source_v, target_v, 3)

    reproduction_errors = []
    reproduction_support = []
    max_gauss_leak = 0.0
    generalized = [[{} for _ in range(2)] for _ in range(2)]
    for i in range(2):
        for j in range(2):
            got, leak = C_volume_component(
                gauss_cov, source_v, target_v, i, j, 3
            )
            generalized[i][j] = got
            max_gauss_leak = max(max_gauss_leak, leak)
            reproduction_errors.append(relerr(got, ref[i][j]))
            reproduction_support.append((len(got), len(ref[i][j])))

    J1state, J1key, original_amp = choose_J1_basis_state(ref)
    # Two hits can raise the already reached j<=3/2 state by another unit, so
    # Jmax=5/2 is a conservative exact wall for this composition identity.
    JMAX2 = 5
    identity_errors = []
    max_J1_identity_leak = 0.0
    for i in range(2):
        for j in range(2):
            got, leak = hh_inverse_component(
                J1state, source_v, target_v, i, j, JMAX2
            )
            max_J1_identity_leak = max(max_J1_identity_leak, leak)
            target = J1state if i == j else {}
            identity_errors.append(relerr(got, target) if target else diff_norm(got, target))

    second = {}
    max_second_leak = 0.0
    for i in range(2):
        for j in range(2):
            got, leak = C_volume_component(
                J1state, source_v, target_v, i, j, JMAX2
            )
            max_second_leak = max(max_second_leak, leak)
            add(second, got)
    second_weights = weight_by_J(second)
    bad_second_weight = sum(
        x for J,x in ((float(k),v) for k,v in second_weights.items())
        if J not in (0.0,1.0,2.0)
    )
    second_norm = math.sqrt(norm2(second))
    max_spin_second = max((max(key[0]) for key in second), default=0)/2

    passed = (
        max(reproduction_errors, default=0.0) < 1e-10
        and all(a == b for a,b in reproduction_support)
        and max_gauss_leak < 1e-10
        and max(identity_errors, default=0.0) < 1e-10
        and max_J1_identity_leak < 1e-10
        and second_norm > 1e-10
        and bad_second_weight < 1e-18
        and max_second_leak < 1e-10
        and max_spin_second <= 2.5 + 1e-12
    )
    return {
        "status":"generalized covariant Peter-Weyl composition gate for C_e(V)",
        "passed":bool(passed),
        "edge":[source_v,target_v],
        "gauss_column_component_relative_errors":reproduction_errors,
        "gauss_column_component_support_pairs":reproduction_support,
        "max_gauss_input_complete_basis_leakage":max_gauss_leak,
        "selected_J1_basis_key":repr(J1key),
        "selected_J1_reference_amplitude":[original_amp.real,original_amp.imag],
        "J1_two_hit_identity_component_errors":identity_errors,
        "J1_two_hit_identity_max_complete_basis_leakage":max_J1_identity_leak,
        "second_CV_combined_norm":second_norm,
        "second_CV_weight_by_source_J":second_weights,
        "second_CV_forbidden_J_weight":bad_second_weight,
        "second_CV_max_complete_basis_leakage":max_second_leak,
        "second_CV_max_spin_reached":max_spin_second,
        "selection_rule":"J_in=1 acted on by a fundamental 2x2 covariant operator may produce only J_out=0,1,2 after target closure.",
        "next_use":"Generalize the same state-to-state bridge to C_e(K), then evaluate the traced two-K one-V Lorentzian triple without multiplying precomputed columns.",
        "scope_note":"Representation-composition gate only; it is not yet H_L and not an HDA closure claim."
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v",type=int,default=0)
    ap.add_argument("--w",type=int,default=1)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args(); out=run(a.v,a.w); text=json.dumps(out,indent=2); print(text)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text+"\n",encoding="utf-8")
    return 0 if out["passed"] else 1

if __name__=="__main__":
    raise SystemExit(main())
