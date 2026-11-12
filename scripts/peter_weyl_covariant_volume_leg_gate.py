#!/usr/bin/env python3
"""Exact matrix-covariant Lorentzian volume leg C_e(V)=h_e[h_e^-1,V_v].

A crucial representation point is kept explicit: C_e(V) is a 2x2 matrix at the
source vertex, not a Gauss-scalar by itself.  After h^-1 and h close the edge at
the target vertex, the source geometry is correlated with the open fundamental
matrix pair

    1/2 tensor 1/2* = J=0 plus J=1.

Therefore the final source tensor is NOT projected to J=0.  It is decomposed
into exact total-J recoupling tensors, while every other K5 node is projected to
its ordinary Gauss singlet.  Only the later trace of three covariant legs is a
full scalar.

Before evaluating C(V), the script checks the stronger two-hit identity

    sum_k h_ik h^-1_kj = delta_ij

in this enlarged covariant basis.  The identity must contain only the J=0
source sector and reproduce the original Gauss state.  The nontrivial C(V) leg
may contain J=0 and J=1, but any J>1 content is a hard covariance failure.

At the intermediate one-hit source state, V acts in the symmetry-adapted
charged J=1/2 blocks constructed by charged_intertwiner_recoupling_gate.py.

This is one genuine factor of Thiemann's Lorentzian kinetic operator; K legs,
the triple trace and HDA remain separate gates.
"""
from __future__ import annotations

import argparse
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
import bqg_pruning_audit as PRUNE
import charged_intertwiner_recoupling_gate as CH
import peter_weyl_lorentzian_K_block_gate as KG


def unorient_local(T, spins_local, v):
    X = T
    for leg, w in enumerate(PW.NEIG[v]):
        if w < v:
            X = PW.apply_axis_np(X, leg, PW.epsilon_j(spins_local[leg]).conj().T)
    return X


def orient_local(T, spins_local, v):
    X = T
    for leg, w in enumerate(PW.NEIG[v]):
        if w < v:
            X = PW.apply_axis_np(X, leg, PW.epsilon_j(spins_local[leg]))
    return X


def canonical_charged_volume(spins_local):
    Mp, Mm = PW.m2vals_t(1)
    _, Qp, _, _ = CH.q_block(tuple(spins_local), 1, Mp)
    _, Qm, _, _ = CH.q_block(tuple(spins_local), 1, Mm)
    Q = 0.5 * (Qp + Qm)
    Q = 0.5 * (Q + Q.conj().T)
    return CH.canonical_volume_block(Q)


def apply_charged_volume_oriented(T, spins_local, v):
    """Apply V on the exact one-hit J=1/2 source sector."""
    X = unorient_local(T, spins_local, v)
    V = canonical_charged_volume(tuple(spins_local))
    Y = np.zeros_like(X)
    recon = np.zeros_like(X)
    labels = CH.allowed_charged_labels(tuple(spins_local), 1)
    for M2 in PW.m2vals_t(1):
        basis = [CH.charged_tensor(tuple(spins_local), a, b, 1, M2) for a, b in labels]
        coeff = np.asarray([np.vdot(B, X) for B in basis], complex)
        for c, B in zip(coeff, basis):
            recon += c * B
        outc = V @ coeff
        for c, B in zip(outc, basis):
            Y += c * B
    leak = float(np.linalg.norm(X - recon) / max(np.linalg.norm(X), 1e-30))
    return orient_local(Y, spins_local, v), leak


def all_total_J2(spins_local):
    maxJ = sum(spins_local)
    return tuple(J for J in range(maxJ + 1) if CH.allowed_charged_labels(tuple(spins_local), J))


def project_covariant_branches(branches, source_v, tol=1e-11):
    """Project nonsource nodes to J=0 and source to every allowed total J.

    Key layout:
      (spins, K_other_tuple, J2, M2, K12, K34)
    where K_other_tuple has -1 at source_v and Gauss recoupling labels elsewhere.
    """
    out = {}
    for spins, tensors, amp in branches:
        other_opts = []
        ok = True
        for u in PW.VERT:
            if u == source_v:
                other_opts.append(((None, 1 + 0j),))
                continue
            ls = PW.local_spins(spins, u)
            opts = []
            for K in PW.allowed_k2_t(*ls):
                c = np.vdot(PW.oriented_intertwiner(u, ls, K), tensors[u])
                if abs(c) > 1e-13:
                    opts.append((K, c))
                else:
                    PRUNE.record_threshold("covariant_volume_leg.gauss_projection.1e-13",(spins,u,K),c,1e-13)
            if not opts:
                ok = False
                break
            other_opts.append(tuple(opts))
        if not ok:
            continue

        ls0 = PW.local_spins(spins, source_v)
        X0 = unorient_local(tensors[source_v], ls0, source_v)
        src_opts = []
        for J2 in all_total_J2(ls0):
            for M2 in PW.m2vals_t(J2):
                for K12, K34 in CH.allowed_charged_labels(tuple(ls0), J2):
                    B = CH.charged_tensor(tuple(ls0), K12, K34, J2, M2)
                    c = np.vdot(B, X0)
                    if abs(c) > 1e-13:
                        src_opts.append((J2, M2, K12, K34, c))
                    else:
                        PRUNE.record_threshold("covariant_volume_leg.charged_projection.1e-13",(spins,u,J2,M2,K12,K34),c,1e-13)

        for chosen in itertools.product(*other_opts):
            base_amp = amp
            Kother = []
            for u, (K, c) in enumerate(chosen):
                if u == source_v:
                    Kother.append(-1)
                else:
                    Kother.append(K)
                    base_amp *= c
            for J2, M2, K12, K34, cs in src_opts:
                val = base_amp * cs
                key = (spins, tuple(Kother), J2, M2, K12, K34)
                if abs(val) > tol:
                    out[key] = out.get(key, 0j) + val
                else:
                    PRUNE.record_threshold("covariant_volume_leg.branch_projection",key,val,tol)
    for key, value in out.items():
        if abs(value) <= tol:
            PRUNE.record_threshold("covariant_volume_leg.branch_output",key,value,tol)
    return {k: v for k, v in out.items() if abs(v) > tol}


def inverse_then_forward(initial, v, w, i, j, Jmax2, with_volume):
    total = {}
    max_charged_projection_leak = 0.0
    for k in range(2):
        branches = []
        for br in [PW.initial_factorized_oriented(initial)]:
            branches.extend(PW.apply_hit_branch(br, w, v, k, j, Jmax2))
        if with_volume:
            vb = []
            for spins, tensors, amp in branches:
                t = list(tensors)
                ls = PW.local_spins(spins, v)
                t[v], leak = apply_charged_volume_oriented(t[v], ls, v)
                max_charged_projection_leak = max(max_charged_projection_leak, leak)
                vb.append((spins, tuple(t), amp))
            branches = vb
        closed = []
        for br in branches:
            closed.extend(PW.apply_hit_branch(br, v, w, i, k, Jmax2))
        projected = project_covariant_branches(closed, v)
        PW.add_dict(total, projected, +1)
    for key, value in total.items():
        if abs(value) <= 1e-10:
            PRUNE.record_threshold("covariant_volume_leg.closed_state_output.1e-10",key,value,1e-10)
    return {k: a for k, a in total.items() if abs(a) > 1e-10}, max_charged_projection_leak


def gauss_to_covariant(state, source_v):
    out = {}
    for (spins, Ks), amp in state.items():
        K = Ks[source_v]
        Kother = tuple(-1 if u == source_v else Ks[u] for u in PW.VERT)
        key = (spins, Kother, 0, 0, K, K)
        out[key] = out.get(key, 0j) + amp
    return out


def covariant_state_norm2(state):
    return float(sum(abs(a) ** 2 for a in state.values()))


def diff_norm(a, b):
    keys = set(a) | set(b)
    return math.sqrt(sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in keys))


def weight_by_J(matrix_states):
    out = {}
    for row in matrix_states:
        for state in row:
            for key, amp in state.items():
                J2 = key[2]
                out[J2] = out.get(J2, 0.0) + abs(amp) ** 2
    return {str(J2 / 2): float(v) for J2, v in sorted(out.items())}


def matrix_state_norm(matrix_states):
    return math.sqrt(sum(covariant_state_norm2(s) for row in matrix_states for s in row))


def run(v=0, w=1):
    JMAX2 = 3
    initial = PW.basis_full_jhalf()[0]
    target_identity = gauss_to_covariant({initial: 1 + 0j}, v)

    ident = [[{} for _ in range(2)] for _ in range(2)]
    max_identity_error = 0.0
    for i in range(2):
        for j in range(2):
            s, _ = inverse_then_forward(initial, v, w, i, j, JMAX2, False)
            ident[i][j] = s
            target = target_identity if i == j else {}
            max_identity_error = max(max_identity_error, diff_norm(s, target))
    ident_weights = weight_by_J(ident)
    ident_nonzero_J_weight = sum(vv for jj, vv in ident_weights.items() if abs(float(jj)) > 1e-15)

    Vgauss = dict(KG.local_volume_column(initial, v))
    Vcov = gauss_to_covariant(Vgauss, v)
    C = [[{} for _ in range(2)] for _ in range(2)]
    max_charged_leak = 0.0
    for i in range(2):
        for j in range(2):
            hVh, cleak = inverse_then_forward(initial, v, w, i, j, JMAX2, True)
            max_charged_leak = max(max_charged_leak, cleak)
            out = {}
            if i == j:
                PW.add_dict(out, Vcov, +1)
            PW.add_dict(out, hVh, -1)
            C[i][j] = {k: a for k, a in out.items() if abs(a) > 1e-10}

    Cnorm = matrix_state_norm(C)
    Cweights = weight_by_J(C)
    C_high_J_weight = sum(vv for jj, vv in Cweights.items() if float(jj) > 1.0 + 1e-15)
    C_total_weight = sum(Cweights.values())
    high_fraction = C_high_J_weight / max(C_total_weight, 1e-30)
    J1_weight = Cweights.get("1.0", 0.0)
    supports = [[len(C[i][j]) for j in range(2)] for i in range(2)]
    max_spin = max((max(key[0]) for row in C for s in row for key in s), default=0) / 2

    passed = (
        max_identity_error < 1e-10
        and ident_nonzero_J_weight < 1e-20
        and max_charged_leak < 1e-10
        and Cnorm > 1e-10
        and J1_weight > 1e-14
        and high_fraction < 1e-20
        and max_spin <= 1.5 + 1e-12
    )
    return {
        "status": "exact matrix-covariant Peter-Weyl volume leg C_e(V)=h[h^-1,V]",
        "passed": bool(passed),
        "edge": [v, w],
        "input": "all ten links j=1/2; all five K=0",
        "Jmax": 1.5,
        "two_hit_identity_max_covariant_state_error": max_identity_error,
        "two_hit_identity_weight_by_source_J": ident_weights,
        "two_hit_identity_nonzero_J_weight": ident_nonzero_J_weight,
        "charged_volume_max_Jhalf_projection_leakage": max_charged_leak,
        "C_matrix_supports": supports,
        "C_matrix_Frobenius_covariant_state_norm": Cnorm,
        "C_weight_by_source_J": Cweights,
        "C_J1_weight": J1_weight,
        "C_J_greater_than_1_weight_fraction": high_fraction,
        "max_spin_reached": max_spin,
        "definition": "C_ij(V)=delta_ij V-sum_k h_ik V h^-1_kj, retaining source J=0+1 covariant geometry",
        "representation_note": (
            "C_ij is not projected to a Gauss scalar. Its source geometry carries the same J=0+1 content as the open fundamental matrix pair; only the later traced triple is gauge scalar."
        ),
        "next_use": "Construct the same matrix-covariant leg with K_v=[V_v,H_E,v], then contract two K legs and one V leg with the frozen epsilon^{ijk} orientation trace.",
        "scope_note": "One exact Lorentzian covariant volume factor only; K legs, triple H_L and HDA remain open.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v", type=int, default=0)
    ap.add_argument("--w", type=int, default=1)
    ap.add_argument("--output", type=Path)
    a = ap.parse_args()
    out = run(a.v, a.w); text = json.dumps(out, indent=2); print(text)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
