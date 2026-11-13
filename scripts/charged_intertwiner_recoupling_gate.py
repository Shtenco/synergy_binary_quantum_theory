#!/usr/bin/env python3
"""Open-charge SU(2) intertwiner blocks needed by h[h^-1,K] Lorentzian legs.

After one fundamental holonomy hit, each affected four-valent geometric node
carries total J=1/2; an external fundamental closes it to a five-valent singlet.
This gate constructs the exact recoupling basis

    ((j1 j2)->K12, (j3 j4)->K34) -> J,M,

verifies its orthonormality, and constructs the volume in the symmetry-adapted
block rather than taking a global magnetic-space square root.

The oriented seed

    Q = J1 . (J2 x J3)

is an SU(2) scalar and must preserve every fixed-(J,M) sector.  We first verify
that statement directly to machine precision and verify that its reduced block
is independent of M.  A single canonical reduced block is then used for the
whole J=1/2 multiplet:

    V_J = sqrt(|Q_J|).

The invariant functional-calculus check is V_J^4 = Q_J^2.  We deliberately do
not compare two independently diagonalised copies of V_J for M=+/-1/2 because
arbitrary eigensolver rotations inside exact-zero eigenspaces can produce a
spurious O(sqrt(machine epsilon)) matrix difference even when Q_J is identical.

This is the compressed representation primitive needed by
C_e(K)=h_e[h_e^-1,K], not H_L or HDA closure.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k5_peter_weyl_safe_hda_column as PW
import bqg_pruning_audit as PRUNE


def allowed_pair2(sa, sb):
    return tuple(range(abs(sa - sb), sa + sb + 1, 2))


def allowed_charged_labels(spins, J2):
    s1, s2, s3, s4 = spins
    out = []
    for K12 in allowed_pair2(s1, s2):
        for K34 in allowed_pair2(s3, s4):
            if abs(K12 - K34) <= J2 <= K12 + K34 and (K12 + K34 + J2) % 2 == 0:
                out.append((K12, K34))
    return tuple(out)


@PRUNE.audit_lru_cache(None)
def charged_tensor(spins, K12, K34, J2, M2):
    """Build the coupled four-leg tensor by a separable CG contraction.

    The previous implementation evaluated the same three Clebsch-Gordan
    factors inside every magnetic-index quadruple.  Precomputing the two
    pair-coupling tensors and the final (K12,K34)->J coupling turns that
    repeated Python loop into one small optimized tensor contraction.
    """
    spins = tuple(spins)
    s1, s2, s3, s4 = spins
    m1, m2, m3, m4 = (PW.m2vals_t(s) for s in spins)
    mk, ml = PW.m2vals_t(K12), PW.m2vals_t(K34)

    c12 = np.zeros((len(m1), len(m2), len(mk)), dtype=float)
    c34 = np.zeros((len(m3), len(m4), len(ml)), dtype=float)
    ctotal = np.zeros((len(mk), len(ml)), dtype=float)

    for i, a in enumerate(m1):
        for j, b in enumerate(m2):
            for k, x in enumerate(mk):
                if a + b == x:
                    c12[i, j, k] = PW.cg2(s1, s2, K12, a, b, x)

    for i, a in enumerate(m3):
        for j, b in enumerate(m4):
            for k, x in enumerate(ml):
                if a + b == x:
                    c34[i, j, k] = PW.cg2(s3, s4, K34, a, b, x)

    for i, x in enumerate(mk):
        for j, y in enumerate(ml):
            if x + y == M2:
                ctotal[i, j] = PW.cg2(K12, K34, J2, x, y, M2)

    return np.einsum(
        "ijx,kly,xy->ijkl", c12, c34, ctotal, optimize=True
    ).astype(complex, copy=False)


def close_with_external_half(TJM, J2, M2):
    if J2 != 1:
        raise ValueError("a single external fundamental closes only J=1/2")
    out = np.zeros(TJM.shape + (2,), complex)
    for ie, me in enumerate(PW.m2vals_t(1)):
        c = PW.cg2(J2, 1, 0, M2, me, 0)
        if c:
            out[..., ie] += c * TJM
    return out


@PRUNE.audit_lru_cache(None)
def q123_matrix(spins3):
    spins3 = tuple(spins3)
    mats = [PW.spin_mats_cached(s) for s in spins3]
    d = int(np.prod([s + 1 for s in spins3]))
    Q = np.zeros((d, d), complex)
    for a, b, c in itertools.product(range(3), repeat=3):
        e = PW.EPS3[a, b, c]
        if e:
            Q += e * np.kron(np.kron(mats[0][a], mats[1][b]), mats[2][c])
    return 0.5 * (Q + Q.conj().T)


def apply_q_four(T, spins):
    d1, d2, d3, d4 = [s + 1 for s in spins]
    Q = q123_matrix(tuple(spins[:3]))
    return (Q @ T.reshape(d1 * d2 * d3, d4)).reshape(d1, d2, d3, d4)


def project_to_basis(X, basis):
    if not basis:
        return np.zeros(0, dtype=complex), np.zeros_like(X)
    # Row-major basis storage turns the magnetic-index projection and
    # reconstruction into two BLAS-backed contractions instead of Python
    # loops over every recoupling channel.
    basis_matrix = np.stack([np.asarray(B).reshape(-1) for B in basis], axis=0)
    flat = np.asarray(X).reshape(-1)
    coeff = basis_matrix.conj() @ flat
    recon = (coeff @ basis_matrix).reshape(np.asarray(X).shape)
    return np.asarray(coeff, dtype=complex), recon


def q_block(spins, J2, M2):
    labels = allowed_charged_labels(spins, J2)
    basis = [charged_tensor(tuple(spins), a, b, J2, M2) for a, b in labels]
    n = len(basis)
    Qb = np.zeros((n, n), complex)
    max_abs_leak = 0.0
    max_rel_leak = 0.0
    for j, B in enumerate(basis):
        QB = apply_q_four(B, spins)
        coeff, recon = project_to_basis(QB, basis)
        Qb[:, j] = coeff
        abs_leak = float(np.linalg.norm(QB - recon))
        rel_leak = abs_leak / max(float(np.linalg.norm(QB)), 1e-30)
        max_abs_leak = max(max_abs_leak, abs_leak)
        max_rel_leak = max(max_rel_leak, rel_leak)
    Qb = 0.5 * (Qb + Qb.conj().T)
    return basis, Qb, max_abs_leak, max_rel_leak


def canonical_volume_block(Qb):
    ev, U = np.linalg.eigh(Qb)
    Vb = (U * np.sqrt(np.abs(ev))) @ U.conj().T
    Vb = 0.5 * (Vb + Vb.conj().T)
    return Vb


def run():
    one_hit = []
    base = [1, 1, 1, 1]
    for leg in range(4):
        for so in (0, 2):
            q = base.copy(); q[leg] = so
            one_hit.append(tuple(q))
    one_hit = tuple(dict.fromkeys(one_hit))

    max_orth = 0.0
    max_closed_orth = 0.0
    max_q_abs_leak = 0.0
    max_q_rel_leak = 0.0
    max_q_herm = 0.0
    max_v_herm = 0.0
    max_q_m_difference = 0.0
    max_volume_function_error = 0.0
    min_v_eig = float("inf")
    rows = []

    for spins in one_hit:
        labels = allowed_charged_labels(spins, 1)
        qblocks = {}
        basis_by_M = {}
        row_q_abs = row_q_rel = row_q_herm = 0.0
        for M2 in PW.m2vals_t(1):
            basis, Qb, q_abs, q_rel = q_block(spins, 1, M2)
            basis_by_M[M2] = basis
            qblocks[M2] = Qb
            G = np.array([[np.vdot(A, B) for B in basis] for A in basis], complex)
            max_orth = max(max_orth, float(np.linalg.norm(G - np.eye(len(basis)))))
            qh = float(np.linalg.norm(Qb - Qb.conj().T))
            row_q_abs = max(row_q_abs, q_abs); row_q_rel = max(row_q_rel, q_rel)
            row_q_herm = max(row_q_herm, qh)
            max_q_abs_leak = max(max_q_abs_leak, q_abs)
            max_q_rel_leak = max(max_q_rel_leak, q_rel)
            max_q_herm = max(max_q_herm, qh)

        Mp, Mm = PW.m2vals_t(1)
        qdiff = float(np.linalg.norm(qblocks[Mp] - qblocks[Mm]))
        max_q_m_difference = max(max_q_m_difference, qdiff)
        # Average only machine-roundoff differences between two proven copies
        # of the same scalar reduced operator; this is not a physics fit.
        Qcanon = 0.5 * (qblocks[Mp] + qblocks[Mm])
        Qcanon = 0.5 * (Qcanon + Qcanon.conj().T)
        Vcanon = canonical_volume_block(Qcanon)
        vh = float(np.linalg.norm(Vcanon - Vcanon.conj().T))
        max_v_herm = max(max_v_herm, vh)
        V2 = Vcanon @ Vcanon
        V4 = V2 @ V2
        Q2 = Qcanon @ Qcanon
        func_err = float(np.linalg.norm(V4 - Q2) / max(np.linalg.norm(Q2), 1.0))
        max_volume_function_error = max(max_volume_function_error, func_err)
        if len(Vcanon):
            min_v_eig = min(min_v_eig, float(np.linalg.eigvalsh(Vcanon).min()))

        closed = []
        for ilabel, _ in enumerate(labels):
            S = np.zeros(tuple(s + 1 for s in spins) + (2,), complex)
            for M2 in PW.m2vals_t(1):
                S += close_with_external_half(basis_by_M[M2][ilabel], 1, M2)
            closed.append(S)
        G5 = np.array([[np.vdot(A, B) for B in closed] for A in closed], complex)
        err5 = float(np.linalg.norm(G5 - np.eye(len(closed)))) if closed else 0.0
        max_closed_orth = max(max_closed_orth, err5)

        rows.append({
            "spins": [s / 2 for s in spins],
            "charged_J": 0.5,
            "charged_dimension": len(labels),
            "labels_doubled": [list(x) for x in labels],
            "five_valent_singlet_orth_error": err5,
            "max_Q_absolute_sector_leakage": row_q_abs,
            "max_Q_relative_sector_leakage": row_q_rel,
            "max_Q_block_hermiticity_error": row_q_herm,
            "Q_block_M_difference": qdiff,
            "V_block_hermiticity_error": vh,
            "V4_minus_Q2_relative_error": func_err,
            "Q_canonical_real": Qcanon.real.tolist(),
            "Q_canonical_imag": Qcanon.imag.tolist(),
            "V_canonical_real": Vcanon.real.tolist(),
            "V_canonical_imag": Vcanon.imag.tolist(),
        })

    dims = sorted({r["charged_dimension"] for r in rows})
    if min_v_eig == float("inf"):
        min_v_eig = 0.0
    passed = (
        len(one_hit) == 8
        and dims == [2, 3]
        and max_orth < 1e-11
        and max_closed_orth < 1e-11
        and max_q_abs_leak < 1e-12
        and max_q_rel_leak < 1e-12
        and max_q_herm < 1e-12
        and max_q_m_difference < 1e-12
        and max_v_herm < 1e-12
        and max_volume_function_error < 1e-12
        and min_v_eig > -1e-12
    )
    return {
        "status": "canonical symmetry-adapted charged-intertwiner Q/V blocks for Lorentzian holonomy commutators",
        "passed": bool(passed),
        "one_fundamental_hit_quartets": len(one_hit),
        "charged_dimensions_observed": dims,
        "rows": rows,
        "max_four_leg_orthonormality_error": max_orth,
        "max_five_valent_singlet_orthonormality_error": max_closed_orth,
        "max_Q_absolute_JM_sector_leakage": max_q_abs_leak,
        "max_Q_relative_JM_sector_leakage": max_q_rel_leak,
        "max_Q_block_hermiticity_error": max_q_herm,
        "max_Q_M_block_difference": max_q_m_difference,
        "max_V_block_hermiticity_error": max_v_herm,
        "max_V4_minus_Q2_relative_error": max_volume_function_error,
        "minimum_V_block_eigenvalue": min_v_eig,
        "representation_statement": (
            "A one-hit endpoint is exactly compressed into a total-J=1/2 four-leg charged block of dimension 2 or 3; the external fundamental closes it to a five-valent singlet."
        ),
        "volume_statement": (
            "Q is first reduced as an SU(2) scalar and verified M-independent. One canonical V_J=sqrt(|Q_J|) is then defined for the entire multiplet and satisfies V_J^4=Q_J^2."
        ),
        "next_use": (
            "Use the canonical charged V_J blocks inside the charged K_v action and close the spectator fundamental with h_e to compute C_e(K_v)=h_e[h_e^-1,K_v]."
        ),
        "scope_note": "Exact representation/volume prerequisite only; C_e(K), the Lorentzian triple and HDA remain open.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    a = ap.parse_args()
    out = run(); text = json.dumps(out, indent=2); print(text)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
