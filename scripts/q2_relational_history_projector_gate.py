#!/usr/bin/env python3
r"""Exact finite-C8 relational-history / combined-projector positive control.

The repository already proved a no-go for treating the history shift alone as
pure gauge: untwisted averaging over U kills every nontrivial history character.
This gate shows the mathematically correct relational alternative in finite
dimension.

Introduce a separate C8 clock shift S and a q=2 system step R=J, where

    J=[[0,-1],[1,0]], J^2=-I, J^8=I.

Define the combined constraint step

    G = S_clock \otimes R_geom.

The history state

    |Psi> = 1/sqrt(8) sum_t |t> \otimes R^t |psi0>

is exactly invariant under G.  The combined rigging/group-average projector

    P_rel = 1/8 sum_tau G^tau

maps a seed |0>\otimes|psi0> onto that relational history. Conditioning on the
clock reading t gives the nontrivial system state R^t|psi0> even though the
global history is gauge invariant.

By contrast, clock-only untwisted averaging kills this nontrivial relational
history exactly. The surviving physical projector is character-correlated:
clock and system phases are paired so their product character is trivial.

R=J is used only as an exact q=2 positive-control system step. This gate does
not derive a physical clock, physical time, or the microscopic gravitational
system evolution operator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def shift(n: int) -> sp.Matrix:
    S = sp.zeros(n)
    for t in range(n):
        S[(t + 1) % n, t] = 1
    return S


def simplify_matrix(M: sp.Matrix) -> sp.Matrix:
    return M.applyfunc(lambda x: sp.simplify(sp.expand_complex(x)))


def zero(M: sp.Matrix) -> bool:
    return all(sp.simplify(sp.expand_complex(x)) == 0 for x in M)


def vec_norm2(v: sp.Matrix):
    return sp.simplify((v.H * v)[0])


def projector_character(S: sp.Matrix, m: int) -> sp.Matrix:
    n = S.rows
    zeta = sp.exp(2 * sp.pi * sp.I / n)
    out = sp.zeros(n)
    St = sp.eye(n)
    for tau in range(n):
        out += zeta ** (-m * tau) * St
        St *= S
    return simplify_matrix(out / n)


def group_average(G: sp.Matrix, order: int) -> sp.Matrix:
    out = sp.zeros(G.rows)
    Gt = sp.eye(G.rows)
    for _ in range(order):
        out += Gt
        Gt *= G
    return simplify_matrix(out / order)


def conditional_system_state(Psi: sp.Matrix, t: int, dsys: int = 2) -> sp.Matrix:
    return sp.Matrix(Psi[t * dsys:(t + 1) * dsys, 0])


def normalized(v: sp.Matrix) -> sp.Matrix:
    n2 = sp.simplify(vec_norm2(v))
    return simplify_matrix(v / sp.sqrt(n2))


def run() -> dict[str, object]:
    n = 8
    S = shift(n)
    I8 = sp.eye(n)
    I2 = sp.eye(2)
    J = sp.Matrix([[0, -1], [1, 0]])
    psi0 = sp.Matrix([1, 0])

    G = sp.kronecker_product(S, J)
    I16 = sp.eye(16)

    # Exact finite history state.
    Psi = sp.zeros(16, 1)
    Rt = sp.eye(2)
    for t in range(n):
        ket_t = sp.zeros(n, 1)
        ket_t[t, 0] = 1
        Psi += sp.kronecker_product(ket_t, Rt * psi0)
        Rt *= J
    Psi = simplify_matrix(Psi / sp.sqrt(n))

    Prel = group_average(G, n)
    Pclock8 = group_average(S, n)
    Pclock = sp.kronecker_product(Pclock8, I2)

    seed_clock = sp.zeros(n, 1)
    seed_clock[0, 0] = 1
    seed = sp.kronecker_product(seed_clock, psi0)
    projected_seed = simplify_matrix(Prel * seed)
    projected_seed_norm = sp.sqrt(sp.simplify(vec_norm2(projected_seed)))
    normalized_projected_seed = simplify_matrix(projected_seed / projected_seed_norm)

    checks = {
        "J_squared_is_minus_I": zero(J * J + I2),
        "J_eighth_power_is_I": zero(J ** 8 - I2),
        "combined_constraint_has_order_8": zero(G ** 8 - I16),
        "relational_history_is_normalized": sp.simplify(vec_norm2(Psi) - 1) == 0,
        "relational_history_is_G_invariant": zero(G * Psi - Psi),
        "combined_group_average_is_projector": zero(Prel * Prel - Prel) and zero(Prel.H - Prel),
        "combined_projector_generates_history_from_seed": zero(normalized_projected_seed - Psi),
        "clock_only_group_average_is_projector": zero(Pclock * Pclock - Pclock) and zero(Pclock.H - Pclock),
        "clock_only_average_kills_nontrivial_relational_history": zero(Pclock * Psi),
    }

    # Conditional states: global state is invariant, but clock-conditioned
    # geometry follows R^t exactly.
    conditional_rows = []
    Rt = sp.eye(2)
    conditional_ok = True
    recurrence_ok = True
    previous = None
    for t in range(n):
        raw = conditional_system_state(Psi, t)
        cond = normalized(raw)
        target = simplify_matrix(Rt * psi0)
        ok = zero(cond - target)
        conditional_ok &= ok
        if previous is not None:
            recurrence_ok &= zero(cond - J * previous)
        conditional_rows.append({
            "t": t,
            "conditional_state": [str(x) for x in cond],
            "target_R_power_state": [str(x) for x in target],
            "matches": bool(ok),
        })
        previous = cond
        Rt *= J

    checks["conditioning_recovers_R_power_system_history"] = bool(conditional_ok)
    checks["conditional_states_obey_one_step_R_recurrence"] = bool(recurrence_ok)

    # Character correlation. S characters use zeta^m. J has eigenvalues
    # +i=zeta^2 and -i=zeta^6. G-invariance requires m_clock+r_system=0 mod 8.
    Pm = [projector_character(S, m) for m in range(n)]
    Qplus = simplify_matrix((I2 - sp.I * J) / 2)   # J eigenvalue +i = zeta^2
    Qminus = simplify_matrix((I2 + sp.I * J) / 2)  # J eigenvalue -i = zeta^6

    matched = simplify_matrix(
        sp.kronecker_product(Pm[6], Qplus)
        + sp.kronecker_product(Pm[2], Qminus)
    )
    checks["combined_projector_equals_matched_character_sum"] = zero(Prel - matched)

    character_rows = []
    nonzero_pairs = []
    for m in range(n):
        for label, Q, r in (("J=+i", Qplus, 2), ("J=-i", Qminus, 6)):
            block = sp.kronecker_product(Pm[m], Q)
            weight = sp.simplify(vec_norm2(block * Psi))
            if sp.simplify(weight) != 0:
                nonzero_pairs.append([m, r, str(weight)])
            character_rows.append({
                "clock_m": m,
                "system_sector": label,
                "system_character_r": r,
                "constraint_sum_mod8": (m + r) % 8,
                "weight": str(weight),
            })

    checks["only_constraint_matched_character_pairs_survive"] = (
        nonzero_pairs == [[2, 6, "1/2"], [6, 2, "1/2"]]
    )

    # Clock-only trivial projector selects m=0, whereas the relational history
    # lives entirely in matched nontrivial m=2,6 sectors.
    clock_m0_weight = sp.simplify(vec_norm2(sp.kronecker_product(Pm[0], I2) * Psi))
    checks["relational_history_has_zero_clock_trivial_character_weight"] = clock_m0_weight == 0

    # Seed projection normalization is exactly 1/sqrt(8), as expected for one
    # representative of an eight-element gauge orbit.
    checks["seed_projected_norm_is_one_over_sqrt8"] = sp.simplify(projected_seed_norm - 1 / sp.sqrt(8)) == 0

    passed = bool(all(checks.values()))
    return {
        "status": "exact finite-C8 relational combined-projector positive control",
        "passed": passed,
        "definitions": {
            "clock_shift": "S|t>=|t+1 mod 8>",
            "system_positive_control": "R=J=[[0,-1],[1,0]]",
            "combined_constraint_step": "G=S_clock tensor R_geom",
            "relational_projector": "P_rel=(1/8) sum_tau G^tau",
            "clock_only_projector": "P_clock=(1/8) sum_tau S^tau tensor I",
        },
        "conditional_history": conditional_rows,
        "character_correlation": character_rows,
        "nonzero_character_pairs_clock_m_system_r_weight": nonzero_pairs,
        "checks": checks,
        "projector_resolution": (
            "Clock-only untwisted averaging kills this nontrivial history, but combined group averaging does not. The physical invariant sector is built from clock/system character pairs whose phases cancel in the total constraint."
        ),
        "relational_time_statement": (
            "The global state is exactly invariant under the combined constraint while conditioning on the clock basis recovers a nontrivial sequence psi(t)=R^t psi(0). This is the finite-dimensional Page-Wootters/rigging-map mechanism in exact algebraic form."
        ),
        "next_physical_gate": (
            "Replace the positive-control R=J and externally declared C8 clock factor by objects derived from the actual q=2 graph-changing constraint/history construction. Then test whether the resulting combined projector yields the physical metric-source generating functional needed for Gamma^(2)_metric."
        ),
        "claim_boundary": (
            "Exact relational-projector positive control only. A separate physical clock degree of freedom has not been derived from the microscopic theory, R=J is not claimed to be the physical gravitational evolution operator, and this does not yet define physical time, the physical inner product, the graviton propagator, or g_YC^gravity."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run()
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
