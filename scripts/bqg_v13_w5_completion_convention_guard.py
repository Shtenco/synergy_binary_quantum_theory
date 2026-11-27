#!/usr/bin/env python3
"""Deterministic operator-valued matrix guard for the v1.3/W5 completion split.

This is an algebraic convention guard, not a physical-operator selector.  It
checks the v1.3 direct-Hermitian 12-word identity against the independently
formed 24-term raw word and displays its difference from the W5 plus-completion
after K_comm -> K_H conversion.  It does not use Peter-Weyl production data.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def ptr_aux(x: np.ndarray, n_aux: int, n_graph: int) -> np.ndarray:
    blocks = x.reshape(n_aux, n_graph, n_aux, n_graph)
    return sum((blocks[i, :, i, :] for i in range(n_aux)),
               np.zeros((n_graph, n_graph), dtype=np.complex128))


def permutation_sign(p: tuple[int, int, int]) -> int:
    inversions = sum(
        p[i] > p[j]
        for i in range(3)
        for j in range(i + 1, 3)
    )
    return -1 if inversions % 2 else 1


def run(seed: int = 20261128, n_aux: int = 2, n_graph: int = 5) -> dict:
    rng = np.random.default_rng(seed)
    dim = n_aux * n_graph
    k_legs = []
    v_legs = []
    for _ in range(4):
        x = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        k_legs.append(x - x.conj().T)
        y = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        v_legs.append(0.5 * (y + y.conj().T))

    # Real orientation coefficients are deliberately explicit and frozen.
    eta = (1.0, -1.0, 1.0, -1.0)
    raw = np.zeros((n_graph, n_graph), dtype=np.complex128)
    direct12 = np.zeros_like(raw)

    for omitted in range(4):
        face = tuple(x for x in range(4) if x != omitted)
        for p in itertools.permutations(face):
            raw += eta[omitted] * permutation_sign(p) * ptr_aux(
                k_legs[p[0]] @ k_legs[p[1]] @ v_legs[p[2]],
                n_aux, n_graph,
            )

        cyclic = (
            (face[0], face[1], face[2]),
            (face[1], face[2], face[0]),
            (face[2], face[0], face[1]),
        )
        for a, b, c in cyclic:
            comm = k_legs[a] @ k_legs[b] - k_legs[b] @ k_legs[a]
            anti = comm @ v_legs[c] + v_legs[c] @ comm
            direct12 += -0.5j * eta[omitted] * ptr_aux(
                anti, n_aux, n_graph
            )

    raw_dagger = raw.conj().T

    # v1.3: K_comm is anti-Hermitian; phase completion is -i/2(B-B†).
    s_v13 = -0.5j * (raw - raw_dagger)

    # W5 ledger: convert each of the two K_comm legs to K_H=(i hbar)^-1 K_comm,
    # giving B_H=-hbar^-2 B_comm, then take the symmetric part.
    hbar = 1.0
    b_h = -raw / hbar**2
    s_w5_plus = 0.5 * (b_h + b_h.conj().T)
    # If the frozen five-bracket/Lorentzian prefactor contributes an explicit
    # i, the Hermitian projection of i*B_H reproduces the v1.3 phase choice.
    s_w5_with_i_phase = 0.5 * (1j * b_h + (1j * b_h).conj().T)

    scale = max(float(np.linalg.norm(s_v13)), 1e-30)
    direct_rel = float(np.linalg.norm(direct12 - s_v13) / scale)
    v13_herm = float(
        np.linalg.norm(s_v13 - s_v13.conj().T)
        / max(float(np.linalg.norm(s_v13)), 1e-30)
    )
    w5_herm = float(
        np.linalg.norm(s_w5_plus - s_w5_plus.conj().T)
        / max(float(np.linalg.norm(s_w5_plus)), 1e-30)
    )
    split_rel = float(
        np.linalg.norm(s_v13 - s_w5_plus)
        / max(float(np.linalg.norm(s_v13)), float(np.linalg.norm(s_w5_plus)), 1e-30)
    )

    checks = {
        "all_K_legs_antihermitian": all(
            np.linalg.norm(x + x.conj().T) < 1e-11 for x in k_legs
        ),
        "all_V_legs_hermitian": all(
            np.linalg.norm(x - x.conj().T) < 1e-11 for x in v_legs
        ),
        "direct_12_word_equals_v13_24_plus_24": direct_rel < 1e-12,
        "v13_completion_is_hermitian": v13_herm < 1e-12,
        "w5_plus_completion_is_hermitian": w5_herm < 1e-12,
        "explicit_i_phase_projection_matches_v13": (
            np.linalg.norm(s_w5_with_i_phase - s_v13) / scale < 1e-12
        ),
        "two_prescriptions_are_distinct_on_witness": split_rel > 1e-6,
    }
    return {
        "status": "PASS_CONVENTION_SPLIT_GUARD" if all(checks.values())
                  else "FAIL_CONVENTION_SPLIT_GUARD",
        "checks": checks,
        "seed": seed,
        "aux_dimension": n_aux,
        "graph_dimension": n_graph,
        "hbar": hbar,
        "eta": list(eta),
        "raw_operator_frobenius_norm": float(np.linalg.norm(raw)),
        "v13_minus_completion_frobenius_norm": float(np.linalg.norm(s_v13)),
        "w5_plus_completion_frobenius_norm": float(np.linalg.norm(s_w5_plus)),
        "w5_i_phase_completion_frobenius_norm": float(np.linalg.norm(s_w5_with_i_phase)),
        "direct12_vs_v13_relative_error": direct_rel,
        "v13_hermiticity_defect": v13_herm,
        "w5_plus_hermiticity_defect": w5_herm,
        "plus_vs_minus_relative_difference": split_rel,
        "v13_definition": "-i/2 * (B_comm - B_comm^dagger)",
        "w5_ledger_definition": "Herm(B_H) with B_H=-hbar^-2*B_comm",
        "phase_restored_definition": "Herm(i*B_H)=i/2*(B_H-B_H^dagger)",
        "interpretation": (
            "The exact 12-word identity reproduces the v1.3 minus-completion. "
            "The W5 plus-completion is also Hermitian but is a different "
            "operator on this generic noncommuting witness. This guard does "
            "not choose the physical convention and does not test W5 "
            "Peter-Weyl component matrices or continuum physics."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261128)
    parser.add_argument("--n-aux", type=int, default=2)
    parser.add_argument("--n-graph", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.seed, args.n_aux, args.n_graph)
    payload = json.dumps(result, indent=2)
    print(payload)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    return 0 if result["status"] == "PASS_CONVENTION_SPLIT_GUARD" else 1


if __name__ == "__main__":
    raise SystemExit(main())
