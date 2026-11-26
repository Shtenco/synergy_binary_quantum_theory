#!/usr/bin/env python3
"""Exact SU(2) volume block for the depth-6 vacuum adjoint intermediate sector.

For doubled local spins (1,1,1,0), the nontrivial part of the four-valent
oriented volume is the triple product on three spin-1/2 legs.  The fourth leg
is spin zero, and the source convention contributes the factor (-1)^3.
This script uses exact Pauli matrices and exact Clebsch--Gordan coefficients;
it does not verify the subsequent holonomy recoupling cancellation.
"""
from __future__ import annotations

import itertools
import json

import sympy as sp
from sympy.physics.wigner import clebsch_gordan


def coupled_basis(k_doubled: int) -> sp.Matrix:
    """Return |((1/2,1/2)K,1/2)J=1/2,M=1/2> in computational order."""
    half = sp.Rational(1, 2)
    k = sp.Rational(k_doubled, 2)
    mvals = (half, -half)
    vector = sp.zeros(8, 1)
    for i, m1 in enumerate(mvals):
        for j, m2 in enumerate(mvals):
            for ell, m3 in enumerate(mvals):
                m12 = m1 + m2
                if m12 + m3 != half:
                    continue
                vector[4 * i + 2 * j + ell] = (
                    clebsch_gordan(half, half, k, m1, m2, m12)
                    * clebsch_gordan(k, half, half, m12, m3, half)
                )
    return vector


def main() -> None:
    pauli_generators = (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )
    q123 = sp.zeros(8)
    for permutation in itertools.permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        epsilon = -1 if inversions % 2 else 1
        q123 += epsilon * sp.kronecker_product(
            *(pauli_generators[a] for a in permutation)
        )

    basis = sp.Matrix.hstack(coupled_basis(0), coupled_basis(2))
    basis_gram = sp.simplify(basis.conjugate().T * basis)
    q123_reduced = sp.simplify(basis.conjugate().T * q123 * basis)
    # In the four-leg convention Q=sum_r (-1)^r Q_(all legs except r).
    # Since j_4=0, all terms except r=3 vanish.
    q_source = sp.simplify(-q123_reduced)
    q2 = sp.simplify(q_source * q_source)
    expected_q2 = sp.Rational(3, 16) * sp.eye(2)
    volume_scalar = sp.simplify(sp.sqrt(sp.sqrt(sp.Rational(3, 16))))
    volume_matrix = volume_scalar * sp.eye(2)

    if basis_gram != sp.eye(2) or q2 != expected_q2:
        raise ArithmeticError("exact coupled-basis or volume-block identity failed")

    print(json.dumps({
        "doubled_spins": [1, 1, 1, 0],
        "doubled_total_spin": 1,
        "coupling_labels_doubled": [[0, 1], [2, 1]],
        "basis_gram": sp.sstr(basis_gram),
        "Q123": sp.sstr(q123_reduced),
        "Q_source": sp.sstr(q_source),
        "Q_source_squared": sp.sstr(q2),
        "volume_matrix": sp.sstr(volume_matrix),
        "scope": "exact local volume block only; no proof of the full adjoint index-channel cancellation",
    }, indent=2))


if __name__ == "__main__":
    main()
