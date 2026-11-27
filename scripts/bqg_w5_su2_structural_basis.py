#!/usr/bin/env python3
"""Threshold-free SU(2) representation-label basis for fixed-spin W5 sectors.

This module enumerates admissible Peter-Weyl recoupling labels from exact
triangle/parity rules. It does not inspect amplitudes, evaluate CG coefficients,
or claim that every admissible label is reached by the BQG Lorentzian operator.
The generated set is a representation-theoretic ambient basis for a *fixed*
ten-link spin assignment.
"""
from __future__ import annotations

from itertools import product
from math import prod
from typing import Iterator

VERTICES = tuple(range(5))
EDGES = tuple((a, b) for a in VERTICES for b in range(a + 1, 5))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
W5_K = (0, 2, 4, 6, 8)


def allowed_pair2(s1: int, s2: int) -> tuple[int, ...]:
    """Doubled SU(2) spins in s1 tensor s2."""
    if min(s1, s2) < 0:
        raise ValueError("doubled spins must be nonnegative")
    return tuple(range(abs(s1 - s2), s1 + s2 + 1, 2))


def local_spins(spins: tuple[int, ...], vertex: int) -> tuple[int, ...]:
    if len(spins) != 10:
        raise ValueError("a K5 spin key must contain ten link labels")
    if vertex not in VERTICES:
        raise ValueError("vertex must lie in 0..4")
    return tuple(spins[EDGE_INDEX[tuple(sorted((vertex, w)))]]
                 for w in VERTICES if w != vertex)


def gauss_labels(spins_local: tuple[int, int, int, int]) -> tuple[int, ...]:
    """Four-valent singlet labels K12=K34, in doubled-spin convention."""
    a, b, c, d = spins_local
    return tuple(sorted(set(allowed_pair2(a, b)) &
                        set(allowed_pair2(c, d))))


def charged_labels(spins_local: tuple[int, int, int, int],
                   J2: int) -> tuple[tuple[int, int], ...]:
    """(K12,K34) labels coupling to the explicit total representation J2."""
    if J2 < 0:
        raise ValueError("J2 must be nonnegative")
    a, b, c, d = spins_local
    left = allowed_pair2(a, b)
    right = allowed_pair2(c, d)
    return tuple((k12, k34) for k12 in left for k34 in right
                 if abs(k12 - k34) <= J2 <= k12 + k34
                 and (k12 + k34 + J2) % 2 == 0)


def magnetic_labels(J2: int) -> tuple[int, ...]:
    if J2 < 0:
        raise ValueError("J2 must be nonnegative")
    return tuple(range(-J2, J2 + 1, 2))


def all_total_J2(spins_local: tuple[int, int, int, int]) -> tuple[int, ...]:
    upper = sum(spins_local)
    return tuple(J2 for J2 in range(upper + 1)
                 if charged_labels(spins_local, J2))


def gauss_node_basis(spins_local: tuple[int, int, int, int]):
    return tuple(("G", k) for k in gauss_labels(spins_local))


def charged_node_basis(spins_local: tuple[int, int, int, int]):
    return tuple(("C", J2, M2, k12, k34)
                 for J2 in all_total_J2(spins_local)
                 for M2 in magnetic_labels(J2)
                 for k12, k34 in charged_labels(spins_local, J2))


def covariant_node_factors(spins: tuple[int, ...], source_v: int):
    """Return four Gauss factors and one complete charged source factor."""
    if len(spins) != 10 or any(s < 0 for s in spins):
        raise ValueError("spins must be ten nonnegative doubled link labels")
    if source_v not in VERTICES:
        raise ValueError("source_v must lie in 0..4")
    factors = []
    for v in VERTICES:
        ls = local_spins(spins, v)
        factors.append(charged_node_basis(ls) if v == source_v
                       else gauss_node_basis(ls))
    return tuple(factors)


def covariant_basis_size(spins: tuple[int, ...], source_v: int) -> int:
    return prod(len(factor) for factor in covariant_node_factors(spins, source_v))


def iter_covariant_keys(spins: tuple[int, ...], source_v: int
                        ) -> Iterator[tuple[tuple[int, ...], tuple[tuple, ...]]]:
    """Lazily enumerate the complete factorized covariant label basis."""
    factors = covariant_node_factors(spins, source_v)
    for labels in product(*factors):
        yield tuple(spins), tuple(labels)


def w5_gauss_seed(K: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if K not in W5_K:
        raise ValueError(f"W5 K must be one of {W5_K}")
    return (4,) * 10, (K, 0, 0, 0, 0)


def w5_seed_as_covariant(K: int, source_v: int = 0):
    """Embed a W5 Gauss seed into the J=0 charged source convention."""
    spins, Ks = w5_gauss_seed(K)
    labels = tuple(("C", 0, 0, Ks[v], Ks[v]) if v == source_v
                   else ("G", Ks[v]) for v in VERTICES)
    return spins, labels


def advance_spin_path(spins: tuple[int, ...], path: tuple[int, ...],
                      Jmax2: int) -> tuple[tuple[int, ...], ...]:
    """Conservative link-spin support of one fundamental hit along a path.

    Each traversed edge undergoes s -> s-1 or s+1 in doubled-spin units.
    Magnetic-index and intertwiner zeros are intentionally not used, so the
    result is an upper support bound, not a claim of nonzero operator action.
    """
    if Jmax2 < 0:
        raise ValueError("Jmax2 must be nonnegative")
    states = {tuple(spins)}
    for x, y in zip(path[:-1], path[1:]):
        edge_index = EDGE_INDEX[tuple(sorted((x, y)))]
        next_states = set()
        for state in states:
            s = state[edge_index]
            for so in (s - 1, s + 1):
                if 0 <= so <= Jmax2:
                    updated = list(state)
                    updated[edge_index] = so
                    next_states.add(tuple(updated))
        states = next_states
    return tuple(sorted(states))


def spin_support_for_sequence(spins: tuple[int, ...], sequence,
                              Jmax2: int) -> tuple[tuple[int, ...], ...]:
    """Propagate only exact edge-spin selection rules through a PW word."""
    states = {tuple(spins)}
    for op in sequence:
        if op[0] == "V":
            continue
        if op[0] != "P":
            raise ValueError(f"unknown Peter-Weyl sequence opcode: {op[0]}")
        path = tuple(op[1])
        updated = set()
        for state in states:
            updated.update(advance_spin_path(state, path, Jmax2))
        states = updated
    return tuple(sorted(states))


def t_sequence_spin_support(spins: tuple[int, ...], v: int,
                            a: int, b: int, c: int, Jmax2: int,
                            include_adjoint: bool = True):
    """Spin-sector superset for the repository's actual T_sequences words.

    The word templates and adjoint reversal are imported from the existing
    Peter-Weyl engine to avoid maintaining a second copy of their definition.
    """
    import sys
    from pathlib import Path
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    import k5_peter_weyl_safe_hda_column as PW

    support = set()
    for _coefficient, sequence in PW.T_sequences(v, a, b, c):
        support.update(spin_support_for_sequence(spins, sequence, Jmax2))
        if include_adjoint:
            reverse = PW.adjoint_sequence(sequence)
            support.update(spin_support_for_sequence(spins, reverse, Jmax2))
    return tuple(sorted(support))


def validate_covariant_key(key, source_v: int) -> bool:
    spins, labels = key
    if len(labels) != 5:
        return False
    factors = covariant_node_factors(tuple(spins), source_v)
    return all(labels[v] in factors[v] for v in VERTICES)
