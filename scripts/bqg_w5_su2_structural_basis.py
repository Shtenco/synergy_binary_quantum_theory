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



def spin_sector_basis_profile(spin_sectors, source_v: int):
    """Return exact ambient covariant-basis sizes for a finite spin-sector set."""
    unique = sorted({tuple(spins) for spins in spin_sectors})
    return tuple((spins, covariant_basis_size(spins, source_v))
                 for spins in unique)


def charged_nodes_factors(spins: tuple[int, ...], charged_nodes=(0,)):
    """Ambient factors with an arbitrary fixed set of charged nodes."""
    if len(spins) != 10 or any(s < 0 for s in spins):
        raise ValueError("spins must be ten nonnegative doubled link labels")
    charged = tuple(sorted(set(charged_nodes)))
    if any(v not in VERTICES for v in charged):
        raise ValueError("charged nodes must lie in 0..4")
    return tuple(charged_node_basis(local_spins(spins, v))
                 if v in charged else gauss_node_basis(local_spins(spins, v))
                 for v in VERTICES)


def charged_nodes_basis_size(spins: tuple[int, ...], charged_nodes=(0,)) -> int:
    return prod(len(factor) for factor in charged_nodes_factors(spins, charged_nodes))


def iter_charged_nodes_keys(spins: tuple[int, ...], charged_nodes=(0,)):
    factors = charged_nodes_factors(spins, charged_nodes)
    for labels in product(*factors):
        yield tuple(spins), tuple(labels)


def charged_nodes_basis_profile(spin_sectors, charged_nodes=(0,)):
    unique = sorted({tuple(spins) for spins in spin_sectors})
    return tuple((spins, charged_nodes_basis_size(spins, charged_nodes))
                 for spins in unique)


def volume_label_targets(spins_local: tuple[int, int, int, int], label):
    """Selection-rule superset for the scalar canonical volume block.

    The repository applies one reduced block independently in each (J2,M2)
    sector. Hence only multiplicity labels (K12,K34) may mix.
    """
    if len(label) != 5 or label[0] != "C":
        raise ValueError("volume requires a charged-node label")
    _, J2, M2, K12, K34 = label
    if M2 not in magnetic_labels(J2):
        raise ValueError("invalid magnetic label for J2")
    rec = charged_labels(tuple(spins_local), J2)
    if (K12, K34) not in rec:
        raise ValueError("input charged label is not admissible")
    return tuple(("C", J2, M2, a, b) for a, b in rec)


def gauss_label_to_charged(label):
    if len(label) != 2 or label[0] != "G":
        raise ValueError("expected a Gauss label")
    return ("C", 0, 0, label[1], label[1])


def lift_key_to_all_charged(key):
    """Lift projected endpoint labels to the unprojected all-node J basis."""
    spins, labels = key
    lifted = tuple(gauss_label_to_charged(lab) if lab[0] == "G" else lab
                   for lab in labels)
    return tuple(spins), lifted


def project_all_charged_key(key, charged_nodes=(0,)):
    """Project final labels to charged sectors at selected nodes and Gauss at others."""
    spins, labels = key
    charged = set(charged_nodes)
    out = []
    for v, lab in enumerate(labels):
        if v in charged:
            if lab[0] != "C":
                raise ValueError("charged endpoint lost its total-J label")
            out.append(lab)
        else:
            if lab[0] != "C" or lab[1] != 0 or lab[2] != 0 or lab[3] != lab[4]:
                return None
            out.append(("G", lab[3]))
    return tuple(spins), tuple(out)


def advance_hit_label_support(key, x: int, y: int, Jmax2: int):
    """Conservative all-node label support for one fundamental edge hit.

    Only the two endpoint link spins and their local representation labels
    change. Endpoint labels are allowed to occupy every admissible total-J
    sector until the final Gauss/charged projection.
    """
    spins, labels = lift_key_to_all_charged(key)
    edge = tuple(sorted((x, y)))
    if edge not in EDGE_INDEX:
        raise ValueError("hit endpoints must define a K5 edge")
    ei = EDGE_INDEX[edge]
    s = spins[ei]
    endpoints = (x, y)
    for so in (s - 1, s + 1):
        if so < 0 or so > Jmax2:
            continue
        new_spins = list(spins)
        new_spins[ei] = so
        new_spins = tuple(new_spins)
        opts = []
        for v in endpoints:
            opts.append(charged_node_basis(local_spins(new_spins, v)))
        for left, right in product(*opts):
            new_labels = list(labels)
            new_labels[x] = left
            new_labels[y] = right
            yield new_spins, tuple(new_labels)


def advance_volume_label_support(key, vertex: int):
    """Apply scalar-volume selection rules at one already charged node."""
    spins, labels = key
    lab = labels[vertex]
    for target in volume_label_targets(local_spins(spins, vertex), lab):
        updated = list(labels)
        updated[vertex] = target
        yield spins, tuple(updated)


def structural_sequence_support(key, sequence, source_v: int, Jmax2: int,
                                charged_nodes=(0,), max_states=None):
    """Conservative label support through one actual Peter-Weyl sequence.

    Intermediate nodes are kept in complete total-J sectors; the requested
    Gauss/charged projection is imposed only after the whole word.
    """
    states = {lift_key_to_all_charged(key)}
    for op in sequence:
        next_states = set()
        if op[0] == "V":
            for state in states:
                next_states.update(advance_volume_label_support(state, source_v))
        elif op[0] == "P":
            path = tuple(op[1])
            for state in states:
                partial = {state}
                for x, y in zip(path[:-1], path[1:]):
                    partial = {nxt for st in partial
                               for nxt in advance_hit_label_support(st, x, y, Jmax2)}
                    if max_states is not None and len(partial) > max_states:
                        raise RuntimeError("structural support exceeded max_states; no truncation performed")
                next_states.update(partial)
                if max_states is not None and len(next_states) > max_states:
                    raise RuntimeError("structural support exceeded max_states; no truncation performed")
        else:
            raise ValueError(f"unknown Peter-Weyl sequence opcode: {op[0]}")
        states = next_states
        if max_states is not None and len(states) > max_states:
            raise RuntimeError("structural support exceeded max_states; no truncation performed")
        if not states:
            break
    projected = set()
    for state in states:
        item = project_all_charged_key(state, charged_nodes)
        if item is not None:
            projected.add(item)
    return projected


def validate_covariant_key(key, source_v: int) -> bool:
    spins, labels = key
    if len(labels) != 5:
        return False
    factors = covariant_node_factors(tuple(spins), source_v)
    return all(labels[v] in factors[v] for v in VERTICES)
