"""Factorized operator-level support envelopes for the frozen W5 KKV word.

The representation stores per-spin-sector marginal label sets at each node.
It is deliberately a Cartesian-product *upper bound*: correlations between
nodes are not represented. No amplitudes or physical non-vanishing claims.
"""
from __future__ import annotations

from itertools import product

import bqg_w5_su2_structural_basis as SB
import k5_peter_weyl_safe_hda_column as PW


def _count(support):
    return sum(len(x) for factors in support.values() for x in factors)


def _guard(support, limit, stage):
    n = _count(support)
    if limit is not None and n > limit:
        raise RuntimeError(
            f"factorized {stage} labels {n} exceed guard {limit}; no truncation"
        )


def _merge(dst, spins, factors):
    factors = tuple(frozenset(x) for x in factors)
    if spins not in dst:
        dst[spins] = factors
    else:
        dst[spins] = tuple(
            frozenset(a | b) for a, b in zip(dst[spins], factors)
        )


def _union(*supports, limit=None, stage="union"):
    out = {}
    for support in supports:
        for spins, factors in support.items():
            _merge(out, spins, factors)
            _guard(out, limit, stage)
    return out


def _lift(support):
    """Replace projected Gauss labels by their unique covariant J=M=0 label."""
    out = {}
    for spins, factors in support.items():
        lifted = []
        for v, labs in enumerate(factors):
            lifted.append(frozenset(
                ("C", 0, 0, lab[1], lab[1]) if lab[0] == "G" else lab
                for lab in labs
            ))
        _merge(out, spins, lifted)
    return out


def _volume(support, vertex, limit):
    out = {}
    for spins, factors in support.items():
        updated = list(factors)
        targets = set()
        for lab in factors[vertex]:
            targets.update(SB.volume_label_targets(
                SB.local_spins(spins, vertex), lab
            ))
        updated[vertex] = frozenset(targets)
        _merge(out, spins, updated)
    _guard(out, limit, "volume")
    return out


def _hit(support, x, y, Jmax2, limit):
    edge = tuple(sorted((x, y)))
    if x == y or edge not in SB.EDGE_INDEX:
        raise ValueError("holonomy hit must connect two distinct K5 nodes")
    ei = SB.EDGE_INDEX[edge]
    out = {}
    for spins, factors in support.items():
        for so in (spins[ei] - 1, spins[ei] + 1):
            if so < 0 or so > Jmax2:
                continue
            new_spins = list(spins)
            new_spins[ei] = so
            new_spins = tuple(new_spins)
            updated = list(factors)
            valid = True
            for v in (x, y):
                all_labels = SB.charged_node_basis(
                    SB.local_spins(new_spins, v)
                )
                possible = set()
                for old in factors[v]:
                    old_j2 = old[1] if old[0] == "C" else 0
                    possible.update(
                        lab for lab in all_labels
                        if abs(lab[1] - old_j2) <= 1
                        and lab[1] + old_j2 >= 1
                    )
                if not possible:
                    valid = False
                    break
                updated[v] = frozenset(possible)
            if valid:
                _merge(out, new_spins, updated)
                _guard(out, limit, "holonomy")
    return out


def _project(support, charged_nodes, limit):
    charged = set(charged_nodes)
    out = {}
    for spins, factors in support.items():
        projected = []
        valid = True
        for v, labs in enumerate(factors):
            if v in charged:
                keep = frozenset(labs)
            else:
                keep = frozenset(
                    ("G", k)
                    for k in SB.gauss_labels(SB.local_spins(spins, v))
                    if ("C", 0, 0, k, k) in labs
                )
                if not keep:
                    valid = False
                    break
            projected.append(keep)
        if valid:
            _merge(out, spins, projected)
            _guard(out, limit, "Gauss projection")
    return out


def _sequence(support, sequence, source_v, Jmax2, charged_nodes, limit):
    work = _lift(support)
    for op in sequence:
        if op[0] == "V":
            work = _volume(work, source_v, limit)
        elif op[0] == "P":
            path = tuple(op[1])
            for x, y in zip(path[:-1], path[1:]):
                work = _hit(work, x, y, Jmax2, limit)
        else:
            raise ValueError(f"unknown Peter-Weyl opcode: {op[0]}")
        if not work:
            break
    return _project(work, charged_nodes, limit)


def _t_sequence(support, v, a, b, c, source_v, Jmax2,
                charged_nodes, limit, include_adjoint=True):
    out = {}
    for _coef, seq in PW.T_sequences(v, a, b, c):
        words = [seq]
        if include_adjoint:
            words.append(PW.adjoint_sequence(seq))
        for word in words:
            part = _sequence(
                support, word, source_v, Jmax2, charged_nodes, limit
            )
            out = _union(out, part, limit=limit, stage="T sequence union")
    return out


def factorized_he(support, source_v, Jmax2, charged_nodes=(0, 1),
                  max_factor_labels=100000):
    """Marginal envelope of the complete oriented H_E primitive."""
    out = {}
    for _sign, spec in PW.oriented_specs(source_v):
        v, a, b, c = spec
        part = _t_sequence(
            support, v, a, b, c, source_v, Jmax2, charged_nodes,
            max_factor_labels, True
        )
        out = _union(out, part, limit=max_factor_labels, stage="H_E")
    return out


def factorized_k(support, source_v, Jmax2, charged_nodes=(0,),
                 max_factor_labels=100000):
    """Envelope of K=[V,H_E], unioning V-after-H_E and H_E-after-V."""
    he = factorized_he(
        support, source_v, Jmax2, charged_nodes, max_factor_labels
    )
    vh = _volume(he, source_v, max_factor_labels)
    vin = _volume(_lift(support), source_v, max_factor_labels)
    vin = _project(vin, charged_nodes, max_factor_labels)
    hv = factorized_he(
        vin, source_v, Jmax2, charged_nodes, max_factor_labels
    )
    return _union(vh, hv, limit=max_factor_labels, stage="K commutator")


def _hit_project(support, x, y, Jmax2, charged_nodes, limit):
    return _project(
        _hit(_lift(support), x, y, Jmax2, limit), charged_nodes, limit
    )


def factorized_c_volume(support, source_v, target_v, Jmax2, i, j,
                        max_factor_labels=100000):
    """Envelope of C_e(V)=V-h V h^-1, including inverse/closure hits."""
    out = {}
    if i == j:
        out = _volume(support, source_v, max_factor_labels)
    inverse = _hit_project(
        support, target_v, source_v, Jmax2, (source_v, target_v),
        max_factor_labels
    )
    middle = _volume(inverse, source_v, max_factor_labels)
    closed = _hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), max_factor_labels
    )
    return _union(out, closed, limit=max_factor_labels, stage="C(V)")


def factorized_c_k(support, source_v, target_v, Jmax2, i, j,
                   max_factor_labels=100000):
    """Envelope of C_e(K)=K-h K h^-1."""
    out = {}
    if i == j:
        out = factorized_k(
            support, source_v, Jmax2, (source_v,), max_factor_labels
        )
    inverse = _hit_project(
        support, target_v, source_v, Jmax2, (source_v, target_v),
        max_factor_labels
    )
    middle = factorized_k(
        inverse, source_v, Jmax2, (source_v, target_v), max_factor_labels
    )
    closed = _hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), max_factor_labels
    )
    return _union(out, closed, limit=max_factor_labels, stage="C(K)")


def _scalar_filter(support, source_v):
    """Apply the frozen source-node J2 in {0,2} filter to each marginal."""
    out = {}
    for spins, factors in support.items():
        keep = frozenset(
            lab for lab in factors[source_v]
            if lab[0] == "C" and lab[1] in (0, 2)
        )
        if not keep:
            continue
        updated = list(factors)
        updated[source_v] = keep
        _merge(out, spins, updated)
    return out


def factorized_triple(support, source_v, a, b, c, Jmax2,
                      adjoint=False, max_factor_labels=100000):
    """Forward or independently reversed-adjoint factorized KKV envelope."""
    out = {}
    for i, j, k in product(range(2), repeat=3):
        if not adjoint:
            q = factorized_c_volume(
                support, source_v, c, Jmax2, k, i, max_factor_labels
            )
            q = factorized_c_k(
                q, source_v, b, Jmax2, j, k, max_factor_labels
            )
            q = _scalar_filter(q, source_v)
            q = factorized_c_k(
                q, source_v, a, Jmax2, i, j, max_factor_labels
            )
        else:
            q = factorized_c_k(
                support, source_v, a, Jmax2, j, i, max_factor_labels
            )
            q = factorized_c_k(
                q, source_v, b, Jmax2, k, j, max_factor_labels
            )
            q = _scalar_filter(q, source_v)
            q = factorized_c_volume(
                q, source_v, c, Jmax2, i, k, max_factor_labels
            )
        out = _union(
            out, q, limit=max_factor_labels, stage="KKV triple"
        )
    return out


def factorized_support_profile(support):
    return SB.factorized_support_profile(support)
