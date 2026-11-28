"""Disjunctive branch-aware factorized support for W5 KKV.

Each alternative is a separate factorized marginal envelope. Distinct
oriented specs, T words, adjoints, commutator branches, and conjugation
branches are not coordinate-wise unioned. Within one alternative the
representation is still a Cartesian-product upper bound and does not retain
all cross-node correlations.
"""
from __future__ import annotations

from itertools import product

import bqg_w5_factorized_operator_support as FO
import bqg_w5_su2_structural_basis as SB
import k5_peter_weyl_safe_hda_column as PW


def _signature(support):
    rows = []
    for spins, factors in support.items():
        rows.append((
            tuple(spins),
            tuple(tuple(sorted(labs, key=repr)) for labs in factors),
        ))
    return tuple(sorted(rows, key=repr))


def _count(support):
    return sum(len(labs) for factors in support.values() for labs in factors)


def _guard(branches, limit, stage):
    n = sum(_count(s) for s in branches)
    if limit is not None and n > limit:
        raise RuntimeError(
            f"disjunctive {stage} storage {n} exceeds guard {limit}; no truncation"
        )


def _dedup(branches, limit, stage):
    seen = set()
    out = []
    for branch in branches:
        if not branch:
            continue
        sig = _signature(branch)
        if sig not in seen:
            seen.add(sig)
            out.append(branch)
    _guard(out, limit, stage)
    return out


def _volume(branches, vertex, limit):
    return _dedup(
        [FO._volume(s, vertex, limit) for s in branches],
        limit, "volume alternatives"
    )


def _hit_project(branches, x, y, Jmax2, charged, limit):
    return _dedup(
        [FO._hit_project(s, x, y, Jmax2, charged, limit) for s in branches],
        limit, "hit/projection alternatives"
    )


def factorized_he_branches(branches, source_v, Jmax2,
                           charged_nodes=(0, 1), limit=100000):
    """Keep each oriented spec, T word, and adjoint word as a separate disjunct."""
    out = []
    for _sign, spec in PW.oriented_specs(source_v):
        v, a, b, c = spec
        for _coef, seq in PW.T_sequences(v, a, b, c):
            for word in (seq, PW.adjoint_sequence(seq)):
                for incoming in branches:
                    part = FO._sequence(
                        incoming, word, source_v, Jmax2, charged_nodes, limit
                    )
                    if part:
                        out.append(part)
                        _guard(out, limit, "H_E alternatives")
    return _dedup(out, limit, "H_E deduplication")


def factorized_k_branches(branches, source_v, Jmax2,
                          charged_nodes=(0,), limit=100000):
    """Disjunctive V H_E and H_E V branches for K=[V,H_E]."""
    out = []
    for he_branch in factorized_he_branches(
        branches, source_v, Jmax2, charged_nodes, limit
    ):
        part = FO._volume([he_branch] if isinstance(he_branch, list) else he_branch,
                          source_v, limit)
        out.extend(part)
    vin = _volume(
        [FO._lift(s) for s in branches], source_v, limit
    )
    vin = _dedup(
        [FO._project(s, charged_nodes, limit) for s in vin],
        limit, "H_E V projection alternatives"
    )
    out.extend(factorized_he_branches(
        vin, source_v, Jmax2, charged_nodes, limit
    ))
    return _dedup(out, limit, "K commutator alternatives")


def factorized_c_volume_branches(branches, source_v, target_v, Jmax2,
                                 i, j, limit=100000):
    out = []
    if i == j:
        out.extend(_volume(branches, source_v, limit))
    inverse = _hit_project(
        [FO._lift(s) for s in branches], target_v, source_v, Jmax2,
        (source_v, target_v), limit
    )
    middle = _volume(inverse, source_v, limit)
    out.extend(_hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), limit
    ))
    return _dedup(out, limit, "C(V) alternatives")


def factorized_c_k_branches(branches, source_v, target_v, Jmax2,
                            i, j, limit=100000):
    out = []
    if i == j:
        out.extend(factorized_k_branches(
            branches, source_v, Jmax2, (source_v,), limit
        ))
    inverse = _hit_project(
        [FO._lift(s) for s in branches], target_v, source_v, Jmax2,
        (source_v, target_v), limit
    )
    middle = factorized_k_branches(
        inverse, source_v, Jmax2, (source_v, target_v), limit
    )
    out.extend(_hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), limit
    ))
    return _dedup(out, limit, "C(K) alternatives")


def _scalar_filter(branches, source_v, limit):
    out = []
    for support in branches:
        filtered = {}
        for spins, factors in support.items():
            keep = frozenset(
                lab for lab in factors[source_v]
                if lab[0] == "C" and lab[1] in (0, 2)
            )
            if not keep:
                continue
            updated = list(factors)
            updated[source_v] = keep
            filtered[spins] = tuple(updated)
        if filtered:
            out.append(filtered)
    return _dedup(out, limit, "middle scalar-filter alternatives")


def factorized_triple_branches(branches, source_v, a, b, c, Jmax2,
                               adjoint=False, limit=100000):
    """Forward KKV or independently reversed adjoint as disjunctive branches."""
    out = []
    for i, j, k in product(range(2), repeat=3):
        if not adjoint:
            q = factorized_c_volume_branches(
                branches, source_v, c, Jmax2, k, i, limit
            )
            q = factorized_c_k_branches(
                q, source_v, b, Jmax2, j, k, limit
            )
            q = _scalar_filter(q, source_v, limit)
            q = factorized_c_k_branches(
                q, source_v, a, Jmax2, i, j, limit
            )
        else:
            q = factorized_c_k_branches(
                branches, source_v, a, Jmax2, j, i, limit
            )
            q = factorized_c_k_branches(
                q, source_v, b, Jmax2, k, j, limit
            )
            q = _scalar_filter(q, source_v, limit)
            q = factorized_c_volume_branches(
                q, source_v, c, Jmax2, i, k, limit
            )
        out.extend(q)
        out = _dedup(out, limit, "KKV disjunction")
    return out


def contains_key(branches, key):
    return any(SB.factorized_contains_key(branch, key) for branch in branches)


def profile(branches):
    return {
        "branch_count": len(branches),
        "stored_marginal_labels": sum(_count(s) for s in branches),
        "spin_sector_count_sum": sum(len(s) for s in branches),
        "largest_branch_labels": max((_count(s) for s in branches), default=0),
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "within_branch_cross_node_correlations_preserved": False,
    }
