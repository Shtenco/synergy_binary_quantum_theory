"""Branch-aware pairwise factor-graph support for W5 KKV.

Each alternative is a pairwise CSP over five nodes. Unlike unary
factorization, it preserves every pair relation within a branch. Unlike the
global pairwise baseline, it does not union distinct oriented-spec branches
before later operator composition. It remains an over-approximate abstraction:
three-/four-/five-node correlations and amplitudes are not represented.
"""
from __future__ import annotations

from itertools import product

import bqg_w5_pairwise_factor_graph as PG
import bqg_w5_su2_structural_basis as SB
import k5_peter_weyl_safe_hda_column as PW


def _signature(support):
    rows = []
    for spins, state in support.items():
        unary = tuple(tuple(sorted(x, key=repr)) for x in state["unary"])
        pairs = tuple(
            (p, tuple(sorted(rel, key=repr)))
            for p, rel in sorted(state["pairs"].items())
        )
        rows.append((tuple(spins), unary, pairs))
    return tuple(sorted(rows, key=repr))


def _size(support):
    return sum(
        sum(map(len, state["unary"])) +
        sum(map(len, state["pairs"].values()))
        for state in support.values()
    )


def _guard(branches, limit, stage):
    n = sum(_size(s) for s in branches)
    if limit is not None and n > limit:
        raise RuntimeError(
            f"disjunctive-pairwise {stage} entries {n} exceed guard {limit}; no truncation"
        )


def _dedup(branches, limit, stage):
    seen, out = set(), []
    for branch in branches:
        if not branch:
            continue
        sig = _signature(branch)
        if sig not in seen:
            seen.add(sig)
            out.append(branch)
    _guard(out, limit, stage)
    return out


def factor_graph_he_spec_branches(branches, source_v, Jmax2,
                                  charged_nodes=(0, 1), limit=100000):
    """Keep each oriented spec as a disjunct; union only its T/adjoint words."""
    out = []
    for incoming in branches:
        for _sign, spec in PW.oriented_specs(source_v):
            v, a, b, c = spec
            spec_union = {}
            for _coef, seq in PW.T_sequences(v, a, b, c):
                for word in (seq, PW.adjoint_sequence(seq)):
                    part = PG._sequence(
                        incoming, word, source_v, Jmax2, charged_nodes, limit
                    )
                    spec_union = PG._union(
                        spec_union, part, limit=limit,
                        stage="oriented-spec T-word union"
                    )
            if spec_union:
                out.append(spec_union)
                _guard(out, limit, "oriented-spec alternatives")
    return _dedup(out, limit, "oriented-spec deduplication")


def factor_graph_k_branches(branches, source_v, Jmax2,
                            charged_nodes=(0,), limit=100000):
    out = []
    for he in factor_graph_he_spec_branches(
        branches, source_v, Jmax2, charged_nodes, limit
    ):
        vh = PG._volume(he, source_v, limit)
        if vh:
            out.append(vh)
    vin = []
    for branch in branches:
        lifted = PG._lift(branch)
        vol = PG._volume(lifted, source_v, limit)
        projected = PG._project(vol, charged_nodes, limit)
        if projected:
            vin.append(projected)
    vin = _dedup(vin, limit, "V-before-H_E alternatives")
    out.extend(factor_graph_he_spec_branches(
        vin, source_v, Jmax2, charged_nodes, limit
    ))
    return _dedup(out, limit, "K commutator alternatives")


def _hit_project(branches, x, y, Jmax2, charged, limit):
    out = []
    for branch in branches:
        part = PG._hit_project(branch, x, y, Jmax2, charged, limit)
        if part:
            out.append(part)
    return _dedup(out, limit, "hit/project alternatives")


def factor_graph_c_volume_branches(branches, source_v, target_v, Jmax2,
                                   i, j, limit=100000):
    out = []
    if i == j:
        out.extend(PG._volume(s, source_v, limit) for s in branches)
    inverse = _hit_project(
        [PG._lift(s) for s in branches], target_v, source_v, Jmax2,
        (source_v, target_v), limit
    )
    middle = [PG._volume(s, source_v, limit) for s in inverse]
    out.extend(_hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), limit
    ))
    return _dedup(out, limit, "C(V) alternatives")


def factor_graph_c_k_branches(branches, source_v, target_v, Jmax2,
                              i, j, limit=100000):
    out = []
    if i == j:
        out.extend(factor_graph_k_branches(
            branches, source_v, Jmax2, (source_v,), limit
        ))
    inverse = _hit_project(
        [PG._lift(s) for s in branches], target_v, source_v, Jmax2,
        (source_v, target_v), limit
    )
    middle = factor_graph_k_branches(
        inverse, source_v, Jmax2, (source_v, target_v), limit
    )
    out.extend(_hit_project(
        middle, source_v, target_v, Jmax2, (source_v,), limit
    ))
    return _dedup(out, limit, "C(K) alternatives")


def _scalar_filter(branches, source_v, limit):
    out = []
    for branch in branches:
        filtered = PG._scalar_filter(branch, source_v)
        if filtered:
            out.append(filtered)
    return _dedup(out, limit, "middle scalar-filter alternatives")


def factor_graph_triple_branches(branches, source_v, a, b, c, Jmax2,
                                 adjoint=False, limit=100000):
    """Forward KKV or independently reversed adjoint, preserving disjuncts."""
    out = []
    for i, j, k in product(range(2), repeat=3):
        if not adjoint:
            q = factor_graph_c_volume_branches(
                branches, source_v, c, Jmax2, k, i, limit
            )
            q = factor_graph_c_k_branches(
                q, source_v, b, Jmax2, j, k, limit
            )
            q = _scalar_filter(q, source_v, limit)
            q = factor_graph_c_k_branches(
                q, source_v, a, Jmax2, i, j, limit
            )
        else:
            q = factor_graph_c_k_branches(
                branches, source_v, a, Jmax2, j, i, limit
            )
            q = factor_graph_c_k_branches(
                q, source_v, b, Jmax2, k, j, limit
            )
            q = _scalar_filter(q, source_v, limit)
            q = factor_graph_c_volume_branches(
                q, source_v, c, Jmax2, i, k, limit
            )
        out.extend(q)
        out = _dedup(out, limit, "KKV disjunction")
    return out


def contains_key(branches, key):
    return any(PG.contains_joint(branch, key) for branch in branches)


def profile(branches):
    return {
        "branch_count": len(branches),
        "stored_pairwise_entries": sum(_size(s) for s in branches),
        "spin_sector_count_sum": sum(len(s) for s in branches),
        "largest_branch_entries": max((_size(s) for s in branches), default=0),
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "within_branch_pair_correlations_preserved": True,
        "higher_order_correlations_preserved": False,
    }
