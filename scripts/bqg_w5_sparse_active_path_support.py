"""Sparse active-path pairwise support for W5 primitive words.

Missing pair factors are interpreted as Cartesian products of their unary
domains. This is a relaxation: it may add spurious combinations but cannot
remove a genuinely joint-supported state. Only pairs whose endpoints occur
in the current primitive word's active node scope are stored.
"""
from __future__ import annotations

from itertools import combinations

import bqg_w5_pairwise_factor_graph as PG
import bqg_w5_su2_structural_basis as SB

NODES = tuple(range(5))


def active_scope(sequence, source_v):
    """Nodes touched by any P path, together with the volume/source vertex."""
    active = {source_v}
    for op in sequence:
        if op[0] == "P":
            active.update(tuple(op[1]))
        elif op[0] != "V":
            raise ValueError(f"unknown operator {op[0]}")
    return tuple(sorted(active))


def _stored_pairs(scope):
    return tuple(combinations(scope, 2))


def compress(support, scope):
    """Drop pair factors outside the active induced subgraph."""
    keep = set(_stored_pairs(scope))
    out = {}
    for spins, state in support.items():
        out[spins] = {
            "unary": tuple(frozenset(x) for x in state["unary"]),
            "pairs": {p: frozenset(rel) for p, rel in state["pairs"].items()
                      if p in keep},
        }
    return out


def expand(support):
    """Materialize omitted factors as unary Cartesian products, never constraints."""
    out = {}
    for spins, state in support.items():
        pairs = {}
        for a, b in combinations(NODES, 2):
            p = (a, b)
            pairs[p] = state["pairs"].get(
                p, frozenset((x, y) for x in state["unary"][a]
                             for y in state["unary"][b])
            )
        out[spins] = {"unary": state["unary"], "pairs": pairs}
    return out


def _guard(support, limit, stage):
    size = sum(
        sum(map(len, s["unary"])) + sum(map(len, s["pairs"].values()))
        for s in support.values()
    )
    if limit is not None and size > limit:
        raise RuntimeError(
            f"sparse-active-path {stage} entries {size} exceed guard {limit}; "
            "no truncation"
        )


def sequence(support, word, source_v, Jmax2, charged_nodes=(0,),
             limit=100000):
    """Propagate one primitive word with active-scope pair storage."""
    scope = {source_v}
    work = compress(PG._lift(support), tuple(sorted(scope)))
    _guard(work, limit, "initial active scope")
    for op in word:
        if op[0] == "V":
            scope.add(source_v)
        elif op[0] == "P":
            scope.update(tuple(op[1]))
        else:
            raise ValueError(f"unknown operator {op[0]}")
        full = expand(work)
        if op[0] == "V":
            full = PG._volume(full, source_v, limit)
        else:
            path = tuple(op[1])
            for x, y in zip(path[:-1], path[1:]):
                full = PG._hit(full, x, y, Jmax2, limit)
        work = compress(full, tuple(sorted(scope)))
        _guard(work, limit, f"after {op[0]}")
        if not work:
            break
    return PG._project(expand(work), charged_nodes, limit)


def contains_joint(support, key):
    spins, labels = key
    if spins not in support:
        return False
    state = support[spins]
    if any(labels[v] not in state["unary"][v] for v in NODES):
        return False
    return all((labels[a], labels[b]) in rel
               for (a, b), rel in state["pairs"].items())


def profile(support):
    unary = sum(len(x) for s in support.values() for x in s["unary"])
    pairs = sum(len(x) for s in support.values() for x in s["pairs"].values())
    scope_sizes = sorted({len(s["pairs"]) for s in support.values()})
    return {
        "spin_sector_count": len(support),
        "unary_entries": unary,
        "stored_pair_entries": pairs,
        "stored_entries": unary + pairs,
        "pair_scopes": scope_sizes,
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "missing_pairs_are_unconstrained_cartesian_products": True,
    }
