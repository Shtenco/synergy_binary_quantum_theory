"""Pairwise factor-graph support propagation for W5 KKV.

Each spin sector stores five unary label domains and ten binary relations.
This retains pair correlations while remaining an over-approximation of the
joint support. Arc consistency is applied after every local update. No
amplitudes, exact joint states, or physical non-vanishing claims are made.
"""
from __future__ import annotations

from itertools import combinations, product

import bqg_w5_su2_structural_basis as SB
import k5_peter_weyl_safe_hda_column as PW

NODES = tuple(range(5))
PAIRS = tuple(combinations(NODES, 2))


def _empty():
    return {"unary": tuple(frozenset() for _ in NODES),
            "pairs": {p: frozenset() for p in PAIRS}}


def _size(state):
    return sum(map(len, state["unary"])) + sum(map(len, state["pairs"].values()))


def _guard(support, limit, stage):
    count = sum(_size(s) for s in support.values())
    if limit is not None and count > limit:
        raise RuntimeError(f"pairwise {stage} entries {count} exceed guard {limit}; no truncation")


def _merge(dst, spins, state):
    if spins not in dst:
        dst[spins] = {
            "unary": tuple(frozenset(x) for x in state["unary"]),
            "pairs": {p: frozenset(x) for p, x in state["pairs"].items()},
        }
        return
    old = dst[spins]
    dst[spins] = {
        "unary": tuple(a | b for a, b in zip(old["unary"], state["unary"])),
        "pairs": {p: old["pairs"][p] | state["pairs"][p] for p in PAIRS},
    }


def _arc_consistency(state):
    """Enforce unary arc consistency and binary path consistency.

    For every node triple (i,j,k), a retained pair (x_i,x_j) must have at
    least one common supporting label at k. This is a safe tightening of the
    pairwise relaxation: every genuinely joint-supported pair survives.
    """
    unary = [set(x) for x in state["unary"]]
    pairs = {p: set(v) for p, v in state["pairs"].items()}

    def oriented(a, b):
        p = tuple(sorted((a, b)))
        rel = pairs[p]
        if a < b:
            return rel
        return {(y, x) for x, y in rel}

    changed = True
    while changed:
        changed = False
        # Unary arc consistency.
        for (a, b), rel in pairs.items():
            na = {x for x, y in rel if y in unary[b]}
            nb = {y for x, y in rel if x in unary[a]}
            if na != unary[a]:
                unary[a] &= na
                changed = True
            if nb != unary[b]:
                unary[b] &= nb
                changed = True
            nr = {(x, y) for x, y in rel
                  if x in unary[a] and y in unary[b]}
            if nr != rel:
                pairs[(a, b)] = nr
                changed = True

        # Binary path consistency across every third node.
        for i, j in PAIRS:
            rel_ij = pairs[(i, j)]
            if not rel_ij:
                continue
            for k in NODES:
                if k == i or k == j:
                    continue
                rik = oriented(i, k)
                rjk = oriented(j, k)
                left = {}
                right = {}
                for x, y in rik:
                    left.setdefault(x, set()).add(y)
                for z, y in rjk:
                    right.setdefault(z, set()).add(y)
                kept = {
                    (x, z) for x, z in rel_ij
                    if left.get(x, set()).intersection(right.get(z, set()))
                }
                if kept != rel_ij:
                    pairs[(i, j)] = kept
                    rel_ij = kept
                    changed = True

    return {"unary": tuple(frozenset(x) for x in unary),
            "pairs": {p: frozenset(v) for p, v in pairs.items()}}


def _union(*supports, limit=None, stage="union"):
    out = {}
    for support in supports:
        for spins, state in support.items():
            _merge(out, spins, state)
            out[spins] = _arc_consistency(out[spins])
            _guard(out, limit, stage)
    return out


def _initial(key):
    spins, labels = SB.lift_key_to_all_charged(key)
    state = {
        "unary": tuple(frozenset((x,)) for x in labels),
        "pairs": {p: frozenset(((labels[p[0]], labels[p[1]]),)) for p in PAIRS},
    }
    return {spins: state}


def _lift(support):
    out = {}
    for spins, state in support.items():
        unary = []
        maps = []
        for labs in state["unary"]:
            mp = {x: (("C", 0, 0, x[1], x[1]) if x[0] == "G" else x)
                  for x in labs}
            maps.append(mp)
            unary.append(frozenset(mp.values()))
        pairs = {
            p: frozenset((maps[p[0]][x], maps[p[1]][y])
                         for x, y in rel)
            for p, rel in state["pairs"].items()
        }
        _merge(out, spins, {"unary": tuple(unary), "pairs": pairs})
    return out


def _project(support, charged_nodes, limit):
    charged = set(charged_nodes)
    out = {}
    for spins, state in support.items():
        maps = []
        unary = []
        for v, labs in enumerate(state["unary"]):
            mp = {}
            for lab in labs:
                if v in charged:
                    if lab[0] == "C":
                        mp[lab] = lab
                elif lab[0] == "C" and lab[1:3] == (0, 0) and lab[3] == lab[4]:
                    mp[lab] = ("G", lab[3])
            maps.append(mp)
            unary.append(frozenset(mp.values()))
        pairs = {
            p: frozenset((maps[p[0]][x], maps[p[1]][y])
                         for x, y in rel if x in maps[p[0]] and y in maps[p[1]])
            for p, rel in state["pairs"].items()
        }
        projected = _arc_consistency({"unary": tuple(unary), "pairs": pairs})
        if all(projected["unary"]):
            _merge(out, spins, projected)
            _guard(out, limit, "projection")
    return out


def _volume(support, vertex, limit):
    out = {}
    for spins, state in support.items():
        mp = {}
        for old in state["unary"][vertex]:
            for new in SB.volume_label_targets(SB.local_spins(spins, vertex), old):
                mp.setdefault(old, set()).add(new)
        unary = list(state["unary"])
        unary[vertex] = frozenset(y for vals in mp.values() for y in vals)
        pairs = {}
        for p, rel in state["pairs"].items():
            if vertex not in p:
                pairs[p] = rel
            elif p[0] == vertex:
                pairs[p] = frozenset((z, y) for x, y in rel
                                     for z in mp.get(x, ()))
            else:
                pairs[p] = frozenset((x, z) for x, y in rel
                                     for z in mp.get(y, ()))
        q = _arc_consistency({"unary": tuple(unary), "pairs": pairs})
        if all(q["unary"]):
            _merge(out, spins, q)
            _guard(out, limit, "volume")
    return out


def _hit(support, x, y, Jmax2, limit):
    edge = tuple(sorted((x, y)))
    if x == y or edge not in SB.EDGE_INDEX:
        raise ValueError("invalid K5 edge")
    ei = SB.EDGE_INDEX[edge]
    out = {}
    for spins, state in support.items():
        for so in (spins[ei] - 1, spins[ei] + 1):
            if so < 0 or so > Jmax2:
                continue
            ns = list(spins)
            ns[ei] = so
            ns = tuple(ns)
            transitions = {}
            for v in (x, y):
                possible_J = set()
                for old in state["unary"][v]:
                    old_j = old[1] if old[0] == "C" else 0
                    possible_J.update(
                        j for j in SB.all_total_J2(SB.local_spins(ns, v))
                        if abs(j - old_j) <= 1 and j + old_j >= 1
                    )
                transitions[v] = {
                    old: frozenset(lab for lab in SB.charged_node_basis(
                        SB.local_spins(ns, v)) if lab[1] in possible_J)
                    for old in state["unary"][v]
                }
            unary = list(state["unary"])
            for v in (x, y):
                unary[v] = frozenset(z for vals in transitions[v].values() for z in vals)
            pairs = {}
            for p, rel in state["pairs"].items():
                a, b = p
                if a in transitions and b in transitions:
                    pairs[p] = frozenset((u, v) for olda, oldb in rel
                        for u in transitions[a].get(olda, ())
                        for v in transitions[b].get(oldb, ()))
                elif a in transitions:
                    pairs[p] = frozenset((u, v) for olda, v in rel
                        for u in transitions[a].get(olda, ()))
                elif b in transitions:
                    pairs[p] = frozenset((u, v) for u, oldb in rel
                        for v in transitions[b].get(oldb, ()))
                else:
                    pairs[p] = rel
            q = _arc_consistency({"unary": tuple(unary), "pairs": pairs})
            if all(q["unary"]):
                _merge(out, ns, q)
                _guard(out, limit, "holonomy")
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
            raise ValueError(f"unknown operator {op[0]}")
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
            out = _union(out, _sequence(support, word, source_v, Jmax2,
                                        charged_nodes, limit),
                         limit=limit, stage="T sequence")
    return out


def factor_graph_he(support, source_v, Jmax2, charged_nodes=(0, 1),
                    limit=100000):
    out = {}
    for _sign, spec in PW.oriented_specs(source_v):
        v, a, b, c = spec
        out = _union(out, _t_sequence(support, v, a, b, c, source_v,
                                      Jmax2, charged_nodes, limit, True),
                     limit=limit, stage="H_E")
    return out


def _hit_project(support, x, y, Jmax2, charged, limit):
    return _project(_hit(_lift(support), x, y, Jmax2, limit), charged, limit)


def factor_graph_k(support, source_v, Jmax2, charged=(0,), limit=100000):
    he = factor_graph_he(support, source_v, Jmax2, charged, limit)
    vh = _volume(he, source_v, limit)
    vin = _project(_volume(_lift(support), source_v, limit), charged, limit)
    hv = factor_graph_he(vin, source_v, Jmax2, charged, limit)
    return _union(vh, hv, limit=limit, stage="K")


def factor_graph_c_k(support, source_v, target_v, Jmax2, i, j, limit=100000):
    out = {}
    if i == j:
        out = factor_graph_k(support, source_v, Jmax2, (source_v,), limit)
    inverse = _hit_project(support, target_v, source_v, Jmax2,
                           (source_v, target_v), limit)
    middle = factor_graph_k(inverse, source_v, Jmax2,
                            (source_v, target_v), limit)
    closed = _hit_project(middle, source_v, target_v, Jmax2,
                          (source_v,), limit)
    return _union(out, closed, limit=limit, stage="C(K)")


def _c_volume(support, source_v, target_v, Jmax2, i, j, limit):
    out = _volume(support, source_v, limit) if i == j else {}
    inverse = _hit_project(support, target_v, source_v, Jmax2,
                           (source_v, target_v), limit)
    middle = _volume(inverse, source_v, limit)
    closed = _hit_project(middle, source_v, target_v, Jmax2,
                          (source_v,), limit)
    return _union(out, closed, limit=limit, stage="C(V)")


def _scalar_filter(support, source_v):
    out = {}
    for spins, state in support.items():
        allowed = frozenset(x for x in state["unary"][source_v]
                            if x[0] == "C" and x[1] in (0, 2))
        if not allowed:
            continue
        unary = list(state["unary"])
        unary[source_v] = allowed
        pairs = {p: frozenset((x, y) for x, y in rel
                              if (p[0] != source_v or x in allowed)
                              and (p[1] != source_v or y in allowed))
                 for p, rel in state["pairs"].items()}
        q = _arc_consistency({"unary": tuple(unary), "pairs": pairs})
        if all(q["unary"]):
            _merge(out, spins, q)
    return out


def factor_graph_triple(support, source_v, a, b, c, Jmax2,
                        adjoint=False, limit=100000):
    out = {}
    for i, j, k in product(range(2), repeat=3):
        if not adjoint:
            q = _c_volume(support, source_v, c, Jmax2, k, i, limit)
            q = factor_graph_c_k(q, source_v, b, Jmax2, j, k, limit)
            q = _scalar_filter(q, source_v)
            q = factor_graph_c_k(q, source_v, a, Jmax2, i, j, limit)
        else:
            q = factor_graph_c_k(support, source_v, a, Jmax2, j, i, limit)
            q = factor_graph_c_k(q, source_v, b, Jmax2, k, j, limit)
            q = _scalar_filter(q, source_v)
            q = _c_volume(q, source_v, c, Jmax2, i, k, limit)
        out = _union(out, q, limit=limit, stage="KKV")
    return out


def contains_joint(support, key):
    spins, labels = key
    if spins not in support:
        return False
    state = support[spins]
    if any(labels[v] not in state["unary"][v] for v in NODES):
        return False
    return all((labels[a], labels[b]) in state["pairs"][(a, b)]
               for a, b in PAIRS)


def profile(support):
    sectors = len(support)
    unary = sum(len(x) for s in support.values() for x in s["unary"])
    pair_entries = sum(len(x) for s in support.values() for x in s["pairs"].values())
    product_bound = sum(
        __import__("math").prod(len(x) for x in s["unary"])
        for s in support.values()
    )
    return {"spin_sector_count": sectors, "unary_entries": unary,
            "pair_relation_entries": pair_entries,
            "stored_entries": unary + pair_entries,
            "pairwise_product_bound": product_bound}
