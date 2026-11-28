"""Triplet-cluster factor graph for bounded W5 KKV support propagation.

Stores unary, pair, and all ten three-node clique relations per spin sector.
This is still an over-approximation (no four-/five-node factors and no
amplitudes), but it preserves selected higher-order compatibility.
"""
from __future__ import annotations
from itertools import combinations, product

import bqg_w5_su2_structural_basis as SB
import bqg_w5_pairwise_factor_graph as PG
import k5_peter_weyl_safe_hda_column as PW

NODES = tuple(range(5))
PAIRS = tuple(combinations(NODES, 2))
TRIPLES = tuple(combinations(NODES, 3))


def _count(s):
    return (sum(map(len, s["unary"])) +
            sum(map(len, s["pairs"].values())) +
            sum(map(len, s["triples"].values())))


def _guard(support, limit, stage):
    n = sum(_count(s) for s in support.values())
    if limit is not None and n > limit:
        raise RuntimeError(f"triplet-cluster {stage} entries {n} exceed guard {limit}; no truncation")


def _merge(dst, spins, state):
    if spins not in dst:
        dst[spins] = {
            "unary": tuple(frozenset(x) for x in state["unary"]),
            "pairs": {p: frozenset(x) for p, x in state["pairs"].items()},
            "triples": {t: frozenset(x) for t, x in state["triples"].items()},
        }
    else:
        old = dst[spins]
        dst[spins] = {
            "unary": tuple(a | b for a, b in zip(old["unary"], state["unary"])),
            "pairs": {p: old["pairs"][p] | state["pairs"][p] for p in PAIRS},
            "triples": {t: old["triples"][t] | state["triples"][t] for t in TRIPLES},
        }


def _consistency(state):
    """Generalized consistency of all unary, pair, and triplet factors."""
    unary = [set(x) for x in state["unary"]]
    pairs = {p: set(v) for p, v in state["pairs"].items()}
    triples = {t: set(v) for t, v in state["triples"].items()}
    changed = True
    while changed:
        changed = False
        # Pairwise arc/path consistency is a safe first pass.
        base = PG._arc_consistency({
            "unary": tuple(frozenset(x) for x in unary),
            "pairs": {p: frozenset(v) for p, v in pairs.items()},
        })
        for i in NODES:
            if set(base["unary"][i]) != unary[i]:
                unary[i] = set(base["unary"][i]); changed = True
        for p in PAIRS:
            if set(base["pairs"][p]) != pairs[p]:
                pairs[p] = set(base["pairs"][p]); changed = True

        # Every triplet tuple must respect all unary and pair factors.
        for t in TRIPLES:
            a, b, c = t
            keep = {(x, y, z) for x, y, z in triples[t]
                    if x in unary[a] and y in unary[b] and z in unary[c]
                    and (x, y) in pairs[(a, b)]
                    and (x, z) in pairs[(a, c)]
                    and (y, z) in pairs[(b, c)]}
            if keep != triples[t]:
                triples[t] = keep; changed = True

        # Pair relations must be supported by every triplet containing them.
        for p in PAIRS:
            a, b = p
            keep = set(pairs[p])
            for t in TRIPLES:
                if a not in t or b not in t:
                    continue
                ia, ib = t.index(a), t.index(b)
                projection = set()
                for row in triples[t]:
                    projection.add((row[ia], row[ib]))
                keep &= projection
            if keep != pairs[p]:
                pairs[p] = keep; changed = True

        # Reproject unary domains from the tightened pair relations.
        for a, b in PAIRS:
            na = {x for x, y in pairs[(a, b)]}
            nb = {y for x, y in pairs[(a, b)]}
            if na != unary[a]:
                unary[a] &= na; changed = True
            if nb != unary[b]:
                unary[b] &= nb; changed = True

    return {"unary": tuple(frozenset(x) for x in unary),
            "pairs": {p: frozenset(v) for p, v in pairs.items()},
            "triples": {t: frozenset(v) for t, v in triples.items()}}


def _initial(key):
    spins, labels = SB.lift_key_to_all_charged(key)
    return {spins: {
        "unary": tuple(frozenset((x,)) for x in labels),
        "pairs": {p: frozenset(((labels[p[0]], labels[p[1]]),)) for p in PAIRS},
        "triples": {t: frozenset(((labels[t[0]], labels[t[1]], labels[t[2]]),))
                    for t in TRIPLES},
    }}


def _lift(support):
    out = {}
    for spins, s in support.items():
        maps = [{x: (("C", 0, 0, x[1], x[1]) if x[0] == "G" else x)
                 for x in labs} for labs in s["unary"]]
        q = {
            "unary": tuple(frozenset(m.values()) for m in maps),
            "pairs": {p: frozenset((maps[p[0]][x], maps[p[1]][y])
                                   for x, y in rel) for p, rel in s["pairs"].items()},
            "triples": {t: frozenset((maps[t[0]][x], maps[t[1]][y], maps[t[2]][z])
                                     for x, y, z in rel)
                        for t, rel in s["triples"].items()},
        }
        _merge(out, spins, q)
    return out


def _project(support, charged, limit):
    charged = set(charged)
    out = {}
    for spins, s in support.items():
        maps = []
        for v, labs in enumerate(s["unary"]):
            mp = {}
            for x in labs:
                if v in charged:
                    if x[0] == "C": mp[x] = x
                elif x[0] == "C" and x[1:3] == (0, 0) and x[3] == x[4]:
                    mp[x] = ("G", x[3])
            maps.append(mp)
        q = {
            "unary": tuple(frozenset(m.values()) for m in maps),
            "pairs": {p: frozenset((maps[p[0]][x], maps[p[1]][y])
                                   for x, y in rel if x in maps[p[0]] and y in maps[p[1]])
                      for p, rel in s["pairs"].items()},
            "triples": {t: frozenset((maps[t[0]][x], maps[t[1]][y], maps[t[2]][z])
                                     for x, y, z in rel
                                     if x in maps[t[0]] and y in maps[t[1]] and z in maps[t[2]])
                        for t, rel in s["triples"].items()},
        }
        q = _consistency(q)
        if all(q["unary"]):
            _merge(out, spins, q); _guard(out, limit, "projection")
    return out


def _volume(support, vertex, limit):
    out = {}
    for spins, s in support.items():
        mp = {}
        for old in s["unary"][vertex]:
            mp[old] = frozenset(SB.volume_label_targets(SB.local_spins(spins, vertex), old))
        unary = list(s["unary"])
        unary[vertex] = frozenset(z for vals in mp.values() for z in vals)
        pairs = {}
        for p, rel in s["pairs"].items():
            if vertex == p[0]:
                pairs[p] = frozenset((z, y) for x, y in rel for z in mp.get(x, ()))
            elif vertex == p[1]:
                pairs[p] = frozenset((x, z) for x, y in rel for z in mp.get(y, ()))
            else: pairs[p] = rel
        triples = {}
        for t, rel in s["triples"].items():
            if vertex not in t: triples[t] = rel; continue
            ix = t.index(vertex)
            rows = set()
            for row in rel:
                for z in mp.get(row[ix], ()):
                    new = list(row); new[ix] = z; rows.add(tuple(new))
            triples[t] = frozenset(rows)
        q = _consistency({"unary": tuple(unary), "pairs": pairs, "triples": triples})
        if all(q["unary"]):
            _merge(out, spins, q); _guard(out, limit, "volume")
    return out


def _hit(support, x, y, Jmax2, limit):
    edge = tuple(sorted((x, y)))
    if x == y or edge not in SB.EDGE_INDEX: raise ValueError("invalid K5 edge")
    ei = SB.EDGE_INDEX[edge]
    out = {}
    for spins, s in support.items():
        for so in (spins[ei]-1, spins[ei]+1):
            if so < 0 or so > Jmax2: continue
            ns = list(spins); ns[ei] = so; ns = tuple(ns)
            trans = {}
            for v in (x, y):
                possible = set()
                for old in s["unary"][v]:
                    oldj = old[1] if old[0] == "C" else 0
                    possible.update(j for j in SB.all_total_J2(SB.local_spins(ns, v))
                                    if abs(j-oldj) <= 1 and j+oldj >= 1)
                basis = SB.charged_node_basis(SB.local_spins(ns, v))
                trans[v] = {old: frozenset(z for z in basis if z[1] in possible)
                            for old in s["unary"][v]}
            unary = list(s["unary"])
            for v in (x, y):
                unary[v] = frozenset(z for vals in trans[v].values() for z in vals)
            pairs = {}
            for p, rel in s["pairs"].items():
                a, b = p
                pairs[p] = frozenset((u, v)
                    for olda, oldb in rel
                    for u in (trans[a].get(olda, ()) if a in trans else (olda,))
                    for v in (trans[b].get(oldb, ()) if b in trans else (oldb,)))
            triples = {}
            for t, rel in s["triples"].items():
                rows = set()
                for row in rel:
                    opts = [trans[v].get(row[i], (row[i],)) if v in trans else (row[i],)
                            for i, v in enumerate(t)]
                    rows.update(product(*opts))
                triples[t] = frozenset(rows)
            q = _consistency({"unary": tuple(unary), "pairs": pairs, "triples": triples})
            if all(q["unary"]):
                _merge(out, ns, q); _guard(out, limit, "holonomy")
    return out


def _sequence(support, sequence, source_v, Jmax2, charged, limit):
    work = _lift(support)
    for op in sequence:
        if op[0] == "V": work = _volume(work, source_v, limit)
        elif op[0] == "P":
            path = tuple(op[1])
            for x, y in zip(path[:-1], path[1:]):
                work = _hit(work, x, y, Jmax2, limit)
        else: raise ValueError(f"unknown operator {op[0]}")
        if not work: break
    return _project(work, charged, limit)


def _he(support, source_v, Jmax2, charged, limit):
    out = {}
    for _sign, spec in PW.oriented_specs(source_v):
        v, a, b, c = spec
        for _coef, seq in PW.T_sequences(v, a, b, c):
            for word in (seq, PW.adjoint_sequence(seq)):
                part = _sequence(support, word, source_v, Jmax2, charged, limit)
                for spins, state in part.items():
                    _merge(out, spins, state)
                    out[spins] = _consistency(out[spins])
                    _guard(out, limit, "H_E")
    return out


def _hit_project(support, x, y, Jmax2, charged, limit):
    return _project(_hit(_lift(support), x, y, Jmax2, limit), charged, limit)


def _k(support, source_v, Jmax2, charged, limit):
    he = _he(support, source_v, Jmax2, charged, limit)
    vh = _volume(he, source_v, limit)
    vin = _project(_volume(_lift(support), source_v, limit), charged, limit)
    hv = _he(vin, source_v, Jmax2, charged, limit)
    out = {}
    for q in (vh, hv):
        for spins, state in q.items():
            _merge(out, spins, state); out[spins] = _consistency(out[spins])
            _guard(out, limit, "K")
    return out


def _c_k(support, source_v, target_v, Jmax2, i, j, limit):
    out = {}
    if i == j: out = _k(support, source_v, Jmax2, (source_v,), limit)
    inverse = _hit_project(support, target_v, source_v, Jmax2,
                           (source_v, target_v), limit)
    middle = _k(inverse, source_v, Jmax2, (source_v, target_v), limit)
    closed = _hit_project(middle, source_v, target_v, Jmax2, (source_v,), limit)
    for q in (out, closed):
        for spins, state in q.items():
            _merge(out, spins, state); out[spins] = _consistency(out[spins])
            _guard(out, limit, "C(K)")
    return out


def _c_v(support, source_v, target_v, Jmax2, i, j, limit):
    out = _volume(support, source_v, limit) if i == j else {}
    inverse = _hit_project(support, target_v, source_v, Jmax2,
                           (source_v, target_v), limit)
    middle = _volume(inverse, source_v, limit)
    closed = _hit_project(middle, source_v, target_v, Jmax2, (source_v,), limit)
    for spins, state in closed.items():
        _merge(out, spins, state); out[spins] = _consistency(out[spins])
        _guard(out, limit, "C(V)")
    return out


def _scalar_filter(support, source_v):
    out = {}
    for spins, s in support.items():
        allowed = frozenset(x for x in s["unary"][source_v]
                            if x[0] == "C" and x[1] in (0, 2))
        if not allowed: continue
        unary = list(s["unary"]); unary[source_v] = allowed
        pairs = {p: frozenset((x, y) for x, y in rel
                              if (p[0] != source_v or x in allowed)
                              and (p[1] != source_v or y in allowed))
                 for p, rel in s["pairs"].items()}
        triples = {t: frozenset(row for row in rel
                                if (t[0] != source_v or row[0] in allowed)
                                and (t[1] != source_v or row[1] in allowed)
                                and (t[2] != source_v or row[2] in allowed))
                   for t, rel in s["triples"].items()}
        q = _consistency({"unary": tuple(unary), "pairs": pairs, "triples": triples})
        if all(q["unary"]): _merge(out, spins, q)
    return out


def factor_graph_triple(support, source_v, a, b, c, Jmax2,
                        adjoint=False, limit=100000):
    out = {}
    for i, j, k in product(range(2), repeat=3):
        if not adjoint:
            q = _c_v(support, source_v, c, Jmax2, k, i, limit)
            q = _c_k(q, source_v, b, Jmax2, j, k, limit)
            q = _scalar_filter(q, source_v)
            q = _c_k(q, source_v, a, Jmax2, i, j, limit)
        else:
            q = _c_k(support, source_v, a, Jmax2, j, i, limit)
            q = _c_k(q, source_v, b, Jmax2, k, j, limit)
            q = _scalar_filter(q, source_v)
            q = _c_v(q, source_v, c, Jmax2, i, k, limit)
        for spins, state in q.items():
            _merge(out, spins, state); out[spins] = _consistency(out[spins])
            _guard(out, limit, "KKV")
    return out


def contains_joint(support, key):
    spins, labels = key
    if spins not in support: return False
    s = support[spins]
    if any(labels[v] not in s["unary"][v] for v in NODES): return False
    if not all((labels[a], labels[b]) in s["pairs"][(a,b)] for a,b in PAIRS): return False
    return all((labels[a], labels[b], labels[c]) in s["triples"][(a,b,c)]
               for a,b,c in TRIPLES)


def profile(support):
    return {
        "spin_sector_count": len(support),
        "unary_entries": sum(len(x) for s in support.values() for x in s["unary"]),
        "pair_entries": sum(len(x) for s in support.values() for x in s["pairs"].values()),
        "triplet_entries": sum(len(x) for s in support.values() for x in s["triples"].values()),
        "stored_entries": sum(_count(s) for s in support.values()),
        "product_bound": sum(__import__("math").prod(len(x) for x in s["unary"])
                             for s in support.values()),
    }
