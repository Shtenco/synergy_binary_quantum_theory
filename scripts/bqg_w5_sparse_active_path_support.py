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


def active_pair_scope(sequence):
    """Directly hit edges in a word; do not complete them to a clique."""
    active = set()
    for op in sequence:
        if op[0] == "P":
            path = tuple(op[1])
            active.update(tuple(sorted(edge)) for edge in
                          zip(path[:-1], path[1:]))
        elif op[0] != "V":
            raise ValueError(f"unknown operator {op[0]}")
    return frozenset(active)


def _stored_pairs(scope):
    return tuple(combinations(scope, 2))


def _scope_of(support, extra=()):
    scope = set(extra)
    for state in support.values():
        for a, b in state["pairs"]:
            scope.update((a, b))
    return tuple(sorted(scope))


def _pair_scope_of(support):
    return frozenset(pair for state in support.values()
                     for pair in state["pairs"])


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


def _compress_pairs(support, pair_scope):
    """Retain nontrivial factors on explicitly active/inherited edges only."""
    keep = set(pair_scope)
    out = {}
    for spins, state in support.items():
        pairs = {}
        for p, rel in state["pairs"].items():
            if p not in keep:
                continue
            relation = frozenset(rel)
            cartesian = frozenset(
                (x, y) for x in state["unary"][p[0]]
                for y in state["unary"][p[1]]
            )
            if relation != cartesian:
                pairs[p] = relation
        out[spins] = {
            "unary": tuple(frozenset(x) for x in state["unary"]),
            "pairs": pairs,
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


def _sparse_arc_consistency(state):
    """Arc/path consistency with omitted pairs interpreted as Cartesian.

    Unlike PG._arc_consistency this never materializes all ten K5 relations.
    """
    unary = [set(x) for x in state["unary"]]
    pairs = {tuple(p): set(rel) for p, rel in state["pairs"].items()}

    def oriented(a, b):
        p = tuple(sorted((a, b)))
        if p in pairs:
            rel = pairs[p]
            return rel if a < b else {(y, x) for x, y in rel}
        return {(x, y) for x in unary[a] for y in unary[b]}

    changed = True
    while changed:
        changed = False
        for (a, b), rel in list(pairs.items()):
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

        for i, j in tuple(pairs):
            rel_ij = pairs[(i, j)]
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


def _sparse_lift(support):
    """Lift Gauss labels without creating omitted pair factors."""
    out = {}
    for spins, state in support.items():
        maps = []
        unary = []
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
        out[spins] = {"unary": tuple(unary), "pairs": pairs}
    return out


def _ensure_pair_scope(support, pair_scope, limit):
    """Materialize only directly hit or inherited pair factors."""
    out = {}
    for spins, state in support.items():
        pairs = dict(state["pairs"])
        for a, b in pair_scope:
            if (a, b) not in pairs:
                pairs[(a, b)] = frozenset(
                    (x, y) for x in state["unary"][a]
                    for y in state["unary"][b]
                )
        out[spins] = {"unary": state["unary"], "pairs": pairs}
    _guard(out, limit, "active-pair materialization")
    return out


def _sparse_volume(support, vertex, limit):
    out = {}
    for spins, state in support.items():
        transitions = {}
        for old in state["unary"][vertex]:
            transitions[old] = frozenset(
                SB.volume_label_targets(SB.local_spins(spins, vertex), old)
            )
        unary = list(state["unary"])
        unary[vertex] = frozenset(
            new for values in transitions.values() for new in values
        )
        pairs = {}
        for p, rel in state["pairs"].items():
            if vertex not in p:
                pairs[p] = rel
            elif p[0] == vertex:
                pairs[p] = frozenset(
                    (new, other) for old, other in rel
                    for new in transitions.get(old, ())
                )
            else:
                pairs[p] = frozenset(
                    (other, new) for other, old in rel
                    for new in transitions.get(old, ())
                )
        q = _sparse_arc_consistency({"unary": tuple(unary), "pairs": pairs})
        if all(q["unary"]):
            out[spins] = q
            _guard(out, limit, "sparse volume")
    return out


def _sparse_hit(support, x, y, Jmax2, limit):
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
                    old: frozenset(
                        lab for lab in SB.charged_node_basis(
                            SB.local_spins(ns, v))
                        if lab[1] in possible_J
                    )
                    for old in state["unary"][v]
                }
            unary = list(state["unary"])
            for v in (x, y):
                unary[v] = frozenset(
                    new for values in transitions[v].values() for new in values
                )
            pairs = {}
            for p, rel in state["pairs"].items():
                a, b = p
                if a in transitions and b in transitions:
                    pairs[p] = frozenset(
                        (u, v) for old_a, old_b in rel
                        for u in transitions[a].get(old_a, ())
                        for v in transitions[b].get(old_b, ())
                    )
                elif a in transitions:
                    pairs[p] = frozenset(
                        (u, v) for old_a, v in rel
                        for u in transitions[a].get(old_a, ())
                    )
                elif b in transitions:
                    pairs[p] = frozenset(
                        (u, v) for u, old_b in rel
                        for v in transitions[b].get(old_b, ())
                    )
                else:
                    pairs[p] = rel
            q = _sparse_arc_consistency(
                {"unary": tuple(unary), "pairs": pairs}
            )
            if all(q["unary"]):
                _merge_sparse(out, ns, q)
                _guard(out, limit, "sparse holonomy")
    return out


def _merge_sparse(dst, spins, state):
    if spins not in dst:
        dst[spins] = {
            "unary": tuple(frozenset(x) for x in state["unary"]),
            "pairs": {p: frozenset(x) for p, x in state["pairs"].items()},
        }
        return
    old = dst[spins]
    unary = tuple(a | b for a, b in zip(old["unary"], state["unary"]))
    keys = set(old["pairs"]) | set(state["pairs"])
    pairs = {}
    for p in keys:
        if p not in old["pairs"] or p not in state["pairs"]:
            pairs[p] = frozenset(
                (x, y) for x in unary[p[0]] for y in unary[p[1]]
            )
        else:
            pairs[p] = old["pairs"][p] | state["pairs"][p]
    dst[spins] = _sparse_arc_consistency({"unary": unary, "pairs": pairs})


def _sparse_project(support, charged_nodes, limit):
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
                elif (lab[0] == "C" and lab[1:3] == (0, 0)
                      and lab[3] == lab[4]):
                    mp[lab] = ("G", lab[3])
            maps.append(mp)
            unary.append(frozenset(mp.values()))
        pairs = {
            p: frozenset(
                (maps[p[0]][x], maps[p[1]][y])
                for x, y in rel
                if x in maps[p[0]] and y in maps[p[1]]
            )
            for p, rel in state["pairs"].items()
        }
        projected = _sparse_arc_consistency(
            {"unary": tuple(unary), "pairs": pairs}
        )
        if all(projected["unary"]):
            _merge_sparse(out, spins, projected)
            _guard(out, limit, "sparse projection")
    return out


def _sparse_scalar_filter(support, source_v, limit):
    out = {}
    for spins, state in support.items():
        allowed = frozenset(
            x for x in state["unary"][source_v]
            if x[0] == "C" and x[1] in (0, 2)
        )
        if not allowed:
            continue
        unary = list(state["unary"])
        unary[source_v] = allowed
        pairs = {
            p: frozenset(
                (x, y) for x, y in rel
                if (p[0] != source_v or x in allowed)
                and (p[1] != source_v or y in allowed)
            )
            for p, rel in state["pairs"].items()
        }
        q = _sparse_arc_consistency(
            {"unary": tuple(unary), "pairs": pairs}
        )
        if all(q["unary"]):
            out[spins] = q
            _guard(out, limit, "sparse scalar filter")
    return out


def _guard(support, limit, stage):
    size = sum(
        sum(map(len, s["unary"])) + sum(map(len, s["pairs"].values()))
        for s in support.values()
    )
    if limit is not None and size > limit:
        unary_entries = sum(
            sum(map(len, state["unary"])) for state in support.values()
        )
        pair_entries = sum(
            sum(map(len, state["pairs"].values())) for state in support.values()
        )
        pair_factor_counts = sorted({
            len(state["pairs"]) for state in support.values()
        })
        pair_entry_counts = {
            f"{a}-{b}": sum(
                len(state["pairs"].get((a, b), ()))
                for state in support.values()
            )
            for a, b in combinations(NODES, 2)
            if any((a, b) in state["pairs"] for state in support.values())
        }
        raise RuntimeError(
            f"sparse-active-path {stage} entries {size} exceed guard {limit}; "
            f"unary_entries={unary_entries}; pair_entries={pair_entries}; "
            f"pair_factor_counts={pair_factor_counts}; "
            f"pair_entry_counts={pair_entry_counts}; no truncation"
        )


def sequence(support, word, source_v, Jmax2, charged_nodes=(0,),
             limit=100000):
    """Propagate one primitive word with active-scope pair storage."""
    pair_scope = set()
    # Preserve inherited nontrivial relations; redundant singleton seed
    # factors do not activate all ten K5 edges.
    for state in support.values():
        for (a, b), relation in state["pairs"].items():
            cartesian = frozenset((x, y) for x in state["unary"][a]
                                  for y in state["unary"][b])
            if relation != cartesian:
                pair_scope.add((a, b))
    work = _compress_pairs(_sparse_lift(support), pair_scope)
    _guard(work, limit, "initial active pair scope")
    for op in word:
        if op[0] == "V":
            pass
        elif op[0] == "P":
            path = tuple(op[1])
            pair_scope.update(tuple(sorted(edge)) for edge in
                              zip(path[:-1], path[1:]))
        else:
            raise ValueError(f"unknown operator {op[0]}")
        work = _ensure_pair_scope(work, pair_scope, limit)
        if op[0] == "V":
            work = _sparse_volume(work, source_v, limit)
        else:
            path = tuple(op[1])
            for x, y in zip(path[:-1], path[1:]):
                try:
                    work = _sparse_hit(work, x, y, Jmax2, limit)
                except RuntimeError as exc:
                    raise RuntimeError(
                        f"{exc}; operator={op!r}; edge={(x, y)!r}; "
                        f"active_pair_scope={tuple(sorted(pair_scope))!r}"
                    ) from exc
        work = _compress_pairs(work, pair_scope)
        _guard(work, limit, f"after {op[0]}")
        if not work:
            break
    projected = _sparse_project(work, charged_nodes, limit)
    return _compress_pairs(projected, pair_scope)


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


def _union(left, right, limit, stage):
    """Conservative disjunctive union in the sparse representation."""
    out = {s: {"unary": tuple(frozenset(x) for x in q["unary"]),
               "pairs": dict(q["pairs"])} for s, q in left.items()}
    for spins, state in right.items():
        if spins not in out:
            out[spins] = {"unary": tuple(frozenset(x) for x in state["unary"]),
                          "pairs": dict(state["pairs"])}
            continue
        old = out[spins]
        unary = tuple(a | b for a, b in zip(old["unary"], state["unary"]))
        pair_keys = set(old["pairs"]) | set(state["pairs"])
        pairs = {}
        for p in pair_keys:
            # A missing relation means unconstrained, so its union with any
            # explicit relation is the full Cartesian product.
            if p not in old["pairs"] or p not in state["pairs"]:
                pairs[p] = frozenset((x, y) for x in unary[p[0]]
                                     for y in unary[p[1]])
            else:
                pairs[p] = old["pairs"][p] | state["pairs"][p]
        out[spins] = {"unary": unary, "pairs": pairs}
    _guard(out, limit, stage)
    return out


def he_spec_branches(branches, source_v, Jmax2, charged_nodes=(0, 1),
                     limit=100000):
    """H_E alternatives: preserve oriented specs; union T/adjoint words per spec."""
    import k5_peter_weyl_safe_hda_column as PW
    out = []
    for incoming in branches:
        for _sign, spec in PW.oriented_specs(source_v):
            v, a, b, c = spec
            spec_union = {}
            for _coef, word in PW.T_sequences(v, a, b, c):
                for candidate in (word, PW.adjoint_sequence(word)):
                    part = sequence(incoming, candidate, source_v, Jmax2,
                                    charged_nodes, limit)
                    # Preserve both the current word's active scope and
                    # every inherited nontrivial pair scope. Dropping an
                    # inherited factor is conservative as a relaxation, but
                    # needlessly destroys correlations across nested words.
                    pair_scope = set(active_pair_scope(candidate))
                    pair_scope.update(_pair_scope_of(part))
                    part = _compress_pairs(part, pair_scope)
                    spec_union = _union(spec_union, part, limit,
                                        "sparse H_E T-word union")
            if spec_union:
                out.append(spec_union)
    return out


def k_branches(branches, source_v, Jmax2, charged_nodes=(0,),
               limit=100000):
    """Disjunctive support for K=[V,H_E] = V H_E - H_E V."""
    out = []
    for he in he_spec_branches(branches, source_v, Jmax2,
                               NODES, limit):
        vh = sequence(he, (("V", source_v),), source_v, Jmax2,
                      NODES, limit)
        if vh:
            out.append(_compress_pairs(vh, _pair_scope_of(vh)))
    vin = []
    for branch in branches:
        lifted = _sparse_lift(branch)
        vol = _sparse_volume(lifted, source_v, limit)
        projected = _sparse_project(vol, NODES, limit)
        if projected:
            vin.append(_compress_pairs(projected, _pair_scope_of(projected)))
    out.extend(he_spec_branches(vin, source_v, Jmax2,
                                NODES, limit))
    return out


def c_k_branches(branches, source_v, target_v, Jmax2, i, j,
                 limit=100000):
    """C_e(K) support with inverse then forward holonomy closure."""
    out = []
    if i == j:
        out.extend(k_branches(branches, source_v, Jmax2,
                              (source_v,), limit))
    inverse = []
    for branch in branches:
        part = sequence(branch, (("P", (target_v, source_v)),),
                        source_v, Jmax2, NODES, limit)
        if part:
            inverse.append(_compress_pairs(part, _pair_scope_of(part)))
    middle = k_branches(inverse, source_v, Jmax2,
                        NODES, limit)
    for branch in middle:
        closed = sequence(branch, (("P", (source_v, target_v)),),
                          source_v, Jmax2, NODES, limit)
        if closed:
            out.append(_compress_pairs(closed, _pair_scope_of(closed)))
    return out


def c_volume_branches(branches, source_v, target_v, Jmax2, i, j,
                      limit=100000):
    """C_e(V) support, retaining the diagonal and conjugated alternatives."""
    out = []
    if i == j:
        for branch in branches:
            part = sequence(branch, (("V", source_v),), source_v,
                            Jmax2, NODES, limit)
            if part:
                out.append(_compress_pairs(part, _pair_scope_of(part)))
    inverse = []
    for branch in branches:
        part = sequence(branch, (("P", (target_v, source_v)),),
                        source_v, Jmax2, NODES, limit)
        if part:
            inverse.append(_compress_pairs(part, _pair_scope_of(part)))
    middle = []
    for branch in inverse:
        part = sequence(branch, (("V", source_v),), source_v,
                         Jmax2, NODES, limit)
        if part:
            middle.append(_compress_pairs(part, _pair_scope_of(part)))
    for branch in middle:
        part = sequence(branch, (("P", (source_v, target_v)),),
                        source_v, Jmax2, NODES, limit)
        if part:
            out.append(_compress_pairs(part, _pair_scope_of(part)))
    return out


def scalar_filter(branches, source_v, limit=100000):
    out = []
    for branch in branches:
        filtered = _sparse_scalar_filter(branch, source_v, limit)
        if filtered:
            out.append(_compress_pairs(filtered, _pair_scope_of(branch)))
    _guard(out, limit, "scalar-source filter")
    return out


def triple_branches(branches, source_v, a, b, c, Jmax2,
                    adjoint=False, limit=100000):
    """Branch-preserving structural envelope for forward/reversed-adjoint KKV."""
    from itertools import product as cartesian_product
    out = []
    for i, j, k in cartesian_product(range(2), repeat=3):
        if not adjoint:
            q = c_volume_branches(branches, source_v, c, Jmax2,
                                  k, i, limit)
            next_q = []
            for b0 in q:
                next_q.extend(c_k_branches([b0], source_v, b, Jmax2,
                                           j, k, limit))
            q = scalar_filter(next_q, source_v, limit)
            next_q = []
            for b0 in q:
                next_q.extend(c_k_branches([b0], source_v, a, Jmax2,
                                           i, j, limit))
            q = next_q
        else:
            q = []
            for b0 in branches:
                q.extend(c_k_branches([b0], source_v, a, Jmax2,
                                      j, i, limit))
            next_q = []
            for b0 in q:
                next_q.extend(c_k_branches([b0], source_v, b, Jmax2,
                                           k, j, limit))
            q = scalar_filter(next_q, source_v, limit)
            next_q = []
            for b0 in q:
                next_q.extend(c_volume_branches([b0], source_v, c,
                                                 Jmax2, i, k, limit))
            q = next_q
        out.extend(q)
        _guard(out, limit, "sparse KKV disjunction")
    return project_final(out, (source_v,), limit)


def contains_key(branches, key):
    return any(contains_joint(branch, key) for branch in branches)


def branches_profile(branches):
    return {
        "branch_count": len(branches),
        "stored_entries": sum(profile(branch)["stored_entries"]
                               for branch in branches),
        "spin_sector_count_sum": sum(len(branch) for branch in branches),
        "joint_keys_enumerated": False,
        "amplitudes_evaluated": False,
        "branch_identity_preserved": True,
        "higher_order_correlations_preserved": False,
    }


def project_final(branches, charged_nodes, limit=100000):
    """Apply the requested Gauss/charged projection only after the full word."""
    out = []
    for branch in branches:
        # Projection is terminal, but expanding all omitted factors before
        # projection must not make the returned representation dense. Keep
        # the branch's inherited active scope; omitted pairs remain Cartesian.
        pair_scope = _pair_scope_of(branch)
        projected = _sparse_project(branch, charged_nodes, limit)
        if projected:
            out.append(_compress_pairs(projected, pair_scope))
    total = sum(
        sum(map(len, state["unary"])) +
        sum(map(len, state["pairs"].values()))
        for branch in out for state in branch.values()
    )
    if limit is not None and total > limit:
        raise RuntimeError(
            f"sparse-active-path final Gauss projection entries {total} "
            f"exceed guard {limit}; no truncation"
        )
    return out
