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


def _scope_of(support, extra=()):
    scope = set(extra)
    for state in support.values():
        for a, b in state["pairs"]:
            scope.update((a, b))
    return tuple(sorted(scope))


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
    # Retain inherited scopes only where a pair relation is genuinely tighter
    # than the current unary Cartesian product. Singleton seed relations are
    # redundant and must not force the initial scope to all five nodes.
    for state in support.values():
        for (a, b), relation in state["pairs"].items():
            cartesian = frozenset((x, y) for x in state["unary"][a]
                                  for y in state["unary"][b])
            if relation != cartesian:
                scope.update((a, b))
    work = compress(PG._lift(expand(support)), tuple(sorted(scope)))
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
    projected = PG._project(expand(work), charged_nodes, limit)
    return compress(projected, tuple(sorted(scope)))


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
                    # Keep only the primitive word's active pair scopes.
                    scope = active_scope(candidate, source_v)
                    part = compress(part, scope)
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
            out.append(compress(vh, _scope_of(vh, (source_v,))))
    vin = []
    for branch in branches:
        lifted = PG._lift(expand(branch))
        vol = PG._volume(lifted, source_v, limit)
        projected = PG._project(vol, NODES, limit)
        if projected:
            vin.append(compress(projected, _scope_of(projected, (source_v,))))
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
            inverse.append(compress(part, _scope_of(part, (source_v, target_v))))
    middle = k_branches(inverse, source_v, Jmax2,
                        NODES, limit)
    for branch in middle:
        closed = sequence(branch, (("P", (source_v, target_v)),),
                          source_v, Jmax2, NODES, limit)
        if closed:
            out.append(compress(closed, _scope_of(closed, (source_v, target_v))))
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
                out.append(compress(part, _scope_of(part, (source_v,))))
    inverse = []
    for branch in branches:
        part = sequence(branch, (("P", (target_v, source_v)),),
                        source_v, Jmax2, NODES, limit)
        if part:
            inverse.append(compress(part, (source_v, target_v)))
    middle = []
    for branch in inverse:
        part = sequence(branch, (("V", source_v),), source_v,
                         Jmax2, NODES, limit)
        if part:
            middle.append(compress(part, _scope_of(part, (source_v, target_v))))
    for branch in middle:
        part = sequence(branch, (("P", (source_v, target_v)),),
                        source_v, Jmax2, NODES, limit)
        if part:
            out.append(compress(part, _scope_of(part, (source_v, target_v))))
    return out


def scalar_filter(branches, source_v, limit=100000):
    out = []
    for branch in branches:
        full = PG._scalar_filter(expand(branch), source_v)
        if full:
            out.append(compress(full, _scope_of(branch, (source_v,))))
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
        projected = PG._project(expand(branch), charged_nodes, limit)
        if projected:
            out.append(projected)
    _guard(out, limit, "final Gauss projection")
    return out
