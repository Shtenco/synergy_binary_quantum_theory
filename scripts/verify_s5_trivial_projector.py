#!/usr/bin/env python3
"""Exact finite-group check of the recovered S3/J4/J5 trivial-sector selector.

All arithmetic is integer set arithmetic on permutations; no floating point,
random vectors, or external packages are used.
"""
from __future__ import annotations

from itertools import permutations

Perm = tuple[int, ...]
N = 5
IDENTITY: Perm = tuple(range(N))


def compose(left: Perm, right: Perm) -> Perm:
    """Composition left ∘ right in the same convention as the source tuples."""
    return tuple(left[right[i]] for i in range(N))


def transposition(i: int, j: int) -> Perm:
    p = list(IDENTITY)
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def symmetric_group(moved: tuple[int, ...]) -> set[Perm]:
    out: set[Perm] = set()
    for image in permutations(moved):
        p = list(IDENTITY)
        for source, target in zip(moved, image):
            p[source] = target
        out.add(tuple(p))
    return out


def left_coset_union(
    subgroup: set[Perm], representatives: tuple[Perm, ...]
) -> tuple[set[Perm], list[set[Perm]]]:
    cosets = [{compose(rep, member) for member in subgroup} for rep in representatives]
    union: set[Perm] = set()
    for coset in cosets:
        if union.intersection(coset):
            raise AssertionError("coset representatives are not disjoint")
        union.update(coset)
    return union, cosets


def verify() -> dict[str, object]:
    # Recovered source: S3 permutes {2,3,4}; S4 fixes 0 and permutes {1,2,3,4}.
    s3 = symmetric_group((2, 3, 4))
    s4 = symmetric_group((1, 2, 3, 4))
    s5 = symmetric_group((0, 1, 2, 3, 4))

    j4_reps = tuple(transposition(1, k) for k in (2, 3, 4))
    j5_reps = tuple(transposition(0, k) for k in (1, 2, 3, 4))

    s4_cosets, _ = left_coset_union(s3, (IDENTITY, *j4_reps))
    s5_cosets, _ = left_coset_union(s4, (IDENTITY, *j5_reps))

    assert (len(s3), len(s4), len(s5)) == (6, 24, 120)
    assert s4_cosets == s4
    assert s5_cosets == s5
    assert set(j4_reps) == {transposition(1, k) for k in (2, 3, 4)}
    assert set(j5_reps) == {transposition(0, k) for k in (1, 2, 3, 4)}

    # The disjoint coset identities give, on S3-invariant vectors,
    # (I+J4)/4 = e_S4 and (I+J5)/5 = e_S5 after the first projection.
    return {
        "valid": True,
        "arithmetic": "exact integer permutation-set identities",
        "group_orders": {"S3": len(s3), "S4": len(s4), "S5": len(s5)},
        "coset_partition_sizes": {"S4_from_S3": 4, "S5_from_S4": 5},
        "projector_identity": "((I+J5)/5)((I+J4)/4)e_S3 = e_S5",
        "selector_indices_zero_based": True,
        "scope": "exact [5] sector selector only; no master-map rank claim",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(verify(), indent=2))

