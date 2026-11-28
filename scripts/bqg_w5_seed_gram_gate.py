#!/usr/bin/env python3
"""Exact Gram-matrix gate for the five W5 seed intertwiners.

The seed has ten links with doubled spin s=4 (j=2), and differs only by
the recoupling channel K at node 0: K2 in {0,2,4,6,8}.  Each four-valent
intertwiner is the normalized singlet coupling of pair (0,1) and pair (2,3).

This gate proves orthonormality of that five-dimensional intertwiner sector
using exact SymPy Clebsch-Gordan coefficients.  The full spin-network seed
Gram matrix is the same identity under the standard normalized Peter-Weyl
link basis and normalized Haar inner product.  This script does not compute
the Gram matrix of a larger dynamically generated closure basis.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp
from sympy.physics.wigner import clebsch_gordan

DOUBLED_LINK_SPIN = 4
SEED_K2_VALUES = (0, 2, 4, 6, 8)


def intertwiner_vector(k2: int) -> dict[tuple[int, int, int, int], sp.Expr]:
    """Return normalized four-valent j=2 intertwiner in magnetic basis."""
    if k2 not in SEED_K2_VALUES:
        raise ValueError(f"unsupported seed K2={k2}")
    s = DOUBLED_LINK_SPIN
    magnetic = tuple(range(s, -s - 1, -2))
    j = sp.Rational(s, 2)
    K = sp.Rational(k2, 2)
    out = {}
    for m1, m2, m3, m4 in itertools.product(magnetic, repeat=4):
        value = sp.S.Zero
        for M2 in range(k2, -k2 - 1, -2):
            M = sp.Rational(M2, 2)
            phase = (-1) ** ((k2 - M2) // 2)
            c_left = clebsch_gordan(
                j, j, K,
                sp.Rational(m1, 2), sp.Rational(m2, 2), M,
            )
            c_right = clebsch_gordan(
                j, j, K,
                sp.Rational(m3, 2), sp.Rational(m4, 2), -M,
            )
            value += phase * c_left * c_right / sp.sqrt(k2 + 1)
        value = sp.simplify(value)
        if value != 0:
            out[(m1, m2, m3, m4)] = value
    return out


def exact_seed_gram() -> sp.Matrix:
    vectors = [intertwiner_vector(k2) for k2 in SEED_K2_VALUES]
    gram = sp.zeros(len(vectors), len(vectors))
    for i, left in enumerate(vectors):
        for j, right in enumerate(vectors):
            keys = set(left) | set(right)
            gram[i, j] = sp.simplify(sum(
                sp.conjugate(left.get(key, 0)) * right.get(key, 0)
                for key in keys
            ))
    return gram


def exact_full_seed_gram() -> sp.Matrix:
    """Assemble the full K5 spin-network Gram from normalized node/link factors.

    Node 0 carries the five-channel intertwiner family. Nodes 1..4 carry the
    same normalized K=0 intertwiner in every seed. Every one of the ten links
    carries the same fixed j=2 Peter-Weyl label. With normalized Haar measure
    and normalized Peter-Weyl matrix elements, edge contractions identify
    matching endpoint magnetic indices unitarily; the resulting spin-network
    inner product factorizes into the node intertwiner inner products.
    """
    node0 = exact_seed_gram()
    fixed_node_gram = sp.eye(1)
    normalized_edge_contraction = sp.eye(1)
    full = node0
    for _node in range(4):
        full = sp.kronecker_product(full, fixed_node_gram)
    for _edge in range(10):
        full = sp.kronecker_product(full, normalized_edge_contraction)
    return sp.simplify(full)


def run_gate() -> dict:
    node0_gram = exact_seed_gram()
    gram = exact_full_seed_gram()
    passed = (
        node0_gram == sp.eye(len(SEED_K2_VALUES))
        and gram == sp.eye(len(SEED_K2_VALUES))
    )
    return {
        "status": "exact W5 five-seed intertwiner Gram gate",
        "passed": bool(passed),
        "doubled_link_spin": DOUBLED_LINK_SPIN,
        "link_spin": "2",
        "seed_K2_values": list(SEED_K2_VALUES),
        "seed_K_values": [str(sp.Rational(k, 2)) for k in SEED_K2_VALUES],
        "node0_intertwiner_gram": [[str(node0_gram[i, j]) for j in range(node0_gram.cols)]
                                   for i in range(node0_gram.rows)],
        "factorization": {
            "node_0": "five-channel normalized intertwiner family",
            "nodes_1_to_4": "same normalized K=0 intertwiner in every seed",
            "links": "ten fixed j=2 normalized Peter-Weyl matrix-element factors",
            "haar_measure": "normalized Haar measure",
            "edge_contraction": "unitary magnetic-index identification; inner product factorizes",
        },
        "input_metric_Q": [[str(gram[i, j]) for j in range(gram.cols)]
                           for i in range(gram.rows)],
        "gram_rank": gram.rank(),
        "claim_scope": (
            "Exact orthonormality of the five seed intertwiners. Extending Q=I "
            "to a larger closure basis additionally relies on the standard "
            "normalized Peter-Weyl link basis and must be checked for the "
            "actual generated basis, including any duplicate or noncanonical states."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_gate()
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
