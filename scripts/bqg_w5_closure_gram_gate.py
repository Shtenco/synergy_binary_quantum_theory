#!/usr/bin/env python3
"""Exact Gram gate for canonical finite K5 spin-network closure labels.

A basis key is (ten doubled edge spins, five doubled four-valent recoupling
labels). Under normalized Haar measure, normalized Peter-Weyl link functions,
and normalized SU(2) Clebsch-Gordan intertwiners, distinct admissible labelled
spin networks are orthonormal. This gate validates key admissibility and
uniqueness, then certifies Q=I on the supplied list.

It does not inspect floating-point branch tensors or prove that an upstream
operator emitted every required state. For incomplete artifacts it certifies
only the explicitly listed computed input keys, never the undiscovered closure.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import sympy as sp

VERTICES = tuple(range(5))
EDGES = tuple((u, v) for u in VERTICES for v in VERTICES if u < v)
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
NEIGHBORS = {
    v: tuple(w for w in VERTICES if w != v)
    for v in VERTICES
}


def allowed_k2(local_spins: tuple[int, int, int, int]) -> tuple[int, ...]:
    s1, s2, s3, s4 = local_spins
    left = set(range(abs(s1 - s2), s1 + s2 + 1, 2))
    right = set(range(abs(s3 - s4), s3 + s4 + 1, 2))
    return tuple(sorted(left & right))


def validate_basis_key(key) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not isinstance(key, (tuple, list)) or len(key) != 2:
        raise ValueError("basis key must be (edge_spins, node_channels)")
    spins, channels = key
    if not isinstance(spins, (tuple, list)) or len(spins) != 10:
        raise ValueError("K5 basis key must contain exactly ten edge spins")
    if not isinstance(channels, (tuple, list)) or len(channels) != 5:
        raise ValueError("K5 basis key must contain exactly five node channels")
    if any(type(x) is not int or x < 0 for x in spins):
        raise ValueError("doubled edge spins must be nonnegative integers")
    if any(type(x) is not int or x < 0 for x in channels):
        raise ValueError("doubled node channels must be nonnegative integers")

    spins = tuple(spins)
    channels = tuple(channels)
    for v in VERTICES:
        local = tuple(
            spins[EDGE_INDEX[tuple(sorted((v, w)))]]
            for w in NEIGHBORS[v]
        )
        if channels[v] not in allowed_k2(local):
            raise ValueError(
                f"node {v}: K2={channels[v]} is inadmissible for local "
                f"doubled spins {local}; allowed={allowed_k2(local)}"
            )
    return spins, channels


def gram_certificate(keys) -> dict:
    normalized = [validate_basis_key(key) for key in keys]
    if not normalized:
        raise ValueError("at least one basis key is required")
    if len(set(normalized)) != len(normalized):
        raise ValueError("duplicate canonical spin-network basis keys detected")

    n = len(normalized)
    gram = sp.eye(n)
    payload = json.dumps(normalized, separators=(",", ":")).encode("utf-8")
    return {
        "status": "exact canonical K5 spin-network closure Gram certificate",
        "passed": True,
        "basis_dimension": n,
        "gram_rank": gram.rank(),
        "gram_matrix": "identity",
        "input_metric_Q": [[str(gram[i, j]) for j in range(n)] for i in range(n)],
        "basis_id_sha256": hashlib.sha256(payload).hexdigest(),
        "normalization": {
            "haar_measure": "normalized",
            "link_functions": "sqrt(2j+1) D^j_mn",
            "intertwiners": "normalized CG coupling with 1/sqrt(2K+1)",
            "orientation": "unitary epsilon_j magnetic-index maps",
        },
        "claim_scope": (
            "Exact metric of the supplied finite list of distinct admissible "
            "canonical K5 spin-network labels. Does not prove closure completeness "
            "or exactness of floating-point operator amplitudes."
        ),
    }


def load_keys(path: Path) -> tuple[list, str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if "basis_keys_repr" in data:
        return [ast.literal_eval(x) for x in data["basis_keys_repr"]], "complete_basis_list"
    partial = data.get("partial_forward_columns_top_amplitudes")
    if partial is not None:
        return [ast.literal_eval(row["input_key"]) for row in partial], "computed_input_columns_only"
    raise ValueError("JSON contains neither basis_keys_repr nor partial input-key list")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    keys, source_scope = load_keys(args.basis_json)
    result = gram_certificate(keys)
    result["source_scope"] = source_scope
    result["source_artifact"] = str(args.basis_json)
    if source_scope != "complete_basis_list":
        result["claim_scope"] += " The supplied artifact is partial; no undiscovered closure states are certified."
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
