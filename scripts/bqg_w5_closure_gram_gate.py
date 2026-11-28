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
import math
from pathlib import Path


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
    payload = json.dumps(normalized, separators=(",", ":")).encode("utf-8")
    return {
        "status": "exact canonical K5 spin-network closure Gram certificate",
        "passed": True,
        "basis_dimension": n,
        "gram_rank": n,
        "gram_matrix": {"type": "identity", "dimension": n},
        "input_metric_Q": {"representation": "implicit_identity", "dimension": n},
        "basis_id_sha256": hashlib.sha256(payload).hexdigest(),
        "graph_convention": "fixed labelled K5 vertices 0..4 and ten labelled edges; no graph-automorphism quotient",
        "normalization": {
            "haar_measure": "normalized",
            "link_functions": "sqrt(2j+1) D^j_mn",
            "intertwiners": "normalized CG coupling with 1/sqrt(2K+1)",
            "orientation": "unitary epsilon_j magnetic-index maps",
        },
        "claim_scope": (
            "Exact metric of the supplied finite list of distinct admissible "
            "canonical K5 spin-network labels on the fixed labelled graph, without "
            "quotienting graph automorphisms. Does not prove closure completeness "
            "or exactness of floating-point operator amplitudes."
        ),
    }


def load_keys(path: Path) -> tuple[list, str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if "basis_keys_repr" in data:
        status = data.get("status")
        if status not in {"PASS_ADJOINT_CONSISTENCY", "FAIL_ADJOINT_OR_LEAKAGE_GATE"}:
            raise ValueError("basis_keys_repr requires a completed closure artifact status")
        raw_keys = data["basis_keys_repr"]
        if (not isinstance(raw_keys, list) or not raw_keys
                or any(not isinstance(x, str) for x in raw_keys)):
            raise ValueError("basis_keys_repr must be a nonempty list of repr strings")
        if data.get("closure_size") != len(raw_keys):
            raise ValueError("closure_size does not match basis_keys_repr length")
        if raw_keys != sorted(raw_keys):
            raise ValueError("completed basis keys are not in canonical repr order")
        basis_payload = "\n".join(raw_keys).encode("utf-8")
        expected_basis_id = hashlib.sha256(basis_payload).hexdigest()
        if data.get("basis_id_sha256") != expected_basis_id:
            raise ValueError("basis_id_sha256 does not match the ordered basis keys")

        def validate_sparse_columns(name, columns):
            if not isinstance(columns, list) or len(columns) != len(raw_keys):
                raise ValueError(f"completed artifact {name} column count does not match basis")
            for col_index, column in enumerate(columns):
                if not isinstance(column, list):
                    raise ValueError(f"{name}[{col_index}] must be a sparse-entry list")
                seen_rows = set()
                for entry in column:
                    if not isinstance(entry, list) or len(entry) != 3:
                        raise ValueError(f"{name}[{col_index}] entries must be [row, real, imag]")
                    row, real, imag = entry
                    if type(row) is not int or not 0 <= row < len(raw_keys):
                        raise ValueError(f"{name}[{col_index}] contains an out-of-basis row index")
                    if row in seen_rows:
                        raise ValueError(f"{name}[{col_index}] contains a duplicate row index")
                    seen_rows.add(row)
                    if (type(real) not in (int, float) or type(imag) not in (int, float)
                            or not math.isfinite(real) or not math.isfinite(imag)):
                        raise ValueError(f"{name}[{col_index}] contains a non-finite or non-numeric amplitude")

        validate_sparse_columns("forward_B_columns_sparse", data.get("forward_B_columns_sparse"))
        validate_sparse_columns("direct_Bdag_columns_sparse", data.get("direct_Bdag_columns_sparse"))
        return [ast.literal_eval(x) for x in raw_keys], "complete_basis_list"
    partial = data.get("partial_forward_columns_top_amplitudes")
    if partial is not None:
        if data.get("status") != "INCOMPLETE_CLOSURE_LIMIT_REACHED":
            raise ValueError("partial input columns require an explicit incomplete-closure status")
        if not isinstance(partial, list) or not partial:
            raise ValueError("partial input-column list must be nonempty")
        if data.get("columns_computed") != len(partial):
            raise ValueError("columns_computed does not match partial input-column count")
        return [ast.literal_eval(row["input_key"]) for row in partial], "computed_input_columns_only"
    raise ValueError("JSON contains neither a validated complete basis nor partial input-key list")


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
