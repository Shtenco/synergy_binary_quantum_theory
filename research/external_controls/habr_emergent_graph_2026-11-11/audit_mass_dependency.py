#!/usr/bin/env python3
"""Dependency-aware audit of pairwise mass ratios.

Usage:
  python audit_mass_dependency.py path/to/mass_relation.json

The script checks the exact dependency structure of the directed ratio table.
It does not treat pair rows as independent observations and does not invent
uncertainties absent from the source JSON.
"""
from __future__ import annotations
import json
import sys
import numpy as np


def main(path: str) -> None:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    particles = data["particles"]
    relations = data["relations"]
    names = [p["name"] for p in particles]
    masses = np.asarray([float(p["mass"]) for p in particles], dtype=float)
    n = len(names)
    index = {name: i for i, name in enumerate(names)}

    # Incidence rows encode log(m_i/m_j) = log(m_i)-log(m_j).
    incidence = []
    exact_ratio_mismatches = []
    pair_rows = set()
    reciprocal_rows = 0
    for row in relations:
        a, b = row["pair"].split(" / ")
        if a not in index or b not in index:
            exact_ratio_mismatches.append((row["pair"], "unknown particle"))
            continue
        i, j = index[a], index[b]
        vec = np.zeros(n)
        vec[i], vec[j] = 1.0, -1.0
        incidence.append(vec)
        pair_rows.add((min(i, j), max(i, j)))
        observed = float(row["experimental_value"])
        reconstructed = masses[i] / masses[j]
        if not np.isclose(observed, reconstructed, rtol=1e-12, atol=0.0):
            exact_ratio_mismatches.append(
                (row["pair"], observed, reconstructed)
            )
    B = np.asarray(incidence, dtype=float)
    rank = int(np.linalg.matrix_rank(B)) if len(B) else 0

    # Directed reciprocal rows are exact sign reversals in log space.
    row_keys = {(tuple(np.flatnonzero(v)), tuple(v[v != 0])) for v in B}
    directed = len(relations)
    unordered = len(pair_rows)
    max_independent_log_ratios = max(n - 1, 0)

    # Graph Laplacian B^T B: for a complete pair table its rank is n-1.
    gram = B.T @ B if len(B) else np.zeros((n, n))
    gram_rank = int(np.linalg.matrix_rank(gram))

    print(f"particles={n}")
    print(f"directed_relation_rows={directed}")
    print(f"unordered_particle_pairs={unordered}")
    print(f"incidence_matrix_rank={rank}")
    print(f"ratio_gram_rank={gram_rank}")
    print(f"max_independent_log_ratios={max_independent_log_ratios}")
    print(f"redundant_directed_rows={directed-rank}")
    print(f"experimental_ratio_reconstruction_mismatches={len(exact_ratio_mismatches)}")
    print("uncertainty_covariance=ABSENT (no valid chi-square from this JSON alone)")
    if exact_ratio_mismatches:
        for item in exact_ratio_mismatches[:20]:
            print("MISMATCH", item)

    # Check teta convention against both directions.
    prose_mismatch = 0
    json_mismatch = 0
    for row in relations:
        model = float(row["model_value"])
        exp = float(row["experimental_value"])
        teta = float(row["teta"])
        if not np.isclose(teta, exp/model, rtol=1e-10, atol=1e-12):
            prose_mismatch += 1
        if not np.isclose(teta, model/exp, rtol=1e-10, atol=1e-12):
            json_mismatch += 1
    print(f"teta_not_matching_article_exp_over_model={prose_mismatch}")
    print(f"teta_not_matching_json_model_over_exp={json_mismatch}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python audit_mass_dependency.py path/to/mass_relation.json")
    main(sys.argv[1])
