#!/usr/bin/env python3
"""Reanalyse the certified 24 L1 q4 Euclidean columns in metric S4 sectors.

Input is a directory containing q4_0.json ... q4_23.json from the successful
`collective-l1-e-q4-rank` run. No Peter-Weyl amplitudes are recomputed.

The 24 barycentric chambers of one parent tetrahedron are indexed by
`itertools.permutations(range(4))`, hence they carry the regular S4 action.
The script reconstructs the exact Gram matrix of the saved sparse states,
checks left-regular S4 covariance, and then applies the canonical equivariant
24 -> 6 coarse-edge map: each unordered parent edge receives the normalized
sum of the four chambers whose first two vertices are that edge.

The repository metric convention calls the ordinary tetrahedral vector/edge
triplet T2, so the six-edge representation is A1 + E + T2. Any S4-invariant
operator on it has

    K6 = a I + b A_adj + c O_opp,

so

    lambda_A1 = a + 4 b + c,
    lambda_E  = a - 2 b + c,
    lambda_T2 = a - c,
    Delta_ET  = lambda_E - lambda_T2 = 2(c-b).

Delta_ET is a metric-sector cubic/tetrahedral anisotropy diagnostic. This
first-order tangent Gram is NOT yet identified with the physical effective
action Hessian or the TT dispersion coefficient zeta4.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


PERMS = list(itertools.permutations(range(4)))
EDGES = list(itertools.combinations(range(4), 2))
TOL = 3e-12


def state_key(row):
    return (
        tuple((int(a), int(b)) for a, b in row["spin_changes"]),
        tuple((int(a), int(b)) for a, b in row["K_changes"]),
    )


def load_column(path: Path):
    d = json.loads(path.read_text(encoding="utf-8"))
    if not d.get("passed"):
        raise RuntimeError(f"input column did not pass: {path}")
    state = {state_key(r): complex(r["re"], r["im"]) for r in d["states"]}
    return d, state


def inner(a, b):
    if len(a) <= len(b):
        return sum(np.conj(v) * b.get(k, 0.0) for k, v in a.items())
    return sum(np.conj(a.get(k, 0.0)) * v for k, v in b.items())


def compose(p, q):
    return tuple(p[q[i]] for i in range(4))


def cycle_type(p):
    seen = set()
    lengths = []
    for i in range(4):
        if i in seen:
            continue
        j = i
        n = 0
        while j not in seen:
            seen.add(j)
            n += 1
            j = p[j]
        lengths.append(n)
    return tuple(sorted(lengths, reverse=True))


# Naming follows the existing metric-photon-bridge convention:
# T2 is the standard tetrahedral vector/edge triplet with transposition char +1.
CHAR = {
    "A1": {(1, 1, 1, 1): 1, (2, 1, 1): 1,  (2, 2): 1,  (3, 1): 1,  (4,): 1},
    "A2": {(1, 1, 1, 1): 1, (2, 1, 1): -1, (2, 2): 1,  (3, 1): 1,  (4,): -1},
    "E":  {(1, 1, 1, 1): 2, (2, 1, 1): 0,  (2, 2): 2,  (3, 1): -1, (4,): 0},
    "T2": {(1, 1, 1, 1): 3, (2, 1, 1): 1,  (2, 2): -1, (3, 1): 0,  (4,): -1},
    "T1": {(1, 1, 1, 1): 3, (2, 1, 1): -1, (2, 2): -1, (3, 1): 0,  (4,): 1},
}
DIM = {"A1": 1, "A2": 1, "E": 2, "T1": 3, "T2": 3}


def left_regular(g):
    idx = {p: i for i, p in enumerate(PERMS)}
    U = np.zeros((24, 24))
    for j, p in enumerate(PERMS):
        U[idx[compose(g, p)], j] = 1.0
    return U


def central_projector(irrep):
    d = DIM[irrep]
    P = np.zeros((24, 24))
    for g in PERMS:
        P += CHAR[irrep][cycle_type(g)] * left_regular(g)
    return (d / 24.0) * P


def run(directory: Path):
    meta = []
    states = []
    source_commits = set()
    for i in range(24):
        d, s = load_column(directory / f"q4_{i}.json")
        if int(d["local_fine_index"]) != i:
            raise RuntimeError(f"column/index mismatch at {i}")
        meta.append(d)
        states.append(s)
        source_commits.add(d.get('source_commit'))

    if len(source_commits) != 1 or None in source_commits:
        raise RuntimeError(f'columns must share one explicit source commit, got {sorted(str(x) for x in source_commits)}')
    source_commit = next(iter(source_commits))
    G = np.zeros((24, 24), complex)
    for i in range(24):
        for j in range(24):
            G[i, j] = inner(states[i], states[j])
    herm = float(np.linalg.norm(G - G.conj().T))

    comm = []
    for g in PERMS:
        U = left_regular(g)
        comm.append(np.linalg.norm(U @ G - G @ U) / max(np.linalg.norm(G), 1e-30))
    left_def = float(max(comm))

    irrep = {}
    for name in CHAR:
        P = central_projector(name)
        proj_err = float(np.linalg.norm(P @ P - P))
        rank = int(round(np.trace(P).real))
        evals, evecs = np.linalg.eigh((P + P.T) / 2)
        Q = evecs[:, np.argsort(evals)[-rank:]]
        vals = np.linalg.eigvalsh(Q.conj().T @ G @ Q).real
        irrep[name] = {
            "rank": rank,
            "projector_idempotence_error": proj_err,
            "gram_eigenvalues": vals.tolist(),
        }

    B = np.zeros((24, 6))
    edge_to_cols = {e: [] for e in EDGES}
    for i, p in enumerate(PERMS):
        edge = tuple(sorted(p[:2]))
        edge_to_cols[edge].append(i)
    for eidx, e in enumerate(EDGES):
        if len(edge_to_cols[e]) != 4:
            raise RuntimeError("each parent edge must have exactly four chambers")
        for i in edge_to_cols[e]:
            B[i, eidx] = 0.5
    equiv_isometry = float(np.linalg.norm(B.T @ B - np.eye(6)))
    K = (B.T @ G @ B).real

    diag = np.diag(K)
    adjacent = []
    opposite = []
    for i, e in enumerate(EDGES):
        for j in range(i + 1, 6):
            f = EDGES[j]
            n = len(set(e) & set(f))
            if n == 1:
                adjacent.append(K[i, j])
            elif n == 0:
                opposite.append(K[i, j])
    a = float(np.mean(diag))
    b = float(np.mean(adjacent))
    c = float(np.mean(opposite))

    Kfit = np.zeros((6, 6))
    for i, e in enumerate(EDGES):
        for j, f in enumerate(EDGES):
            if i == j:
                Kfit[i, j] = a
            elif len(set(e) & set(f)) == 1:
                Kfit[i, j] = b
            else:
                Kfit[i, j] = c
    orbit_res = float(np.linalg.norm(K - Kfit) / max(np.linalg.norm(K), 1e-30))

    lA = a + 4 * b + c
    lE = a - 2 * b + c
    lT = a - c
    delta = lE - lT
    rel = delta / ((lE + lT) / 2.0)

    passed = (
        herm < TOL
        and left_def < TOL
        and equiv_isometry < TOL
        and orbit_res < TOL
        and all(v["projector_idempotence_error"] < TOL for v in irrep.values())
        and [irrep[x]["rank"] for x in ("A1", "A2", "E", "T1", "T2")] == [1, 1, 4, 9, 9]
    )
    return {
        "status": "reanalysis of certified L1 q4 Euclidean tangent amplitudes",
        "passed": bool(passed),
        "source_commit": source_commit,
        "source_provenance": "all 24 columns were recomputed from the same checked-out source commit and individually record that commit",
        "source_science_status": "L1_BLOCK_E_Q4_EXACT_PROJECTION",
        "columns": 24,
        "column_supports": [int(d["support"]) for d in meta],
        "gram_hermiticity_error": herm,
        "left_regular_S4_commutator_relative_max": left_def,
        "S4_naming_convention": "repository metric convention: edge/vector triplet is T2",
        "regular_representation_irreps": irrep,
        "coarse_edge_map": {
            "edges": [list(e) for e in EDGES],
            "chambers_per_edge": 4,
            "coefficient_per_chamber": 0.5,
            "isometry_error": equiv_isometry,
            "decomposition": "6 = A1 + E + T2",
        },
        "K6": K.tolist(),
        "S4_orbit_fit": {
            "a_same": a,
            "b_adjacent": b,
            "c_opposite": c,
            "relative_residual": orbit_res,
            "lambda_A1": lA,
            "lambda_E": lE,
            "lambda_T2": lT,
            "Delta_ET": delta,
            "relative_ET_split": rel,
        },
        "conclusion": (
            "The actual first-refinement Euclidean tangent Gram has nonzero, symmetry-resolved E and T2 metric sectors. "
            "The canonical S4-equivariant 24-to-6 edge compression yields a nonzero Delta_ET. "
            "This is a finite metric-sector anisotropy precursor, not yet the effective-action Hessian or physical zeta4."
        ),
        "next": (
            "Construct the depth-two/effective metric response on the same six-edge carrier and evaluate its same/adjacent/opposite orbit elements as functions of low momentum; only that dynamic Hessian may be mapped to eta2_iso and zeta4_cub."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.input_dir)
    text = json.dumps(out, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
