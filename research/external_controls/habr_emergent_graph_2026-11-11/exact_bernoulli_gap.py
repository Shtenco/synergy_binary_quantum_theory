#!/usr/bin/env python3
"""Exact Bernoulli control for the cubic-torus + inverse-square shortcut graph.

For every unordered nonlocal pair {u,v}, independently include the edge with
q_uv = p / r(u,v)^2, p = 1/(6L), where r is periodic Euclidean distance.
The six nearest-neighbour torus edges are deterministic.

Sampling is exact: for each displacement class with M eligible unordered
pairs and common probability q, draw K~Binomial(M,q), then choose K of the M
pairs uniformly without replacement. This is distributionally identical to
M independent Bernoulli(q) trials; it is not Poissonization.

Usage:
  python exact_bernoulli_gap.py --sizes 8 10 12 14 16 20 24 28 32 --reps 20
  python exact_bernoulli_gap.py --sizes 40 --reps 10 --seed 3000
"""
from __future__ import annotations
import argparse
import csv
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components


def vertex_id(x: int, y: int, z: int, L: int) -> int:
    return (x % L) * L * L + (y % L) * L + (z % L)


def build_graph(L: int, seed: int, pscale: float = 1.0):
    if L < 4 or L % 2:
        raise ValueError("Use even L >= 4 for the stated torus displacement enumeration")
    rng = np.random.default_rng(seed)
    N = L ** 3
    edges: set[tuple[int, int]] = set()

    # Exactly 3N undirected nearest-neighbour edges on the periodic cubic lattice.
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = vertex_id(x, y, z, L)
                for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    j = vertex_id(x + dx, y + dy, z + dz, L)
                    edges.add((min(i, j), max(i, j)))

    p = pscale / (6.0 * L)
    nonlocal_added = 0

    # One representative from each pair of opposite displacement classes.
    for dx in range(L):
        for dy in range(L):
            for dz in range(L):
                d = (dx, dy, dz)
                inv = ((-dx) % L, (-dy) % L, (-dz) % L)
                if d == (0, 0, 0) or not d < inv:
                    continue

                md = (min(dx, L-dx), min(dy, L-dy), min(dz, L-dz))
                r2 = sum(v*v for v in md)
                if r2 == 1:  # nearest-neighbour edges are already deterministic
                    continue
                q = p / r2
                if q <= 0.0:
                    continue
                if q > 1.0:
                    raise ValueError(f"Invalid Bernoulli probability q={q}")

                self_inverse = d == inv
                if self_inverse:
                    starts = []
                    for i in range(N):
                        x, y, z = i // (L*L), (i // L) % L, i % L
                        j = vertex_id(x + dx, y + dy, z + dz, L)
                        if i < j:
                            starts.append(i)
                    starts = np.asarray(starts, dtype=np.int64)
                    M = len(starts)
                else:
                    starts = None
                    M = N

                # Exact independent Bernoulli(q) sample via Binomial + uniform subset.
                k = int(rng.binomial(M, q))
                if k == 0:
                    continue
                selected = rng.choice(M, size=k, replace=False)
                chosen_starts = starts[selected] if self_inverse else selected
                for i0 in chosen_starts:
                    i = int(i0)
                    x, y, z = i // (L*L), (i // L) % L, i % L
                    j = vertex_id(x + dx, y + dy, z + dz, L)
                    edge = (min(i, j), max(i, j))
                    if edge not in edges:
                        edges.add(edge)
                        nonlocal_added += 1

    e = np.asarray(list(edges), dtype=np.int64)
    rows = np.concatenate((e[:, 0], e[:, 1]))
    cols = np.concatenate((e[:, 1], e[:, 0]))
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N))
    return A, nonlocal_added


def measure(L: int, seed: int):
    A, n_shortcuts = build_graph(L, seed)
    degree = np.asarray(A.sum(axis=1)).ravel()
    lap = sp.diags(degree) - A
    vals = sla.eigsh(lap, k=2, which="SM", return_eigenvectors=False,
                     tol=1e-7, maxiter=100000)
    vals = np.sort(vals)
    ncomp = int(connected_components(A, directed=False, return_labels=False))
    return {
        "L": L, "N": L**3, "seed": seed, "lambda2": float(vals[1]),
        "mean_degree": float(degree.mean()), "shortcuts": n_shortcuts,
        "components": ncomp,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", nargs="+", type=int, default=[8,10,12,14,16,20,24,28,32])
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1000)
    ap.add_argument("--csv", default="exact_bernoulli_gap.csv")
    args = ap.parse_args()
    rows = []
    for L in args.sizes:
        start = time.time()
        batch = [measure(L, args.seed + 1000*L + r) for r in range(args.reps)]
        rows.extend(batch)
        x = np.asarray([r["lambda2"] for r in batch])
        mean, sd = float(x.mean()), float(x.std(ddof=1)) if len(x)>1 else (float(x.mean()), float("nan"))
        cv = sd/mean if mean else float("nan")
        print(f"L={L:3d} N={L**3:7d} reps={len(x):3d} "
              f"mean_gap={mean:.9f} sd={sd:.9f} CV={cv:.6%} "
              f"mean_degree={np.mean([r['mean_degree'] for r in batch]):.6f} "
              f"connected={sum(r['components']==1 for r in batch)}/{len(batch)} "
              f"elapsed={time.time()-start:.1f}s")
    with open(args.csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} graph realizations to {args.csv}")


if __name__ == "__main__":
    main()
