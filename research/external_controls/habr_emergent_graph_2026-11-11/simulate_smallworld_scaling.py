#!/usr/bin/env python3
"""Finite-size control for the article's p_n = p/r^2, p=1/(K L) graph law.

Operationalization (explicit because the article does not fully specify one):
- periodic cubic L×L×L lattice;
- six nearest-neighbour local edges;
- for each nonlocal unordered vertex pair at toroidal Euclidean distance r,
  add an edge independently with probability p/r^2;
- p = 1/(K L^alpha), K=6;
- report sampled mean shortest-path length and combinatorial-Laplacian gap.

Dependencies: networkx, numpy, scipy.
This tests one natural implementation, not every possible interpretation.
"""
from __future__ import annotations
import argparse
import math
import random
import time
import networkx as nx
import numpy as np
from scipy.sparse.csgraph import laplacian as cs_laplacian
from scipy.sparse.linalg import eigsh


def make_graph(L: int, seed: int, alpha: float, K: float = 6.0):
    rng = random.Random(seed)
    coords = [(x, y, z) for x in range(L) for y in range(L) for z in range(L)]
    index = {c: i for i, c in enumerate(coords)}
    graph = nx.Graph()
    graph.add_nodes_from(range(len(coords)))

    for (x, y, z), i in index.items():
        for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            j = index[((x + dx) % L, (y + dy) % L, (z + dz) % L)]
            graph.add_edge(i, j)

    p = 1.0 / (K * L**alpha)
    for i, (x, y, z) in enumerate(coords):
        for j in range(i + 1, len(coords)):
            if graph.has_edge(i, j):
                continue
            a, b, c = coords[j]
            dx = min(abs(x-a), L-abs(x-a))
            dy = min(abs(y-b), L-abs(y-b))
            dz = min(abs(z-c), L-abs(z-c))
            r2 = dx*dx + dy*dy + dz*dz
            if r2 <= 1:
                continue
            if rng.random() < min(1.0, p/r2):
                graph.add_edge(i, j)
    return graph, p


def sampled_mean_path(graph, seed: int, n_sources: int = 24):
    rng = random.Random(seed + 100000)
    nodes = list(graph.nodes)
    sources = rng.sample(nodes, min(n_sources, len(nodes)))
    distances = []
    for source in sources:
        distances.extend(nx.single_source_shortest_path_length(graph, source).values())
    return float(np.mean(distances))


def spectral_gap(graph):
    adjacency = nx.to_scipy_sparse_array(
        graph, nodelist=sorted(graph.nodes), format="csr", dtype=float
    )
    lap = cs_laplacian(adjacency, normed=False)
    values = eigsh(lap, k=2, sigma=1e-8, which="LM",
                   return_eigenvectors=False, tol=1e-7)
    values = np.sort(values)
    return float(values[1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", type=int, default=[8, 12, 16])
    parser.add_argument("--alphas", nargs="+", type=float, default=[0.5, 1.0, 1.5, 2.0])
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22])
    args = parser.parse_args()

    print("alpha,L,N,seed,p,mean_degree,mean_path,gap,gap_lnN2,seconds")
    for alpha in args.alphas:
        for L in args.sizes:
            for seed in args.seeds:
                start = time.time()
                graph, p = make_graph(L, seed, alpha)
                path = sampled_mean_path(graph, seed)
                gap = spectral_gap(graph)
                n = L**3
                print(f"{alpha},{L},{n},{seed},{p:.12g},"
                      f"{2*graph.number_of_edges()/n:.9g},{path:.9g},"
                      f"{gap:.9g},{gap*math.log(n)**2:.9g},"
                      f"{time.time()-start:.3f}")


if __name__ == "__main__":
    main()
