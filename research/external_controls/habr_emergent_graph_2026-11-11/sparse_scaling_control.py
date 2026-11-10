#!/usr/bin/env python3
"""Larger finite-size control for p_r=p/r^2 on a periodic cubic lattice.

The random graph is Poissonized for efficient sparse sampling:
each endpoint emits Poisson(p/(2 r^2)) directed proposals; the undirected
edge exists if either endpoint proposes it. Thus P(edge)=1-exp(-p/r^2),
which differs from Bernoulli(p/r^2) only at second order.

Also computes the exact spectrum of the expected Laplacian by FFT. This is
NOT the same quantity as E[lambda_2(L_random)]; both are reported separately.

Dependencies: numpy, scipy.
"""
from __future__ import annotations
import argparse
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def kernel(L: int):
    d = np.minimum(np.arange(L), L - np.arange(L)).astype(float)
    dx, dy, dz = np.meshgrid(d, d, d, indexing="ij")
    r2 = dx*dx + dy*dy + dz*dz
    w = np.zeros((L, L, L), dtype=float)
    w[r2 > 1] = 1.0 / r2[r2 > 1]
    return w


def expected_operator(L: int, alpha: float = 1.0):
    """Exact eigenvalues of the expected combinatorial Laplacian."""
    w = kernel(L)
    p = 1.0 / (6.0 * L**alpha)
    nonlocal = p * (w.sum() - np.fft.fftn(w).real)
    k = np.arange(L)
    one_dim = 2.0 * (1.0 - np.cos(2.0*np.pi*k/L))
    local = (one_dim[:, None, None] + one_dim[None, :, None]
             + one_dim[None, None, :])
    eig = nonlocal + local
    eig[0, 0, 0] = np.inf
    return float(eig.min()), float(6.0 + p*w.sum()), float(w.sum())


def sample_graph(L: int, seed: int, alpha: float = 1.0):
    rng = np.random.default_rng(seed)
    n = L**3
    ids = np.arange(n).reshape((L, L, L))
    rows, cols = [], []

    # Six deterministic local neighbours, inserted as undirected edges.
    for axis in range(3):
        u = ids.ravel()
        v = np.roll(ids, -1, axis=axis).ravel()
        rows.extend(np.concatenate((u, v)))
        cols.extend(np.concatenate((v, u)))

    w = kernel(L)
    flat = w.ravel()
    S = float(flat.sum())
    probs = flat / S
    p = 1.0 / (6.0 * L**alpha)
    # Each endpoint proposes at half-rate. The union gives
    # P(edge)=1-exp(-p/r^2), close to p/r^2.
    counts = rng.poisson(p*S/2.0, size=n)
    src = np.repeat(np.arange(n), counts)
    choice = rng.choice(n, size=len(src), p=probs)

    ix = choice // (L*L)
    iy = (choice // L) % L
    iz = choice % L
    x = src // (L*L)
    y = (src // L) % L
    z = src % L
    dst = (((x+ix) % L)*L + ((y+iy) % L))*L + ((z+iz) % L)

    valid = src != dst
    a = np.minimum(src[valid], dst[valid])
    b = np.maximum(src[valid], dst[valid])
    # Encode and deduplicate unordered proposals.
    pairs = np.unique(a*n+b)
    a, b = pairs//n, pairs % n
    rows.extend(np.concatenate((a, b)))
    cols.extend(np.concatenate((b, a)))

    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A.data[:] = 1.0
    return A, len(pairs)


def sampled_gap(A):
    degree = np.asarray(A.sum(axis=1)).ravel()
    lap = sp.diags(degree) - A
    vals = eigsh(lap, k=2, which="SM", return_eigenvectors=False,
                 tol=1e-5, maxiter=8000)
    vals.sort()
    return float(vals[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", nargs="+", type=int, default=[16, 24, 32, 48])
    ap.add_argument("--seeds", nargs="+", type=int, default=[101, 202, 303])
    args = ap.parse_args()
    print("L,N,seed,mean_degree,sampled_gap,sampled_gap_lnN2,"
          "expected_operator_gap,expected_gap_lnN2,expected_mean_degree,"
          "expected_shortcut_degree")
    for L in args.sizes:
        exp_gap, exp_degree, shell_sum = expected_operator(L)
        p = 1.0/(6.0*L)
        for seed in args.seeds:
            A, shortcuts = sample_graph(L, seed)
            gap = sampled_gap(A)
            n = L**3
            mean_degree = A.nnz/n
            print(f"{L},{n},{seed},{mean_degree:.10g},{gap:.10g},"
                  f"{gap*math.log(n)**2:.10g},{exp_gap:.10g},"
                  f"{exp_gap*math.log(n)**2:.10g},{exp_degree:.10g},"
                  f"{p*shell_sum:.10g}")


if __name__ == "__main__":
    main()
