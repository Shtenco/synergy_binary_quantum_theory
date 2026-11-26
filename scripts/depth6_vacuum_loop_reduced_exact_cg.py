#!/usr/bin/env python3
"""Exact reduced-basis CG evaluator for the recovered K5 node constraint on vacuum.

This mirrors the conventions of bqg_block_apply_invariants:
endpoint-specific local_hit_outputs, oriented coupled tensors, reversed-edge
adjoints, path midpoint contractions, radial-adjoint/volume ordering, and
final Gauss projection. Charged intermediate sectors are retained until the
source-matched final projection.
"""
from __future__ import annotations

import itertools
import json
import argparse
from pathlib import Path
from functools import lru_cache

import sympy as sp
from sympy.physics.wigner import clebsch_gordan

VERT = tuple(range(5))
SOURCE_PROVENANCE = {
    "repository": "Shtenco/synergy_binary_quantum_theory",
    "commit": "f6bfcd775b8c3cd33c6c235d48293e40657baf9d",
    "bundle_path": "mixed_parallel/bqg_mixed_master_bundle.tar.gz.b64",
    "bundle_sha256": "2a42c6cdbf82b3b65a2ada02116e31cc9dd30c514589bb2885992050d1d16d32",
    "source_module_sha256": {
        "bqg_block_apply_invariants.py": "199a73d45645981f6a6c195bce8b595519172a71a0555004d14b5b7a4afb39f5",
        "bqg_depth6_generic_master_row.py": "0435668bb05595598a020fbb1b4dde471b6cb09b634fbe078a5976336f8b505e",
        "bqg_safe_boundary_master_gate.py": "9570b696e241ab1b01e5e327f86798ad6099f90b39d8db5cd5ceeac04a268c05",
        "bqg_reduced_transfer_hamiltonian.py": "7179b328431a31ffec5ce1e66ec94e97d02f78ec49323b9d4e87ab6cb881fe8f",
        "bqg_covariant_volume_recoupling.py": "45ab0a5befc51ad8323e3c9d658662318bcb2d58460387e8fe0b8b1f6e1afba4",
    },
    "comparison_scope": "endpoint local_hit_outputs, oriented coupled tensors, reversed-edge adjoints, path ordering, radial adjoint order, oriented simplex signs, and final local-J=0 Gauss projection",
    "identity_with_historical_branch_sum": "not established",
}
EDGES = tuple(itertools.combinations(VERT, 2))
EIDX = {e: i for i, e in enumerate(EDGES)}
NEIG = {v: tuple(w for w in VERT if w != v) for v in VERT}
LEGIDX = {(v, w): NEIG[v].index(w) for v in VERT for w in VERT if w != v}


def mvals(s):
    return tuple(range(s, -s - 1, -2))


@lru_cache(None)
def cg(s1, s2, S, m1, m2, M):
    return clebsch_gordan(*(sp.Rational(x, 2) for x in (s1, s2, S, m1, m2, M)))


def allowed_pair(s1, s2):
    return tuple(range(abs(s1 - s2), s1 + s2 + 1, 2))


@lru_cache(None)
def allowed_total_J(spins):
    s1, s2, s3, s4 = spins
    vals = set()
    for K in allowed_pair(s1, s2):
        for L in allowed_pair(s3, s4):
            vals.update(range(abs(K - L), K + L + 1, 2))
    return tuple(sorted(vals))


@lru_cache(None)
def multiplicity_labels(spins, J):
    s1, s2, s3, s4 = spins
    return tuple((K, L) for K in allowed_pair(s1, s2)
                 for L in allowed_pair(s3, s4)
                 if abs(K - L) <= J <= K + L and (J - abs(K - L)) % 2 == 0)


def apply_axis(tensor, axis, matrix):
    out = {}
    for index, val in tensor.items():
        for row in range(matrix.rows):
            z = matrix[row, index[axis]] * val
            if z == 0:
                continue
            key = index[:axis] + (row,) + index[axis + 1:]
            out[key] = out.get(key, 0) + z
    return {k: sp.simplify(v) for k, v in out.items() if v != 0}


@lru_cache(None)
def coupled_tensor(spins, K, L, J, M):
    spins = tuple(spins)
    dims = tuple(s + 1 for s in spins)
    mvs = [mvals(s) for s in spins]
    out = {}
    for index in itertools.product(*(range(d) for d in dims)):
        m1, m2, m3, m4 = [mvs[a][index[a]] for a in range(4)]
        MK, ML = m1 + m2, m3 + m4
        if MK + ML != M or abs(MK) > K or abs(ML) > L:
            continue
        value = (cg(spins[0], spins[1], K, m1, m2, MK)
                 * cg(spins[2], spins[3], L, m3, m4, ML)
                 * cg(K, L, J, MK, ML, M))
        if value != 0:
            out[index] = sp.simplify(value)
    return out


@lru_cache(None)
def oriented_coupled_tensor(v, spins, K, L, J, M):
    tensor = coupled_tensor(tuple(spins), K, L, J, M)
    for leg, w in enumerate(NEIG[v]):
        if w >= v or spins[leg] == 0:
            continue
        values = mvals(spins[leg])
        pos = {m: i for i, m in enumerate(values)}
        eps = sp.zeros(len(values))
        for m in values:
            eps[pos[m], pos[-m]] = (-1) ** ((spins[leg] - m) // 2)
        tensor = apply_axis(tensor, leg, eps)
    return tensor


@lru_cache(None)
def cg_leg_matrix(sin, sout, a2):
    rows, cols = mvals(sout), mvals(sin)
    mat = sp.zeros(len(rows), len(cols))
    for col, m in enumerate(cols):
        for row, mout in enumerate(rows):
            mat[row, col] = cg(sin, 1, sout, m, a2, mout)
    return mat


@lru_cache(None)
def local_hit_outputs(v, spins, leg, s_out, kind, color, label):
    spins = tuple(spins)
    K, L, J, M = label
    s_in = spins[leg]
    a2 = 1 if color == 0 else -1
    if kind == "cg":
        A = cg_leg_matrix(s_in, s_out, a2)
    else:
        A = cg_leg_matrix(s_out, s_in, a2).conjugate().T
    spo = list(spins)
    spo[leg] = s_out
    spo = tuple(spo)
    Bin = oriented_coupled_tensor(v, spins, K, L, J, M)
    Y = apply_axis(Bin, leg, A)
    outputs = []
    for Jo in allowed_total_J(spo):
        if abs(Jo - J) != 1:
            continue
        for Mo in mvals(Jo):
            for Ko, Lo in multiplicity_labels(spo, Jo):
                Bout = oriented_coupled_tensor(v, spo, Ko, Lo, Jo, Mo)
                amp = sp.simplify(sum(sp.conjugate(Bout.get(k, 0)) * y for k, y in Y.items()))
                if amp != 0:
                    outputs.append(((Ko, Lo, Jo, Mo), amp))
    return tuple(outputs)


def local_spins(spins, v):
    return tuple(spins[EIDX[tuple(sorted((v, w)))]] for w in NEIG[v])


def apply_hit(state, x, y, i, j):
    out = {}
    ei = EIDX[tuple(sorted((x, y)))]
    u, w = sorted((x, y))
    for (spins, locs), amp in state.items():
        s = spins[ei]
        for so in (s - 1, s + 1):
            if so < 0 or so > 1:
                continue
            if x < y:
                kind_u, col_u, kind_w, col_w = "cg", i, "cg", j
                norm = sp.sqrt(sp.Rational(s + 1, so + 1))
            else:
                kind_u, col_u, kind_w, col_w = "adj", j, "adj", i
                norm = sp.sqrt(sp.Rational(so + 1, s + 1))
            lsu, lsw = local_spins(spins, u), local_spins(spins, w)
            ou = local_hit_outputs(u, lsu, LEGIDX[(u, w)], so, kind_u, col_u, locs[u])
            ow = local_hit_outputs(w, lsw, LEGIDX[(w, u)], so, kind_w, col_w, locs[w])
            if not ou or not ow:
                continue
            new_spins = list(spins)
            new_spins[ei] = so
            new_spins = tuple(new_spins)
            for lu, au in ou:
                for lw, aw in ow:
                    new_locs = list(locs)
                    new_locs[u], new_locs[w] = lu, lw
                    key = (new_spins, tuple(new_locs))
                    out[key] = sp.simplify(out.get(key, 0) + amp * norm * au * aw)
    return {k: v for k, v in out.items() if v != 0}


def apply_path(path, iout, iin):
    initial_spins = (0,) * len(EDGES)
    initial_locs = tuple((0, 0, 0, 0) for _ in VERT)
    total = {}
    for mids in itertools.product(range(2), repeat=len(path) - 2):
        cols = (iout,) + tuple(mids) + (iin,)
        state = {(initial_spins, initial_locs): sp.Integer(1)}
        for r, (x, y) in enumerate(zip(path[:-1], path[1:])):
            state = apply_hit(state, x, y, cols[r], cols[r + 1])
            if not state:
                break
        for key, amp in state.items():
            total[key] = sp.simplify(total.get(key, 0) + amp)
    # Keep charged local J,M sectors: this is an intermediate path, before
    # volume/radial action and before the final H_v Gauss projection.
    return {k: sp.simplify(v) for k, v in total.items() if v != 0}


def merge_state(dst, src, scale=1):
    for key, value in src.items():
        z = sp.simplify(dst.get(key, 0) + scale * value)
        if z == 0:
            dst.pop(key, None)
        else:
            dst[key] = z
    return dst


def project_gauss_output(state):
    """Mirror the recovered H_node_apply final local-J=0 projection."""
    return {
        key: value for key, value in state.items()
        if all(J == 0 and M == 0 and K == L for K, L, J, M in key[1])
    }


def apply_volume_at_v(state, v):
    """Exact volume block for the only nonzero local sector reached here."""
    out = {}
    for (spins, locs), amp in state.items():
        K, L, J, M = locs[v]
        ls = local_spins(spins, v)
        if sum(s != 0 for s in ls) < 3:
            continue
        if tuple(sorted(ls)) != (0, 1, 1, 1):
            raise ArithmeticError(("unexpected local volume sector", v, ls, locs[v]))
        if J == 1:
            # Exact Q^2=(3/16)I in this 2D multiplicity space, hence V=3^(1/4)/2 I.
            z = sp.Pow(3, sp.Rational(1, 4)) / 2
        elif J == 3:
            # The multiplicity-one block has Q=0.
            z = sp.Integer(0)
        else:
            raise ArithmeticError(("unexpected total J in volume sector", J, locs[v]))
        if z != 0:
            key = (spins, locs)
            out[key] = sp.simplify(out.get(key, 0) + amp * z)
    return out


def radial_adjoint(state, v, c, i, j):
    """Exact reduced-basis radial-adjoint block in the recovered source order."""
    out = {}
    for k in range(2):
        first = apply_hit(state, c, v, k, j)
        # V P(v,c) P(c,v), applied right-to-left: two hits, then volume.
        plus = apply_hit(first, v, c, i, k)
        merge_state(out, apply_volume_at_v(plus, v), +1)
        # P(v,c) V P(c,v): first hit, volume, final hit.
        middle_volume = apply_volume_at_v(first, v)
        minus = apply_hit(middle_volume, v, c, i, k)
        merge_state(out, minus, -1)
    return out


def triangle_adjoint(i, j, v, a, b):
    pfd = apply_path((v, b, a, v), j, i)
    prd = apply_path((v, a, b, v), j, i)
    out = {}
    for key in set(pfd) | set(prd):
        value = sp.simplify(pfd.get(key, 0) - prd.get(key, 0))
        if value != 0:
            out[key] = value
    return out


def full_adjoint_term(i, j, v, a, b, c):
    loop = triangle_adjoint(i, j, v, a, b)
    return radial_adjoint(loop, v, c, i, j)


def oriented_specs(v):
    neighbors = NEIG[v]
    specs = []
    for r in range(4):
        tri = tuple(neighbors[i] for i in range(4) if i != r)
        sign = (-1) ** r
        a, b, c = tri
        specs.extend(((sign, (v, a, b, c)),
                      (sign, (v, b, c, a)),
                      (sign, (v, c, a, b))))
    return tuple(specs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--node", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    term_results = []
    node_total = {}
    for sign, (v, a, b, c) in oriented_specs(args.node):
        rows = [[full_adjoint_term(i, j, v, a, b, c) for j in range(2)] for i in range(2)]
        keys = set().union(*(set(rows[i][j]) for i in range(2) for j in range(2)))
        matrices = {str(key): sp.Matrix(2, 2, lambda i, j: rows[i][j].get(key, 0))
                    for key in keys}
        term_total = {}
        for i in range(2):
            for j in range(2):
                merge_state(term_total, rows[i][j], sp.Rational(1, 2))
        merge_state(node_total, term_total, sign)
        term_results.append({
            "spec": [v, a, b, c],
            "orientation_sign": sign,
            "adjoint_channel_matrices_by_output_state": {k: sp.sstr(m) for k, m in matrices.items()},
            "hermitian_half_sum_output": {str(k): sp.sstr(val) for k, val in term_total.items()},
            "term_annihilates_vacuum": not bool(term_total),
        })
    gauss_projected = project_gauss_output(node_total)
    result = {
        "node": args.node,
        "source_provenance": SOURCE_PROVENANCE,
        "oriented_term_count": len(term_results),
        "terms": term_results,
        "node_constraint_pre_gauss_output": {str(k): sp.sstr(v) for k, v in node_total.items()},
        "node_constraint_output": {str(k): sp.sstr(v) for k, v in gauss_projected.items()},
        "source_gauss_projection_applied": True,
        "node_constraint_annihilates_vacuum": not bool(gauss_projected),
        "exact_cg_arithmetic": True,
        "exact_volume_input": "For three doubled-spin-1 legs and one spin-zero leg, Q^2=(3/16)I on J=1 and Q=0 on J=3; hence V=3^(1/4)/2 I on J=1.",
        "forward_radial_piece": "zero analytically: volume acts on vacuum or after only one fundamental hit, leaving fewer than three nonzero incident legs",
        "scope": "exact reduced-basis forward-plus-adjoint H_node across all 12 oriented terms, with source-matched final Gauss projection; forward radial part is zero analytically on vacuum; source node formula assembled directly",
    }
    text = json.dumps(result, indent=2, ensure_ascii=False)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
