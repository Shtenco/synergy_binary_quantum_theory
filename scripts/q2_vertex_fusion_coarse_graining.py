#!/usr/bin/env python3
"""Exact vertex-block fusion and boundary-QFI equivalence for finite BQG graphs."""
from __future__ import annotations
import numpy as np
from scripts.q2_relational_higher_spin_graph_qfi_gate import graph_fixtures, invariant_singlet_tensor, random_intertwiner, spin_generators
from scripts.q2_relational_higher_spin_network_qfi_gate import four_valent_intertwiner_basis


def _partition(n, blocks):
    raw = tuple(tuple(b) for b in blocks)
    if any(not b for b in raw):
        raise ValueError("vertex blocks must be non-empty")
    flat = [v for b in raw for v in b]
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) for v in flat):
        raise ValueError("vertex blocks must contain integer indices")
    blocks = tuple(tuple(int(v) for v in b) for b in raw)
    flat = [v for b in blocks for v in b]
    if sorted(flat) != list(range(n)) or len(flat) != len(set(flat)):
        raise ValueError("vertex blocks must partition every vertex exactly once")
    return blocks


def _contract(tensors, edges, outputs):
    labels = [[None] * t.ndim for t in tensors]
    operands, nxt = [], 0
    for (ta, aa), (tb, ab), _s in edges:
        if labels[ta][aa] is not None or labels[tb][ab] is not None:
            raise ValueError("a tensor leg is used more than once")
        labels[ta][aa], labels[tb][ab] = nxt, nxt + 1
        nxt += 2
    out_labels = []
    for ti, axis in outputs:
        if labels[ti][axis] is not None:
            raise ValueError("output leg is contracted")
        labels[ti][axis] = nxt
        out_labels.append(nxt)
        nxt += 1
    for tensor, labs in zip(tensors, labels):
        if any(x is None for x in labs):
            raise ValueError("unclassified tensor axis")
        operands.extend([tensor, labs])
    for (ta, aa), (tb, ab), singlet in edges:
        operands.extend([singlet, [labels[ta][aa], labels[tb][ab]]])
    operands.append(out_labels)
    return np.einsum(*operands, optimize=True)


def _apply_one_leg(tensor, axis, generator):
    moved = np.moveaxis(tensor, axis, 0)
    transformed = np.tensordot(generator, moved, axes=(1, 0))
    return np.moveaxis(transformed, 0, axis)


def fuse_vertex_blocks(node_states, edges, singlet, blocks):
    """Contract all intra-block singlets, preserving every open/cut half-link."""
    n = len(node_states)
    blocks = _partition(n, blocks)
    owner = {v: bi for bi, block in enumerate(blocks) for v in block}
    used, normalized = set(), []
    for ei, edge in enumerate(edges):
        if len(edge) != 4:
            raise ValueError("edges must be (u, leg_u, v, leg_v)")
        u, lu, v, lv = edge
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer))
               for x in (u, lu, v, lv)):
            raise ValueError("edge indices must be integers")
        if u == v or not (0 <= u < n and 0 <= v < n and 0 <= lu < 4 and 0 <= lv < 4):
            raise ValueError("invalid graph edge")
        if (u, lu) in used or (v, lv) in used:
            raise ValueError("a half-link is used more than once")
        used.update(((u, lu), (v, lv)))
        normalized.append((int(u), int(lu), int(v), int(lv), ei))
    fused, boundaries = [], []
    d = singlet.shape[0]
    for bi, block in enumerate(blocks):
        local = sorted(block)
        local_idx = {v: i for i, v in enumerate(local)}
        tensors = [np.asarray(node_states[v], complex).reshape((d,) * 4) for v in local]
        internal = [((local_idx[u], lu), (local_idx[v], lv), singlet)
                    for u, lu, v, lv, _ in normalized if owner[u] == bi == owner[v]]
        cut = {(local[ta], aa) for (ta, aa), (tb, ab), _ in internal}
        cut |= {(local[tb], ab) for (ta, aa), (tb, ab), _ in internal}
        boundary = [(v, leg) for v in local for leg in range(4) if (v, leg) not in cut]
        fused.append(_contract(tensors, internal,
                               [(local_idx[v], leg) for v, leg in boundary]))
        boundaries.append(boundary)
    lookup = {(v, leg): (bi, axis) for bi, b in enumerate(boundaries)
              for axis, (v, leg) in enumerate(b)}
    coarse_edges = [(lookup[(u, lu)], lookup[(v, lv)], ei, (u, lu, v, lv))
                    for u, lu, v, lv, ei in normalized if owner[u] != owner[v]]
    return fused, boundaries, coarse_edges


def _micro_state(states, edges, singlet, twist=None):
    d = singlet.shape[0]
    tensors = [np.asarray(s, complex).reshape((d,) * 4).copy() for s in states]
    if twist is not None:
        ei, gen = twist
        u, lu, _v, _lv = edges[ei]
        tensors[u] = -1j * _apply_one_leg(tensors[u], lu, gen)
    specs = [((u, lu), (v, lv), singlet) for u, lu, v, lv in edges]
    occupied = {(u, lu) for u, lu, v, lv in edges} | {(v, lv) for u, lu, v, lv in edges}
    outputs = [(v, a) for v in range(len(states)) for a in range(4) if (v, a) not in occupied]
    return _contract(tensors, specs, outputs).reshape(-1)


def _coarse_state(fused, boundaries, coarse_edges, singlet, twist=None):
    tensors = [t.copy() for t in fused]
    if twist is not None:
        ei, gen = twist
        (bi, axis), _other, _idx, _original = next(e for e in coarse_edges if e[2] == ei)
        tensors[bi] = -1j * _apply_one_leg(tensors[bi], axis, gen)
    specs = [(a, b, singlet) for a, b, _i, _o in coarse_edges]
    occupied = {a for a, _b, _i, _o in coarse_edges} | {b for _a, b, _i, _o in coarse_edges}
    axis_map = {(v, leg): (bi, axis) for bi, boundary in enumerate(boundaries)
                for axis, (v, leg) in enumerate(boundary)}
    outputs = [axis_map[k] for k in sorted(axis_map) if axis_map[k] not in occupied]
    return _contract(tensors, specs, outputs).reshape(-1)


def compare_micro_and_fused_qfi(node_states, edges, singlet, blocks, generators=None):
    if generators is None:
        generators = spin_generators(singlet.shape[0] - 1)
    fused, boundaries, coarse_edges = fuse_vertex_blocks(node_states, edges, singlet, blocks)
    micro = _micro_state(node_states, edges, singlet)
    coarse = _coarse_state(fused, boundaries, coarse_edges, singlet)
    nm, nc = np.linalg.norm(micro), np.linalg.norm(coarse)
    if nm <= 1e-14 or nc <= 1e-14:
        raise ValueError("zero-support boundary state")
    pm, pc = micro / nm, coarse / nc
    overlap = np.vdot(pc, pm)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-14 else 1.0
    amp_error = float(np.linalg.norm(pm - phase * pc))
    dm_all, dc_all = [], []
    for _a, _b, ei, _original in coarse_edges:
        for gen in generators:
            dm = _micro_state(node_states, edges, singlet, (ei, gen)) / nm
            dc = _coarse_state(fused, boundaries, coarse_edges, singlet, (ei, gen)) / nc
            dm -= pm * np.vdot(pm, dm)
            dc -= pc * np.vdot(pc, dc)
            dm_all.append(dm)
            dc_all.append(dc)
    def make_qfi(ds):
        return 4.0 * np.real([[np.vdot(a, b) for b in ds] for a in ds])
    fm, fc = make_qfi(dm_all), make_qfi(dc_all)
    error = float(np.max(np.abs(fm - fc))) if fm.size else 0.0
    return {"amplitude_error_up_to_global_phase": amp_error, "micro_qfi": fm.tolist(),
            "coarse_qfi": fc.tolist(), "qfi_max_abs_error": error,
            "retained_edge_indices": [e[2] for e in coarse_edges],
            "coarse_parameter_dimension": len(dm_all),
            "passed": amp_error < 2e-9 and error < 2e-8}


def deterministic_case(name="path_4", two_j=1, seed=20261118, blocks=None):
    nodes, edges = graph_fixtures()[name]
    rng = np.random.default_rng(seed)
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    if blocks is None:
        blocks = tuple((v,) for v in range(nodes))
    return compare_micro_and_fused_qfi(states, edges, singlet, blocks)
