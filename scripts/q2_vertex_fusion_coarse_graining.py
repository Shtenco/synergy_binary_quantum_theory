#!/usr/bin/env python3
"""Exact vertex-block fusion and boundary-QFI equivalence for finite BQG graphs."""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np

from q2_relational_higher_spin_graph_qfi_gate import (
    graph_fixtures,
    invariant_singlet_tensor,
    node_generator,
    random_intertwiner,
    spin_generators,
)


def _check_partition(node_count: int, blocks: Iterable[Iterable[int]]):
    blocks = tuple(tuple(int(v) for v in block) for block in blocks)
    flat = [v for block in blocks for v in block]
    if any(not block for block in blocks):
        raise ValueError("vertex blocks must be non-empty")
    if sorted(flat) != list(range(node_count)) or len(flat) != len(set(flat)):
        raise ValueError("vertex blocks must partition every vertex exactly once")
    return blocks


def _network_contract(tensors, edge_specs, output_legs):
    """Contract tensors with pairwise singlet edges; output_legs fixes free-axis order."""
    labels_by_tensor = [[None] * tensor.ndim for tensor in tensors]
    operands = []
    next_label = 0
    for edge in edge_specs:
        (ta, aa), (tb, ab), _singlet = edge
        if labels_by_tensor[ta][aa] is not None or labels_by_tensor[tb][ab] is not None:
            raise ValueError("a tensor leg is used by multiple edges")
        labels_by_tensor[ta][aa] = next_label
        labels_by_tensor[tb][ab] = next_label + 1
        next_label += 2
    output_labels = []
    for tensor_index, axis in output_legs:
        if labels_by_tensor[tensor_index][axis] is not None:
            raise ValueError("an output leg is already contracted")
        labels_by_tensor[tensor_index][axis] = next_label
        output_labels.append(next_label)
        next_label += 1
    for tensor, labels in zip(tensors, labels_by_tensor):
        if any(label is None for label in labels):
            raise ValueError("unclassified tensor axis")
        operands.extend([tensor, labels])
    for (ta, aa), (tb, ab), singlet in edge_specs:
        operands.extend([singlet, [labels_by_tensor[ta][aa], labels_by_tensor[tb][ab]]])
    operands.append(output_labels)
    return np.einsum(*operands, optimize=True)


def fuse_vertex_blocks(node_states, edges, singlet, blocks):
    """Fuse each vertex block, contracting only edges internal to that block.

    Returns fused tensors, their ordered original (vertex, leg) boundary labels,
    coarse edges, and a canonical global boundary ordering.
    """
    node_count = len(node_states)
    blocks = _check_partition(node_count, blocks)
    owner = {vertex: block_index for block_index, block in enumerate(blocks)
             for vertex in block}
    normalized_edges = []
    used = set()
    for edge_index, edge in enumerate(edges):
        if len(edge) != 4:
            raise ValueError("edges must be (u, leg_u, v, leg_v)")
        u, lu, v, lv = map(int, edge)
        if u == v or not (0 <= u < node_count and 0 <= v < node_count):
            raise ValueError("invalid edge endpoints")
        if not (0 <= lu < 4 and 0 <= lv < 4):
            raise ValueError("invalid half-link slots")
        for endpoint in ((u, lu), (v, lv)):
            if endpoint in used:
                raise ValueError("a half-link is used by multiple edges")
            used.add(endpoint)
        normalized_edges.append((u, lu, v, lv, edge_index))

    fused_tensors = []
    fused_boundaries = []
    for block_index, block in enumerate(blocks):
        local_nodes = sorted(block)
        local_index = {node: i for i, node in enumerate(local_nodes)}
        tensors = [np.asarray(node_states[node], dtype=complex).reshape(
            (singlet.shape[0],) * 4) for node in local_nodes]
        internal = []
        for u, lu, v, lv, _edge_index in normalized_edges:
            if owner[u] == block_index and owner[v] == block_index:
                internal.append(((local_index[u], lu), (local_index[v], lv), singlet))
        contracted = set()
        for (ta, aa), (tb, ab), _ in internal:
            contracted.add((local_nodes[ta], aa))
            contracted.add((local_nodes[tb], ab))
        boundary = [(node, leg) for node in local_nodes for leg in range(4)
                    if (node, leg) not in contracted]
        output = [(local_index[node], leg) for node, leg in boundary]
        fused = _network_contract(tensors, internal, output)
        fused_tensors.append(fused)
        fused_boundaries.append(boundary)

    boundary_owner = {
        original_leg: (block_index, axis)
        for block_index, boundary in enumerate(fused_boundaries)
        for axis, original_leg in enumerate(boundary)
    }
    coarse_edges = []
    for u, lu, v, lv, edge_index in normalized_edges:
        bu, bv = owner[u], owner[v]
        if bu != bv:
            coarse_edges.append((
                (bu, boundary_owner[(u, lu)][1]),
                (bv, boundary_owner[(v, lv)][1]),
                edge_index,
                (u, lu, v, lv),
            ))
    global_boundary = sorted(
        leg for boundary in fused_boundaries for leg in boundary
        if all(not ((u == leg[0] and lu == leg[1]) or
                    (v == leg[0] and lv == leg[1]))
               for u, lu, v, lv, _ in normalized_edges)
    )
    return fused_tensors, fused_boundaries, coarse_edges, global_boundary


def _contract_micro(node_states, edges, singlet, twist=None):
    tensors = [np.asarray(state, complex).reshape((singlet.shape[0],) * 4).copy()
               for state in node_states]
    if twist is not None:
        edge_index, generator = twist
        u, lu, _v, _lv = edges[edge_index]
        tensors[u] = ( -1j * (node_generator(generator, lu, singlet.shape[0])
                              @ tensors[u].reshape(-1))).reshape(tensors[u].shape)
    specs = [((u, lu), (v, lv), singlet) for u, lu, v, lv in edges]
    contracted = {(u, lu) for u, lu, v, lv in edges} | {
        (v, lv) for u, lu, v, lv in edges
    }
    output = [(node, leg) for node in range(len(tensors)) for leg in range(4)
              if (node, leg) not in contracted]
    raw = _network_contract(tensors, specs, output).reshape(-1)
    return raw, output


def _contract_coarse(fused, boundaries, coarse_edges, singlet, twist=None):
    tensors = [tensor.copy() for tensor in fused]
    if twist is not None:
        edge_index, generator = twist
        edge = next(item for item in coarse_edges if item[2] == edge_index)
        (block, axis), _other, _idx, _original = edge
        tensors[block] = np.moveaxis(tensors[block], axis, 0)
        tensors[block] = (-1j * (node_generator(generator, 0, singlet.shape[0])
                                  @ tensors[block].reshape(singlet.shape[0], -1))).reshape(
                                      tensors[block].shape)
        tensors[block] = np.moveaxis(tensors[block], 0, axis)
    specs = [((a[0], a[1]), (b[0], b[1]), singlet)
             for a, b, _idx, _original in coarse_edges]
    contracted = {(a[0], a[1]) for a, _b, _idx, _orig in coarse_edges} | {
        (b[0], b[1]) for _a, b, _idx, _orig in coarse_edges
    }
    output = [(bi, axis) for bi, tensor in enumerate(tensors)
              for axis in range(tensor.ndim) if (bi, axis) not in contracted]
    # Canonical order is original vertex/leg order, not block storage order.
    axis_lookup = {(node, leg): (bi, axis)
                   for bi, boundary in enumerate(boundaries)
                   for axis, (node, leg) in enumerate(boundary)}
    canonical = [axis_lookup[key] for key in sorted(
        key for boundary in boundaries for key in boundary
        if key not in {(u, lu) for _a, _b, _i, (u, lu, v, lv) in coarse_edges}
        and key not in {(v, lv) for _a, _b, _i, (u, lu, v, lv) in coarse_edges}
    )]
    raw = _network_contract(tensors, specs, canonical).reshape(-1)
    return raw


def compare_micro_and_fused_qfi(node_states, edges, singlet, blocks, generators=None):
    """Compare normalized boundary amplitudes and QFI for all surviving edge twists."""
    if generators is None:
        generators = spin_generators(singlet.shape[0] - 1)
    fused, boundaries, coarse_edges, _ = fuse_vertex_blocks(
        node_states, edges, singlet, blocks
    )
    micro_raw, micro_order = _contract_micro(node_states, edges, singlet)
    coarse_raw = _contract_coarse(fused, boundaries, coarse_edges, singlet)
    if len(micro_order) != coarse_raw.size.bit_length() - 1:
        # Dimension check without assuming a particular spin representation.
        if micro_raw.size != coarse_raw.size:
            raise AssertionError("micro/coarse boundary dimensions differ")
    norm_micro = np.linalg.norm(micro_raw)
    norm_coarse = np.linalg.norm(coarse_raw)
    if norm_micro <= 1e-14 or norm_coarse <= 1e-14:
        raise ValueError("zero-support boundary state")
    psi_micro, psi_coarse = micro_raw / norm_micro, coarse_raw / norm_coarse
    amplitude_error = float(min(
        np.linalg.norm(psi_micro - psi_coarse),
        np.linalg.norm(psi_micro + psi_coarse),
    ))
    surviving = [item for item in coarse_edges]
    micro_derivatives, coarse_derivatives = [], []
    for _a, _b, edge_index, _original in surviving:
        for generator in generators:
            dm, _ = _contract_micro(node_states, edges, singlet, (edge_index, generator))
            dc = _contract_coarse(fused, boundaries, coarse_edges, singlet,
                                  (edge_index, generator))
            dm = dm / norm_micro
            dc = dc / norm_coarse
            dm -= psi_micro * np.vdot(psi_micro, dm)
            dc -= psi_coarse * np.vdot(psi_coarse, dc)
            micro_derivatives.append(dm)
            coarse_derivatives.append(dc)
    def qfi(derivatives):
        return 4.0 * np.real(np.array([
            [np.vdot(a, b) for b in derivatives] for a in derivatives
        ]))
    fm, fc = qfi(micro_derivatives), qfi(coarse_derivatives)
    return {
        "amplitude_error_up_to_global_sign": amplitude_error,
        "micro_qfi": fm,
        "coarse_qfi": fc,
        "qfi_max_abs_error": float(np.max(np.abs(fm - fc))) if fm.size else 0.0,
        "retained_edge_indices": [item[2] for item in surviving],
        "coarse_parameter_dimension": len(micro_derivatives),
        "passed": amplitude_error < 2e-9 and (
            not fm.size or np.allclose(fm, fc, atol=2e-8, rtol=2e-8)
        ),
    }


def deterministic_case(name="path_4", two_j=1, seed=20261118, blocks=None):
    nodes, edges = graph_fixtures()[name]
    rng = np.random.default_rng(seed)
    # Local import avoids duplicating the intertwiner construction.
    from q2_relational_higher_spin_network_qfi_gate import (
        four_valent_intertwiner_basis,
    )
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    if blocks is None:
        blocks = tuple((node,) for node in range(nodes))
    return compare_micro_and_fused_qfi(states, edges, singlet, blocks)
