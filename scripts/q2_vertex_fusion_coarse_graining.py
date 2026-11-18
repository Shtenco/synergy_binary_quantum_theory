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



def fused_block_gauss_residuals(fused_tensors, generators):
    """Relative residuals of the diagonal SU(2) Gauss generators on each fused block."""
    residuals = []
    for tensor in fused_tensors:
        tensor = np.asarray(tensor, complex)
        norm = float(np.linalg.norm(tensor))
        block_residual = 0.0
        for generator in generators:
            total = np.zeros_like(tensor, dtype=complex)
            for axis in range(tensor.ndim):
                total += _apply_one_leg(tensor, axis, generator)
            relative = float(np.linalg.norm(total)) / max(norm, 1e-14)
            block_residual = max(block_residual, relative)
        residuals.append(block_residual)
    return residuals

def compare_micro_and_fused_qfi(node_states, edges, singlet, blocks, generators=None):
    if generators is None:
        generators = spin_generators(singlet.shape[0] - 1)
    fused, boundaries, coarse_edges = fuse_vertex_blocks(node_states, edges, singlet, blocks)
    gauss_residuals = fused_block_gauss_residuals(fused, generators)
    gauss_max = max(gauss_residuals, default=0.0)
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
            "fused_block_gauss_relative_residuals": gauss_residuals,
            "gauss_max_relative_residual": gauss_max,
            "passed": amp_error < 2e-9 and error < 2e-8 and gauss_max < 2e-9}


def deterministic_case(name="path_4", two_j=1, seed=20261118, blocks=None):
    nodes, edges = graph_fixtures()[name]
    rng = np.random.default_rng(seed)
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    if blocks is None:
        blocks = tuple((v,) for v in range(nodes))
    return compare_micro_and_fused_qfi(states, edges, singlet, blocks)


def fuse_fused_blocks(fused_tensors, boundaries, coarse_edges, blocks, singlet):
    """Fuse already-coarse tensors once more, preserving original boundary labels."""
    n = len(fused_tensors)
    if len(boundaries) != n:
        raise ValueError("coarse boundary metadata must match tensor count")
    blocks = _partition(n, blocks)
    owner = {v: bi for bi, block in enumerate(blocks) for v in block}
    normalized, used = [], set()
    for edge in coarse_edges:
        if len(edge) != 4:
            raise ValueError("coarse edges must contain endpoints, original index, and metadata")
        a, b, edge_index, original = edge
        if len(a) != 2 or len(b) != 2:
            raise ValueError("coarse edge endpoints must be (tensor_index, axis)")
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer))
               for x in (*a, *b, edge_index)):
            raise ValueError("coarse edge indices must be integers")
        a, b = (int(a[0]), int(a[1])), (int(b[0]), int(b[1]))
        if a[0] == b[0] or not (0 <= a[0] < n and 0 <= b[0] < n):
            raise ValueError("invalid coarse edge tensor endpoints")
        if not (0 <= a[1] < np.ndim(fused_tensors[a[0]])
                and 0 <= b[1] < np.ndim(fused_tensors[b[0]])):
            raise ValueError("invalid coarse edge tensor axes")
        if a in used or b in used:
            raise ValueError("a coarse half-link is used more than once")
        used.update((a, b))
        normalized.append((a, b, int(edge_index), original))

    fused, new_boundaries, remap = [], [], {}
    for block_index, block in enumerate(blocks):
        local = sorted(block)
        local_index = {v: i for i, v in enumerate(local)}
        tensors = [np.asarray(fused_tensors[v], complex) for v in local]
        internal = [
            ((local_index[a[0]], a[1]), (local_index[b[0]], b[1]), singlet)
            for a, b, _edge_index, _original in normalized
            if owner[a[0]] == block_index == owner[b[0]]
        ]
        contracted_axes = {
            endpoint
            for a, b, _edge_index, _original in normalized
            if owner[a[0]] == block_index == owner[b[0]]
            for endpoint in (a, b)
        }
        open_axes = [
            (v, axis)
            for v in local
            for axis in range(np.ndim(fused_tensors[v]))
            if (v, axis) not in contracted_axes
        ]
        fused.append(_contract(
            tensors, internal,
            [(local_index[v], axis) for v, axis in open_axes],
        ))
        new_boundaries.append([boundaries[v][axis] for v, axis in open_axes])
        for new_axis, old_endpoint in enumerate(open_axes):
            remap[old_endpoint] = (block_index, new_axis)

    retained = [
        (remap[a], remap[b], edge_index, original)
        for a, b, edge_index, original in normalized
        if owner[a[0]] != owner[b[0]]
    ]
    return fused, new_boundaries, retained

def _coarse_qfi(fused, boundaries, coarse_edges, singlet, generators):
    raw = _coarse_state(fused, boundaries, coarse_edges, singlet)
    norm = np.linalg.norm(raw)
    if norm <= 1e-14:
        raise ValueError("zero-support coarse boundary state")
    state = raw / norm
    derivatives = []
    for _a, _b, edge_index, _original in coarse_edges:
        for generator in generators:
            derivative = _coarse_state(
                fused, boundaries, coarse_edges, singlet, (edge_index, generator)
            ) / norm
            derivative -= state * np.vdot(state, derivative)
            derivatives.append(derivative)
    qfi = 4.0 * np.real([[np.vdot(a, b) for b in derivatives] for a in derivatives])
    return state, qfi


def compare_staged_and_direct_fusion(node_states, edges, singlet,
                                     first_blocks, second_blocks, direct_blocks):
    """Check that two successive blocking steps equal one direct blocking step."""
    first = fuse_vertex_blocks(node_states, edges, singlet, first_blocks)
    staged = fuse_fused_blocks(*first, second_blocks, singlet)
    direct = fuse_vertex_blocks(node_states, edges, singlet, direct_blocks)
    generators = spin_generators(singlet.shape[0] - 1)
    staged_state, staged_qfi = _coarse_qfi(*staged, singlet, generators)
    direct_state, direct_qfi = _coarse_qfi(*direct, singlet, generators)
    overlap = np.vdot(direct_state, staged_state)
    phase = overlap / abs(overlap) if abs(overlap) > 1e-14 else 1.0
    amplitude_error = float(np.linalg.norm(staged_state - phase * direct_state))
    qfi_error = float(np.max(np.abs(staged_qfi - direct_qfi))) if staged_qfi.size else 0.0
    staged_edges = [edge[2] for edge in staged[2]]
    direct_edges = [edge[2] for edge in direct[2]]
    return {
        "amplitude_error_up_to_global_phase": amplitude_error,
        "qfi_max_abs_error": qfi_error,
        "staged_retained_edge_indices": staged_edges,
        "direct_retained_edge_indices": direct_edges,
        "passed": (
            amplitude_error < 2e-9 and qfi_error < 2e-8
            and staged_edges == direct_edges
        ),
    }


def _set_partitions(items):
    """Generate every set partition exactly once."""
    items = tuple(items)
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in _set_partitions(rest):
        yield ((first,),) + partition
        for index in range(len(partition)):
            yield (
                partition[:index]
                + ((first,) + partition[index],)
                + partition[index + 1:]
            )


def exhaustive_graph_partition_control(graph_name="path_5", two_j=1, seed=20261120):
    """Check every vertex partition of a fixture and invariance under block ordering."""
    if graph_name not in graph_fixtures():
        raise ValueError(f"unknown graph fixture: {graph_name}")
    nodes, edges = graph_fixtures()[graph_name]
    rng = np.random.default_rng(seed)
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    partitions = list(_set_partitions(range(nodes)))
    max_amplitude_error = 0.0
    max_qfi_error = 0.0
    max_gauss_residual = 0.0
    for blocks in partitions:
        result = compare_micro_and_fused_qfi(states, edges, singlet, blocks)
        reordered_blocks = tuple(tuple(reversed(block)) for block in reversed(blocks))
        reordered = compare_micro_and_fused_qfi(
            states, edges, singlet, reordered_blocks
        )
        if not result["passed"] or not reordered["passed"]:
            return {
                "partition_count": len(partitions),
                "expected_partition_count": len(partitions),
                "passed": False,
                "failed_partition": blocks,
                "reordered_partition": reordered_blocks,
            }
        if result["retained_edge_indices"] != reordered["retained_edge_indices"]:
            return {
                "partition_count": len(partitions),
                "expected_partition_count": len(partitions),
                "passed": False,
                "failed_partition": blocks,
                "reason": "retained edge identity depends on block ordering",
            }
        max_amplitude_error = max(
            max_amplitude_error,
            result["amplitude_error_up_to_global_phase"],
            reordered["amplitude_error_up_to_global_phase"],
        )
        max_qfi_error = max(
            max_qfi_error,
            result["qfi_max_abs_error"],
            reordered["qfi_max_abs_error"],
        )
        max_gauss_residual = max(
            max_gauss_residual,
            result["gauss_max_relative_residual"],
            reordered["gauss_max_relative_residual"],
        )
    return {
        "graph": graph_name,
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "partition_count": len(partitions),
        "expected_partition_count": 52,
        "includes_block_order_reversal": True,
        "max_amplitude_error_up_to_global_phase": max_amplitude_error,
        "max_qfi_error": max_qfi_error,
        "max_fused_block_gauss_relative_residual": max_gauss_residual,
        "passed": (
            len(partitions) == 52
            and max_amplitude_error < 2e-9
            and max_qfi_error < 2e-8
            and max_gauss_residual < 2e-9
        ),
    }


def exhaustive_path5_partition_control(two_j=1, seed=20261120):
    """Backward-compatible path_5 exhaustive partition control."""
    return exhaustive_graph_partition_control("path_5", two_j=two_j, seed=seed)

def run_coarse_graining_suite(two_j=1):
    """Run the deterministic block-fusion controls on all supported graph fixtures."""
    cases = {
        "path_3": ((0, 1), (2,)),
        "path_4": ((0, 1), (2, 3)),
        "path_5": ((0, 1), (2, 3, 4)),
        "star_4": ((0, 1), (2, 3)),
        "cycle_4": ((0, 1), (2, 3)),
    }
    results = {
        name: deterministic_case(name, two_j=two_j, blocks=blocks)
        for name, blocks in cases.items()
    }
    nodes, edges = graph_fixtures()["path_5"]
    rng = np.random.default_rng(20261119 + two_j)
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    staged = compare_staged_and_direct_fusion(
        states, edges, invariant_singlet_tensor(two_j),
        first_blocks=((0, 1), (2,), (3, 4)),
        second_blocks=((0, 1), (2,)),
        direct_blocks=((0, 1, 2), (3, 4)),
    )
    exhaustive = {
        name: exhaustive_graph_partition_control(
            name, two_j=two_j, seed=20261120 + two_j
        )
        for name in graph_fixtures()
    }
    return {
        "status": "exact finite vertex-fusion coarse-graining",
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "passed": (
            all(row["passed"] for row in results.values())
            and staged["passed"]
            and all(row["passed"] for row in exhaustive.values())
        ),
        "exhaustive_graph_partition_controls": exhaustive,
        "staged_vs_direct_two_level_fusion": staged,
        "claim_boundary": (
            "Exact equality of finite tensor contractions and retained-edge pure-state "
            "QFI under vertex-block fusion. This is not an RG fixed-point or continuum result."
        ),
        "cases": results,
    }


def main():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-j", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.two_j < 1:
        parser.error("--two-j must be a positive integer")
    result = run_coarse_graining_suite(args.two_j)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
