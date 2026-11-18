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



def _reduced_density_and_derivative(state, derivative, dimension, keep_axes):
    """Trace all but selected tensor axes from a normalized pure state."""
    axis_count = int(round(np.log(state.size) / np.log(dimension)))
    keep_axes = tuple(keep_axes)
    traced_axes = tuple(axis for axis in range(axis_count) if axis not in keep_axes)
    permutation = keep_axes + traced_axes
    shape = (dimension,) * axis_count
    matrix = state.reshape(shape).transpose(permutation).reshape(
        dimension ** len(keep_axes), -1
    )
    derivative_matrix = derivative.reshape(shape).transpose(permutation).reshape(
        dimension ** len(keep_axes), -1
    )
    rho = matrix @ matrix.conj().T
    drho = derivative_matrix @ matrix.conj().T + matrix @ derivative_matrix.conj().T
    return rho, drho


def _sld_mixed_qfi(rho, derivatives, eigenvalue_tolerance=1e-12):
    """Symmetric-logarithmic-derivative QFI for a finite density matrix."""
    eigenvalues, eigenvectors = np.linalg.eigh(rho)
    eigenvalues = np.maximum(eigenvalues.real, 0.0)
    denominators = eigenvalues[:, None] + eigenvalues[None, :]
    weights = np.zeros_like(denominators)
    support = denominators > eigenvalue_tolerance
    weights[support] = 2.0 / denominators[support]
    transformed = [
        eigenvectors.conj().T @ derivative @ eigenvectors
        for derivative in derivatives
    ]
    qfi = np.zeros((len(derivatives), len(derivatives)), dtype=float)
    for a in range(len(derivatives)):
        for b in range(len(derivatives)):
            qfi[a, b] = float(np.real(np.sum(
                weights * transformed[a] * transformed[b].T
            )))
    return 0.5 * (qfi + qfi.T)


def compare_micro_and_fused_reduced_qfi(
    node_states, edges, singlet, blocks, keep_boundary_keys, generators=None
):
    """Compare mixed-boundary SLD QFI before/after exact vertex fusion.

    All boundary legs not listed in keep_boundary_keys are traced out. The
    retained parameters are the same inter-block edge twists used by the
    pure-boundary comparison. Also checks QFI monotonicity under this partial
    trace (quantum data-processing inequality).
    """
    if generators is None:
        generators = spin_generators(singlet.shape[0] - 1)
    fused, boundaries, coarse_edges = fuse_vertex_blocks(
        node_states, edges, singlet, blocks
    )
    dimension = singlet.shape[0]
    occupied = {
        (u, lu) for u, lu, _v, _lv in edges
    } | {
        (v, lv) for _u, _lu, v, lv in edges
    }
    boundary_keys = [
        (v, leg)
        for v in range(len(node_states))
        for leg in range(4)
        if (v, leg) not in occupied
    ]
    keep_boundary_keys = tuple(tuple(key) for key in keep_boundary_keys)
    if not keep_boundary_keys or len(set(keep_boundary_keys)) != len(keep_boundary_keys):
        raise ValueError("kept boundary keys must be non-empty and unique")
    if any(key not in boundary_keys for key in keep_boundary_keys):
        raise ValueError("kept boundary key is not an external half-link")
    keep_axes = tuple(boundary_keys.index(key) for key in keep_boundary_keys)

    micro_raw = _micro_state(node_states, edges, singlet)
    coarse_raw = _coarse_state(fused, boundaries, coarse_edges, singlet)
    micro_norm, coarse_norm = np.linalg.norm(micro_raw), np.linalg.norm(coarse_raw)
    if micro_norm <= 1e-14 or coarse_norm <= 1e-14:
        raise ValueError("zero-support boundary state")
    micro_state, coarse_state = micro_raw / micro_norm, coarse_raw / coarse_norm

    micro_derivatives, coarse_derivatives = [], []
    for _a, _b, edge_index, _original in coarse_edges:
        for generator in generators:
            dm = _micro_state(
                node_states, edges, singlet, (edge_index, generator)
            ) / micro_norm
            dc = _coarse_state(
                fused, boundaries, coarse_edges, singlet, (edge_index, generator)
            ) / coarse_norm
            dm -= micro_state * np.vdot(micro_state, dm)
            dc -= coarse_state * np.vdot(coarse_state, dc)
            micro_derivatives.append(dm)
            coarse_derivatives.append(dc)

    zero_micro = np.zeros_like(micro_state)
    zero_coarse = np.zeros_like(coarse_state)
    micro_rho, _ = _reduced_density_and_derivative(
        micro_state, zero_micro, dimension, keep_axes
    )
    coarse_rho, _ = _reduced_density_and_derivative(
        coarse_state, zero_coarse, dimension, keep_axes
    )
    micro_drhos, coarse_drhos = [], []
    for derivative in micro_derivatives:
        _rho, drho = _reduced_density_and_derivative(
            micro_state, derivative, dimension, keep_axes
        )
        micro_drhos.append(drho)
    for derivative in coarse_derivatives:
        _rho, drho = _reduced_density_and_derivative(
            coarse_state, derivative, dimension, keep_axes
        )
        coarse_drhos.append(drho)

    micro_qfi = _sld_mixed_qfi(micro_rho, micro_drhos)
    coarse_qfi = _sld_mixed_qfi(coarse_rho, coarse_drhos)
    pure = np.asarray(compare_micro_and_fused_qfi(
        node_states, edges, singlet, blocks, generators
    )["micro_qfi"], dtype=float).reshape(len(micro_derivatives), len(micro_derivatives))
    qfi_gap = 0.5 * ((pure - micro_qfi) + (pure - micro_qfi).T)
    data_processing_min_eigenvalue = (
        float(np.linalg.eigvalsh(qfi_gap)[0]) if qfi_gap.size else 0.0
    )
    density_error = float(np.max(np.abs(micro_rho - coarse_rho)))
    derivative_error = max(
        (float(np.max(np.abs(a - b))) for a, b in zip(micro_drhos, coarse_drhos)),
        default=0.0,
    )
    qfi_error = float(np.max(np.abs(micro_qfi - coarse_qfi))) if micro_qfi.size else 0.0
    return {
        "kept_boundary_keys": [list(key) for key in keep_boundary_keys],
        "kept_boundary_dimension": int(dimension ** len(keep_axes)),
        "traced_boundary_leg_count": len(boundary_keys) - len(keep_axes),
        "retained_edge_indices": [edge[2] for edge in coarse_edges],
        "micro_reduced_density_trace_error": abs(float(np.trace(micro_rho).real) - 1.0),
        "coarse_reduced_density_trace_error": abs(float(np.trace(coarse_rho).real) - 1.0),
        "reduced_density_max_abs_error": density_error,
        "reduced_derivative_max_abs_error": derivative_error,
        "micro_mixed_qfi": micro_qfi.tolist(),
        "coarse_mixed_qfi": coarse_qfi.tolist(),
        "mixed_qfi_max_abs_error": qfi_error,
        "pure_minus_mixed_qfi_min_eigenvalue": data_processing_min_eigenvalue,
        "passed": (
            density_error < 2e-9
            and derivative_error < 2e-8
            and qfi_error < 2e-8
            and data_processing_min_eigenvalue > -2e-8
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



def _bell_number(n):
    """Return the nth Bell number using an integer-only recurrence."""
    row = [1]
    for _ in range(n):
        next_row = [row[-1]]
        for k in range(len(row)):
            next_row.append(next_row[-1] + row[k])
        row = next_row
    return row[0]


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
        "expected_partition_count": len(partitions),
        "includes_block_order_reversal": True,
        "max_amplitude_error_up_to_global_phase": max_amplitude_error,
        "max_qfi_error": max_qfi_error,
        "max_fused_block_gauss_relative_residual": max_gauss_residual,
        "passed": (
            len(partitions) == _bell_number(nodes)
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





def exhaustive_reduced_boundary_partition_control(
    graph_name, two_j, keep_boundary_keys, seed=20261202
):
    """Exhaust all graph partitions for a fixed reduced boundary subsystem."""
    if graph_name not in graph_fixtures():
        raise ValueError(f"unknown graph fixture: {graph_name}")
    nodes, edges = graph_fixtures()[graph_name]
    rng = np.random.default_rng(seed + two_j)
    basis, _ = four_valent_intertwiner_basis(two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    partitions = list(_set_partitions(range(nodes)))
    maxima = {
        "density": 0.0,
        "derivative": 0.0,
        "qfi": 0.0,
    }
    minimum_data_processing_eigenvalue = float("inf")
    for blocks in partitions:
        result = compare_micro_and_fused_reduced_qfi(
            states, edges, singlet, blocks, keep_boundary_keys
        )
        reversed_blocks = tuple(tuple(reversed(block)) for block in reversed(blocks))
        reordered = compare_micro_and_fused_reduced_qfi(
            states, edges, singlet, reversed_blocks, keep_boundary_keys
        )
        if not result["passed"] or not reordered["passed"]:
            return {
                "graph": graph_name,
                "two_j": int(two_j),
                "partition_count": len(partitions),
                "expected_partition_count": _bell_number(nodes),
                "passed": False,
                "failed_partition": blocks,
                "reordered_partition": reversed_blocks,
            }
        maxima["density"] = max(
            maxima["density"],
            result["reduced_density_max_abs_error"],
            reordered["reduced_density_max_abs_error"],
        )
        maxima["derivative"] = max(
            maxima["derivative"],
            result["reduced_derivative_max_abs_error"],
            reordered["reduced_derivative_max_abs_error"],
        )
        maxima["qfi"] = max(
            maxima["qfi"],
            result["mixed_qfi_max_abs_error"],
            reordered["mixed_qfi_max_abs_error"],
        )
        minimum_data_processing_eigenvalue = min(
            minimum_data_processing_eigenvalue,
            result["pure_minus_mixed_qfi_min_eigenvalue"],
            reordered["pure_minus_mixed_qfi_min_eigenvalue"],
        )
    return {
        "graph": graph_name,
        "two_j": int(two_j),
        "spin": float(two_j) / 2.0,
        "partition_count": len(partitions),
        "expected_partition_count": _bell_number(nodes),
        "includes_block_order_reversal": True,
        "kept_boundary_keys": [list(key) for key in keep_boundary_keys],
        "max_reduced_density_error": maxima["density"],
        "max_reduced_derivative_error": maxima["derivative"],
        "max_mixed_qfi_error": maxima["qfi"],
        "minimum_pure_minus_mixed_qfi_eigenvalue": minimum_data_processing_eigenvalue,
        "passed": (
            len(partitions) == _bell_number(nodes)
            and maxima["density"] < 2e-9
            and maxima["derivative"] < 2e-8
            and maxima["qfi"] < 2e-8
            and minimum_data_processing_eigenvalue > -2e-8
        ),
    }


def run_reduced_boundary_suite(two_j_values=(1, 2)):
    """Compare mixed-state QFI after tracing boundary legs on two graph families."""
    specifications = {
        "path_3": {
            "keep_boundary_keys": ((0, 0), (0, 1)),
        },
        "cycle_4": {
            "keep_boundary_keys": ((0, 1), (0, 3)),
        },
    }
    controls = {}
    for two_j in two_j_values:
        controls[str(two_j / 2.0)] = {}
        for graph_name, spec in specifications.items():
            nodes, _edges = graph_fixtures()[graph_name]
            controls[str(two_j / 2.0)][graph_name] = (
                exhaustive_reduced_boundary_partition_control(
                    graph_name, two_j, spec["keep_boundary_keys"],
                    seed=20261201 + nodes
                )
            )
    return {
        "status": "mixed reduced-boundary QFI under exact vertex fusion",
        "tested_spins": [two_j / 2.0 for two_j in two_j_values],
        "graphs": list(specifications),
        "passed": all(
            row["passed"]
            for spin_rows in controls.values()
            for row in spin_rows.values()
        ),
        "controls": controls,
        "claim_boundary": (
            "Exact reduced density matrix and SLD-QFI equivalence after tracing "
            "the same external boundary legs, plus the QFI data-processing "
            "inequality for the tested finite states. This is not a physical "
            "renormalization-group flow."
        ),
    }


def run_high_spin_small_graph_suite(two_j_values=(3, 4), graph_names=("path_3", "cycle_4")):
    """Exhaust all vertex partitions for small-boundary graphs at higher spins."""
    controls = {}
    for two_j in two_j_values:
        controls[str(two_j / 2.0)] = {
            name: exhaustive_graph_partition_control(
                name, two_j=two_j, seed=20261130 + two_j
            )
            for name in graph_names
        }
    return {
        "status": "higher-spin exhaustive vertex-fusion partition control",
        "tested_spins": [two_j / 2.0 for two_j in two_j_values],
        "graphs": list(graph_names),
        "partition_order_evaluations": sum(
            2 * row["partition_count"]
            for spin_rows in controls.values()
            for row in spin_rows.values()
        ),
        "passed": all(
            row["passed"]
            for spin_rows in controls.values()
            for row in spin_rows.values()
        ),
        "controls": controls,
        "claim_boundary": (
            "Finite exact vertex fusion and retained-edge pure-state QFI for "
            "the listed graph fixtures and spin labels only; no RG fixed point "
            "or continuum statement is implied."
        ),
    }


def main():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-j", type=int, default=1)
    parser.add_argument("--high-spin-small-graphs", action="store_true")
    parser.add_argument("--reduced-boundary-qfi", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.two_j < 1:
        parser.error("--two-j must be a positive integer")
    if args.high_spin_small_graphs and args.reduced_boundary_qfi:
        parser.error("choose only one specialized control mode")
    if args.high_spin_small_graphs:
        result = run_high_spin_small_graph_suite()
    elif args.reduced_boundary_qfi:
        result = run_reduced_boundary_suite()
    else:
        result = run_coarse_graining_suite(args.two_j)
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
