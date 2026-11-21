#!/usr/bin/env python3
"""Build the exact 30-column strict-interior five-block BQG P carrier.

This is a production bridge, not a new physical assumption.  It reuses the
frozen L1 Peter-Weyl Euclidean dynamics from research/bcqg-core-candidate-v1.
For the canonical parent block 0 and its four coarse-face neighbours it:

  1. constructs the 24 strict-interior q=4 source columns in the FULL global
     L1 Gauss basis (no boundary contraction and no compact local relabelling),
  2. forms the existing six parity/coset coarse-edge lifts per parent,
  3. verifies the true microscopic Hilbert Gram in each parent,
  4. assembles all 5 x 6 = 30 normalized columns and verifies their global Gram,
  5. serializes the actual global basis labels and amplitudes to compressed NPZ.

Because the strict channel changes only links internal to one parent and leaves
all boundary links/exterior labels exactly at background, columns belonging to
different parents must have disjoint global Gauss-basis support.  Therefore
shared-face SU(2) recoupling is NOT inserted into P.  Complete face recoupling
is required later for E/S images that actually hit/cross shared faces.

No ADM/GR ratio, TT projector, Wilson coefficient, dispersion target or fitted
physical scale enters this producer.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import itertools
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Mapping, MutableMapping

import numpy as np

EXPECTED_UPSTREAM_SHA = "2b2c9f623544f5d38f7ffd7f37617f91f4dae306"
JMAX2 = 3
TOL = 1e-10
GRAM_TOL = 1e-11


def add(dst: MutableMapping, src: Mapping, scale: complex = 1.0,
        tol: float = TOL) -> None:
    for key, amp in src.items():
        z = dst.get(key, 0j) + scale * amp
        if abs(z) > tol:
            dst[key] = z
        elif key in dst:
            del dst[key]


def inner(a: Mapping, b: Mapping) -> complex:
    if len(a) > len(b):
        return np.conjugate(inner(b, a))
    return sum((np.conjugate(z) * b.get(k, 0j) for k, z in a.items()), 0j)


def parity(p) -> int:
    inv = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inv % 2 else 1


def verify_upstream(root: Path, expected_sha: str) -> str:
    got = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if got != expected_sha:
        raise RuntimeError(
            f"frozen upstream SHA mismatch: expected {expected_sha}, got {got}"
        )
    return got


def load_upstream(root: Path):
    scripts = root / "scripts"
    if not scripts.is_dir():
        raise RuntimeError(f"missing frozen upstream scripts directory: {scripts}")
    sys.path.insert(0, str(scripts))
    return {
        "ZVM": importlib.import_module(
            "peter_weyl_zeroaware_volume_migration_experiment"
        ),
        "BE": importlib.import_module(
            "collective_barycentric_E_boundary_support_gate"
        ),
        "PL": importlib.import_module("pl_dual_complex"),
        "E": importlib.import_module("pl_peter_weyl_euclidean_local"),
    }


def face_neighbours(D, parent, parent_id: int):
    inside = {v for v, p in enumerate(parent) if p == parent_id}
    out = set()
    for a, b in D.dual_edges():
        if (a in inside) ^ (b in inside):
            outside = b if a in inside else a
            out.add(int(parent[outside]))
    return sorted(out)


def build_parent_sources(D, G, parent, parent_id: int, seed):
    inside = {v for v, p in enumerate(parent) if p == parent_id}
    nodes = sorted(inside)
    internal_edges = sorted(
        e for e in G.EDGES if e[0] in inside and e[1] in inside
    )
    boundary_edges = sorted(
        e for e in G.EDGES if (e[0] in inside) ^ (e[1] in inside)
    )
    internal_global = {G.EIDX[e] for e in internal_edges}
    boundary_global = {G.EIDX[e] for e in boundary_edges}

    if len(nodes) != 24 or len(internal_edges) != 36 or len(boundary_edges) != 24:
        raise RuntimeError(
            ("unexpected barycentric parent geometry", parent_id,
             len(nodes), len(internal_edges), len(boundary_edges))
        )

    source_columns = []
    rows = []
    for local_index, u in enumerate(nodes):
        strict = []
        q4_total = 0
        for sign, spec in G.oriented_specs(u):
            v, ra, rb, rc = spec
            path = D.plaquette_path(v, ra, rb)
            if len(path) - 1 != 4:
                continue
            q4_total += 1
            source_neighbor = D.neighbor[(v, rc)]
            if set(path[:-1]) <= inside and source_neighbor in inside:
                strict.append((sign, spec, path))
        if len(strict) != 1:
            raise RuntimeError(
                ("expected one strict-interior q4 spec", parent_id,
                 local_index, u, len(strict))
            )

        sign, spec, _ = strict[0]
        raw = {}
        add(raw, dict(G.T_items(seed, *spec, JMAX2, False)), -0.5j * sign)
        add(raw, dict(G.T_items(seed, *spec, JMAX2, True)), +0.5j * sign)
        p4 = {
            k: a for k, a in raw.items()
            if abs(a) > TOL and sum(s != 1 for s in k[0]) == 4
        }
        if not p4:
            raise RuntimeError(("empty strict q4 P4 column", parent_id, local_index))

        boundary_clean = True
        exterior_clean = True
        full_label_lengths = True
        for (spins, Ks), amp in p4.items():
            full_label_lengths &= len(spins) == len(G.EDGES) and len(Ks) == D.n_tets
            changed = {i for i, s in enumerate(spins) if s != 1}
            if changed & boundary_global:
                boundary_clean = False
            if not changed <= internal_global:
                exterior_clean = False
            if any(Ks[v] != 0 for v in range(D.n_tets) if v not in inside):
                exterior_clean = False
            if not np.isfinite([amp.real, amp.imag]).all():
                raise RuntimeError(("non-finite amplitude", parent_id, local_index))

        source_columns.append(p4)
        rows.append({
            "parent": parent_id,
            "local_index": local_index,
            "global_node": int(u),
            "q4_specs_total": q4_total,
            "strict_specs": len(strict),
            "support": len(p4),
            "norm": math.sqrt(float(sum(abs(a) ** 2 for a in p4.values()))),
            "boundary_spins_unchanged": bool(boundary_clean),
            "exterior_exactly_background_by_labels": bool(exterior_clean),
            "full_global_label_lengths": bool(full_label_lengths),
        })

        # Same cache discipline as the frozen strict-boundary producer.
        G.primitive_items.cache_clear()
        G.T_items.cache_clear()
        G.oriented_intertwiner.cache_clear()

    return nodes, source_columns, rows


def six_edge_lifts(coarse_tet, source_columns):
    perms = tuple(itertools.permutations(range(4)))
    if len(source_columns) != len(perms):
        raise RuntimeError("source ordering no longer matches the 24 S4 chambers")

    local_groups = [tuple(sorted(p[:2])) for p in perms]
    edges = sorted(set(local_groups))
    if len(edges) != 6 or any(local_groups.count(e) != 4 for e in edges):
        raise RuntimeError("invalid six-edge S4 coset partition")

    lifts = [dict() for _ in edges]
    for p, e, col in zip(perms, local_groups, source_columns):
        add(lifts[edges.index(e)], col, parity(p) / 2.0)

    raw_gram = np.empty((6, 6), complex)
    for i in range(6):
        for j in range(6):
            raw_gram[i, j] = inner(lifts[i], lifts[j])
    raw_gram = 0.5 * (raw_gram + raw_gram.conj().T)
    nu = float(np.mean(np.diag(raw_gram).real))
    if not nu > 0:
        raise RuntimeError(f"non-positive six-edge common norm: {nu}")
    defect = float(
        np.linalg.norm(raw_gram - nu * np.eye(6)) /
        max(np.linalg.norm(raw_gram), 1e-300)
    )
    if defect >= GRAM_TOL:
        raise RuntimeError(f"six-edge microscopic Gram is not common-norm I6: {defect}")

    scale = 1.0 / math.sqrt(nu)
    normalized = [{k: scale * z for k, z in c.items()} for c in lifts]
    actual_edges = [
        tuple(sorted((int(coarse_tet[e[0]]), int(coarse_tet[e[1]]))))
        for e in edges
    ]
    return normalized, actual_edges, nu, defect


def full_gram(columns):
    n = len(columns)
    K = np.empty((n, n), complex)
    for i in range(n):
        K[i, i] = inner(columns[i], columns[i])
        for j in range(i + 1, n):
            z = inner(columns[i], columns[j])
            K[i, j] = z
            K[j, i] = np.conjugate(z)
    return 0.5 * (K + K.conj().T)


def save_npz(path: Path, columns):
    states = sorted(set().union(*(set(c) for c in columns)), key=repr)
    if not states:
        raise RuntimeError("empty five-block P support")
    sidx = {s: i for i, s in enumerate(states)}
    nspin = len(states[0][0])
    nK = len(states[0][1])
    spins = np.empty((len(states), nspin), dtype=np.int16)
    Ks = np.empty((len(states), nK), dtype=np.int16)
    P = np.zeros((len(states), len(columns)), complex)
    for i, st in enumerate(states):
        spins[i] = np.asarray(st[0], dtype=np.int16)
        Ks[i] = np.asarray(st[1], dtype=np.int16)
    for j, col in enumerate(columns):
        for st, amp in col.items():
            P[sidx[st], j] = amp
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        basis_spins2=spins,
        basis_K2=Ks,
        P_real=P.real,
        P_imag=P.imag,
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return states, P, digest


def run(upstream_root: Path, expected_sha: str, npz_path: Path):
    upstream_sha = verify_upstream(upstream_root, expected_sha)
    M = load_upstream(upstream_root)
    ZVM, BE, PL, E = M["ZVM"], M["BE"], M["PL"], M["E"]

    ZVM.patch_and_clear()
    coarse = PL.seed_16cell_boundary()
    fine, parent = BE.barycentric_with_parent(coarse)
    D = PL.DualComplex(fine)
    G = E.LocalPLPeterWeylEuclidean(D)
    seed = ((1,) * len(G.EDGES), (0,) * D.n_tets)

    neighbours = face_neighbours(D, parent, 0)
    blocks = [0] + neighbours
    if len(neighbours) != 4 or len(set(blocks)) != 5:
        raise RuntimeError(("center parent does not have four face neighbours", neighbours))
    center_vertices = set(coarse[0])
    face_intersections = [len(center_vertices & set(coarse[p])) for p in neighbours]

    columns = []
    labels = []
    block_rows = []
    source_rows = []
    block_column_ranges = []
    for parent_id in blocks:
        nodes, source_cols, rows = build_parent_sources(D, G, parent, parent_id, seed)
        six, edge_labels, nu, defect = six_edge_lifts(coarse[parent_id], source_cols)
        start = len(columns)
        for edge_index, (col, edge) in enumerate(zip(six, edge_labels)):
            columns.append(col)
            labels.append({
                "column": len(columns) - 1,
                "parent": int(parent_id),
                "coarse_parent_tetra": [int(x) for x in coarse[parent_id]],
                "local_edge_index": edge_index,
                "coarse_edge_vertices": [int(x) for x in edge],
            })
        stop = len(columns)
        block_column_ranges.append((start, stop))
        block_rows.append({
            "parent": int(parent_id),
            "coarse_parent_tetra": [int(x) for x in coarse[parent_id]],
            "fine_chambers": len(nodes),
            "source_columns": len(source_cols),
            "P_columns": len(six),
            "microscopic_edge_lift_common_norm_square": nu,
            "microscopic_edge_Gram_relative_I6_defect": defect,
            "P_support_sizes": [len(c) for c in six],
        })
        source_rows.extend(rows)

    if len(columns) != 30:
        raise RuntimeError(f"expected 30 P columns, got {len(columns)}")

    K = full_gram(columns)
    evals = np.linalg.eigvalsh(K)
    maxeval = float(np.max(np.abs(evals)))
    rank_tol = max(1e-12, maxeval * 1e-10)
    rank = int(np.sum(evals > rank_tol))
    gram_I_defect = float(
        np.linalg.norm(K - np.eye(30)) / max(np.linalg.norm(K), 1e-300)
    )
    herm_defect = float(
        np.linalg.norm(K - K.conj().T) / max(np.linalg.norm(K), 1e-300)
    )

    cross_max = 0.0
    cross_support_collisions = 0
    for bi, (a0, a1) in enumerate(block_column_ranges):
        for bj, (b0, b1) in enumerate(block_column_ranges):
            if bj <= bi:
                continue
            cross_max = max(cross_max, float(np.max(np.abs(K[a0:a1, b0:b1]))))
            for i in range(a0, a1):
                for j in range(b0, b1):
                    cross_support_collisions += len(set(columns[i]) & set(columns[j]))

    states, P, npz_sha256 = save_npz(npz_path, columns)
    serialized_gram = P.conj().T @ P
    serialized_defect = float(
        np.linalg.norm(serialized_gram - K) / max(np.linalg.norm(K), 1e-300)
    )

    checks = {
        "frozen_upstream_sha_exact": upstream_sha == expected_sha,
        "L1_closed_nodes_384": D.n_tets == 384,
        "L1_dual_links_768": len(G.EDGES) == 768,
        "center_has_exactly_four_face_neighbours": len(neighbours) == 4,
        "each_neighbour_shares_one_coarse_face": face_intersections == [3, 3, 3, 3],
        "five_parent_blocks": len(blocks) == 5,
        "twenty_four_sources_per_block": all(r["source_columns"] == 24 for r in block_rows),
        "six_P_columns_per_block": all(r["P_columns"] == 6 for r in block_rows),
        "thirty_P_columns_total": len(columns) == 30,
        "one_strict_q4_spec_per_source": all(r["strict_specs"] == 1 for r in source_rows),
        "six_total_q4_specs_per_source": all(r["q4_specs_total"] == 6 for r in source_rows),
        "strict_boundary_links_unchanged": all(r["boundary_spins_unchanged"] for r in source_rows),
        "strict_exterior_exact_background": all(r["exterior_exactly_background_by_labels"] for r in source_rows),
        "full_global_labels_preserved": all(r["full_global_label_lengths"] for r in source_rows),
        "each_block_true_microscopic_I6": all(
            r["microscopic_edge_Gram_relative_I6_defect"] < GRAM_TOL for r in block_rows
        ),
        "cross_block_basis_support_disjoint": cross_support_collisions == 0,
        "cross_block_overlap_zero": cross_max < 1e-14,
        "global_P_Gram_Hermitian": herm_defect < 1e-14,
        "global_P_Gram_I30": gram_I_defect < GRAM_TOL,
        "global_P_rank_30": rank == 30,
        "serialized_P_reproduces_Gram": serialized_defect < 1e-14,
    }

    return {
        "status": "exact five-block strict-interior 30-column microscopic P carrier",
        "passed": bool(all(checks.values())),
        "science_status": "FIVEBLOCK_STRICT_P_CARRIER_PRECURSOR",
        "upstream_sha": upstream_sha,
        "checks": checks,
        "center_parent": 0,
        "face_neighbours": neighbours,
        "selected_parents": blocks,
        "center_neighbour_coarse_vertex_intersections": face_intersections,
        "block_rows": block_rows,
        "columns": labels,
        "source_rows": source_rows,
        "global_basis_states": len(states),
        "global_P_shape": list(P.shape),
        "P_column_support_sizes": [len(c) for c in columns],
        "P_Gram_eigenvalues": [float(x) for x in evals],
        "P_Gram_rank_tolerance": rank_tol,
        "P_Gram_rank": rank,
        "P_Gram_I30_relative_defect": gram_I_defect,
        "P_Gram_Hermiticity_relative_defect": herm_defect,
        "cross_block_max_abs_overlap": cross_max,
        "cross_block_basis_support_collisions": cross_support_collisions,
        "serialized_P_Gram_relative_defect": serialized_defect,
        "npz_path": str(npz_path),
        "npz_sha256": npz_sha256,
        "npz_arrays": {
            "basis_spins2": [len(states), len(G.EDGES)],
            "basis_K2": [len(states), D.n_tets],
            "P_real": [len(states), 30],
            "P_imag": [len(states), 30],
        },
        "construction_note": (
            "P uses only strict-interior q4 dynamics. Shared-face recoupling is deliberately absent because boundary links are unchanged and cross-parent supports are exactly disjoint in the full global Gauss basis."
        ),
        "next_step": (
            "Apply frozen H_E^sine and Hermitian S to these 30 actual columns; only the resulting boundary-touching images require complete six-link face recoupling before target-independent Q whitening."
        ),
        "scope_note": (
            "No GR/ADM/TT/Wilson target information is used. This is a spatial microscopic carrier precursor, not a physical spacetime TT kernel."
        ),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--upstream-root", type=Path, required=True)
    p.add_argument("--expected-upstream-sha", default=EXPECTED_UPSTREAM_SHA)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--npz", type=Path)
    a = p.parse_args()

    npz_path = a.npz or a.output.with_suffix(".npz")
    try:
        result = run(a.upstream_root.resolve(), a.expected_upstream_sha, npz_path)
    except Exception as exc:
        result = {
            "status": "exact five-block strict-interior 30-column microscopic P carrier",
            "passed": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(text, end="")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(text, encoding="utf-8")
    return 0 if result.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
