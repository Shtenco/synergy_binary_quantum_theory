#!/usr/bin/env python3
"""Assemble a finite graph-changing K5 constraint shell and its full pullback master.

For a finite input shell P and the complete regulated output columns H_v P,
the code keeps both the compressed action C_v=P H_v P and the coherent
outside-shell map L_v=(I-P)H_vP. It therefore reports

    M_projected = sum_v C_v^dagger C_v
    M_full      = M_projected + sum_v L_v^dagger L_v.

The second matrix is the exact pullback of the full regulated H_v actions to
the declared input shell, relative to the supplied orthonormal Peter-Weyl keys.
No amplitude threshold is applied in the default no_prune mode. This remains a
finite-spin, finite-input-shell diagnostic, not a continuum physical projector.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as K5
from bqg_sparse_leakage_gram import assemble_leakage_gram, classify_psd_spectrum


def graph_signature(key):
    """Abstract labelled-edge signature after cylindrical removal of j=0 links."""
    spins, _Ks = key
    return tuple(K5.EDGES[i] for i, s2 in enumerate(spins) if s2 != 0)


def apply_H_to_key(key, v: int, jmax2: int, mode: str = "no_prune"):
    """Apply the production node constraint without adding a second prune layer."""
    if mode not in ("legacy", "no_prune"):
        raise ValueError("mode must be 'legacy' or 'no_prune'")
    out = K5.apply_H_cached_state(
        {key: 1.0 + 0.0j}, v, jmax2, mode=mode
    )
    return dict(out)


def generate_shell(depth: int, jmax2: int, mode: str = "no_prune"):
    if not isinstance(depth, int) or isinstance(depth, bool) or depth < 0:
        raise ValueError("depth must be a nonnegative integer")
    if not isinstance(jmax2, int) or isinstance(jmax2, bool) or jmax2 < 1:
        raise ValueError("jmax2 must be a positive integer")
    initial = K5.basis_full_jhalf()[0]
    levels = [{initial}]
    seen = {initial}
    frontier = {initial}
    for _ in range(depth):
        nxt = set()
        for key in frontier:
            for v in K5.VERT:
                nxt.update(apply_H_to_key(key, v, jmax2, mode).keys())
        nxt -= seen
        levels.append(nxt)
        seen |= nxt
        frontier = nxt
        if not frontier:
            break
    basis = sorted(seen, key=repr)
    return initial, levels, basis


def compressed_constraints(basis, jmax2: int, mode: str = "no_prune"):
    """Return C_v and coherent L_v^dagger L_v for every node.

    Each output column is first fully and coherently assembled by the K5
    producer. The leakage helper then groups equal canonical outside keys
    across *different input columns*, retaining off-diagonal interference.
    """
    mats = []
    leakage_grams = []
    leakage = []
    basis_index = {key: i for i, key in enumerate(basis)}
    if len(basis_index) != len(basis):
        raise ValueError("input shell basis contains duplicate keys")

    for v in K5.VERT:
        output_columns = [
            apply_H_to_key(key, v, jmax2, mode) for key in basis
        ]
        result = assemble_leakage_gram(basis, output_columns)
        C = result.shell_action
        G_leak = result.leakage_gram
        mats.append(C)
        leakage_grams.append(G_leak)

        per_column = []
        for column in output_columns:
            kept2 = math.fsum(
                abs(amp) ** 2 for key, amp in column.items()
                if key in basis_index
            )
            outside2 = math.fsum(
                abs(amp) ** 2 for key, amp in column.items()
                if key not in basis_index
            )
            per_column.append({
                "kept_norm2": kept2,
                "outside_shell_norm2": outside2,
                "outside_nonzero_count": sum(
                    1 for key, amp in column.items()
                    if key not in basis_index and amp != 0j
                ),
            })

        total_kept2 = math.fsum(x["kept_norm2"] for x in per_column)
        total_outside2 = math.fsum(x["outside_shell_norm2"] for x in per_column)
        leakage.append({
            "vertex": v,
            "kept_frobenius_norm2": total_kept2,
            "outside_shell_frobenius_norm2": total_outside2,
            "outside_fraction": (
                total_outside2 / (total_kept2 + total_outside2)
                if total_kept2 + total_outside2 else 0.0
            ),
            "per_column": per_column,
            "coherent_leakage_gram_frobenius_norm": float(
                np.linalg.norm(G_leak, ord="fro")
            ),
            "coherent_offdiagonal_leakage_norm": result.diagnostics[
                "coherent_offdiagonal_leakage_norm"
            ],
            "leakage_hermiticity_defect": result.diagnostics[
                "leakage_hermiticity_defect"
            ],
            "leakage_min_eigenvalue": result.diagnostics[
                "leakage_min_eigenvalue"
            ],
            "outside_nonzero_column_entries": result.diagnostics[
                "outside_nonzero_column_entries"
            ],
            "pruning_performed_here": False,
            "exactness_condition": result.diagnostics["exactness_condition"],
        })
    return mats, leakage_grams, leakage


def relative_hermiticity_diagnostics(matrix):
    A = np.asarray(matrix, dtype=np.complex128)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("matrix must be square")
    if not (np.all(np.isfinite(A.real)) and np.all(np.isfinite(A.imag))):
        raise ValueError("matrix contains non-finite entries")
    defect = float(np.linalg.norm(A - A.conj().T, ord="fro"))
    scale = float(np.linalg.norm(A, ord="fro"))
    relative = defect / scale if scale else (0.0 if defect == 0.0 else float("inf"))
    return {
        "absolute_frobenius_defect": defect,
        "matrix_frobenius_scale": scale,
        "relative_frobenius_defect": relative,
    }


def master_matrices(constraints, leakage_grams):
    """Assemble raw projected/full masters and measure defects before cleanup."""
    if not constraints or len(constraints) != len(leakage_grams):
        raise ValueError("one leakage Gram is required for every node constraint")
    n = np.asarray(constraints[0]).shape[1]
    projected_raw = np.zeros((n, n), dtype=np.complex128)
    leakage_raw = np.zeros((n, n), dtype=np.complex128)
    for C, G in zip(constraints, leakage_grams):
        C = np.asarray(C, dtype=np.complex128)
        G = np.asarray(G, dtype=np.complex128)
        if C.shape != (n, n) or G.shape != (n, n):
            raise ValueError("constraint and leakage matrices must match shell dimension")
        if not (np.all(np.isfinite(C.real)) and np.all(np.isfinite(C.imag))):
            raise ValueError("constraint matrix contains non-finite entries")
        if not (np.all(np.isfinite(G.real)) and np.all(np.isfinite(G.imag))):
            raise ValueError("leakage Gram contains non-finite entries")
        projected_raw += C.conj().T @ C
        leakage_raw += G
    full_raw = projected_raw + leakage_raw

    # Diagnostics are intentionally computed on the raw assembled matrices.
    diagnostics = {
        "projected_master_raw_hermiticity": relative_hermiticity_diagnostics(projected_raw),
        "leakage_master_raw_hermiticity": relative_hermiticity_diagnostics(leakage_raw),
        "full_master_raw_hermiticity": relative_hermiticity_diagnostics(full_raw),
        "full_minus_projected_identity_frobenius": float(
            np.linalg.norm(full_raw - projected_raw - leakage_raw, ord="fro")
        ),
    }
    # Only now form Hermitian representatives for Hermitian eigensolvers.
    projected = 0.5 * projected_raw + 0.5 * projected_raw.conj().T
    leakage = 0.5 * leakage_raw + 0.5 * leakage_raw.conj().T
    full = 0.5 * full_raw + 0.5 * full_raw.conj().T
    return {
        "projected_raw": projected_raw,
        "leakage_raw": leakage_raw,
        "full_raw": full_raw,
        "projected": projected,
        "leakage": leakage,
        "full": full,
        "diagnostics": diagnostics,
    }


def projector_from_master(M, reltol=1e-10):
    if (not isinstance(reltol, (int, float, np.integer, np.floating))
            or isinstance(reltol, (bool, np.bool_))
            or not math.isfinite(float(reltol)) or reltol < 0.0):
        raise ValueError("reltol must be finite and nonnegative")
    result = classify_psd_spectrum(np.asarray(M, dtype=np.complex128), rel_tol=float(reltol))
    ev = result["eigenvalues"]
    threshold = result["relative_zero_threshold"]
    positive = ev[ev > threshold]
    return {
        "eigenvalues": ev,
        "zero_projector": result["zero_projector"],
        "zero_rank": result["nullity"],
        "positive_gap": float(np.min(positive)) if len(positive) else None,
        "spectral_scale": result["spectral_scale"],
        "zero_threshold": threshold,
        "psd_pass": result["psd_pass"],
        "significant_negative_eigenvalue": result["significant_negative_eigenvalue"],
        "hermiticity_defect": result["hermiticity_defect"],
    }


def graph_blocks(basis, P0):
    groups = defaultdict(list)
    for i, key in enumerate(basis):
        groups[graph_signature(key)].append(i)
    sigs = list(groups)
    rows = []
    for ga in sigs:
        ia = groups[ga]
        for gb in sigs:
            ib = groups[gb]
            block = P0[np.ix_(ib, ia)]
            fro = float(np.linalg.norm(block, ord="fro"))
            if fro > 1e-12:
                rows.append({
                    "source_active_edge_count": len(ga),
                    "target_active_edge_count": len(gb),
                    "source_sector_dim": len(ia),
                    "target_sector_dim": len(ib),
                    "frobenius_norm": fro,
                    "changes_graph": ga != gb,
                    "source_signature": [list(e) for e in ga],
                    "target_signature": [list(e) for e in gb],
                })
    rows.sort(key=lambda r: r["frobenius_norm"], reverse=True)
    return groups, rows


def run(depth: int = 1, jmax2: int = 5, mode: str = "no_prune",
        zero_rtol: float = 1e-10):
    initial, levels, basis = generate_shell(depth, jmax2, mode)
    constraints, leakage_grams, leakage = compressed_constraints(basis, jmax2, mode)
    masters = master_matrices(constraints, leakage_grams)
    full_spec = projector_from_master(masters["full"], zero_rtol)
    projected_spec = projector_from_master(masters["projected"], zero_rtol)
    groups, blocks = graph_blocks(basis, full_spec["zero_projector"])
    _projected_groups, projected_blocks = graph_blocks(
        basis, projected_spec["zero_projector"]
    )

    init_sig = graph_signature(initial)
    changed_basis = sum(1 for key in basis if graph_signature(key) != init_sig)
    changed_blocks = [item for item in blocks if item["changes_graph"]]
    raw_herm = masters["diagnostics"]
    herm_ok = all(
        item["relative_frobenius_defect"] <= 1e-10
        for item in (
            raw_herm["projected_master_raw_hermiticity"],
            raw_herm["leakage_master_raw_hermiticity"],
            raw_herm["full_master_raw_hermiticity"],
        )
    )
    decomposition_scale = float(np.linalg.norm(masters["full_raw"], ord="fro"))
    decomposition_tolerance = (
        64 * np.finfo(float).eps * max(decomposition_scale, np.finfo(float).tiny)
    )
    decomposition_ok = (
        raw_herm["full_minus_projected_identity_frobenius"] <= decomposition_tolerance
    )
    leakage_psd = [
        classify_psd_spectrum(G, rel_tol=zero_rtol)
        for G in leakage_grams
    ]
    leakage_psd_ok = all(item["psd_pass"] for item in leakage_psd)
    finite_ok = all(
        np.all(np.isfinite(matrix.real)) and np.all(np.isfinite(matrix.imag))
        for matrix in (masters["projected_raw"], masters["leakage_raw"], masters["full_raw"])
    )
    structural_ok = bool(
        finite_ok and herm_ok and decomposition_ok and leakage_psd_ok
        and full_spec["psd_pass"] and projected_spec["psd_pass"]
    )

    return {
        "status": "K5 Peter-Weyl graph-sector full-pullback master-shell diagnostic",
        "depth": depth,
        "Jmax": jmax2 / 2,
        "pruning_mode": mode,
        "basis_dimension": len(basis),
        "level_sizes": [len(level) for level in levels],
        "graph_sector_count": len(groups),
        "basis_states_outside_initial_graph": changed_basis,
        "constraint_shell_leakage": leakage,
        "projected_master_raw_hermiticity": raw_herm[
            "projected_master_raw_hermiticity"
        ],
        "leakage_master_raw_hermiticity": raw_herm[
            "leakage_master_raw_hermiticity"
        ],
        "full_master_raw_hermiticity": raw_herm[
            "full_master_raw_hermiticity"
        ],
        "full_minus_projected_identity_frobenius": raw_herm[
            "full_minus_projected_identity_frobenius"
        ],
        "projected_master_zero_rank": projected_spec["zero_rank"],
        "projected_master_gap": projected_spec["positive_gap"],
        "full_master_zero_rank": full_spec["zero_rank"],
        "full_master_gap_on_shell_input": full_spec["positive_gap"],
        "full_master_min_eigenvalue": float(np.min(full_spec["eigenvalues"])),
        "full_master_spectral_scale": full_spec["spectral_scale"],
        "full_master_zero_threshold": full_spec["zero_threshold"],
        "full_master_psd_pass": full_spec["psd_pass"],
        "projected_master_psd_pass": projected_spec["psd_pass"],
        "per_node_leakage_gram_psd_pass": [item["psd_pass"] for item in leakage_psd],
        "per_node_leakage_gram_min_eigenvalue": [
            float(np.min(item["eigenvalues"])) for item in leakage_psd
        ],
        "assembly_identity_tolerance": decomposition_tolerance,
        "projected_zero_projector_graph_blocks": projected_blocks,
        "full_zero_projector_graph_blocks": blocks,
        "full_nonzero_cross_graph_block_count": len(changed_blocks),
        "passed_structural_checks": structural_ok,
        "calculation_completed": True,
        "scientific_acceptance": None,
        "scope": (
            "No-prune or explicitly selected legacy finite-spin K5/Peter-Weyl H_v "
            "columns on a finite input shell. The full pullback includes the coherent "
            "outside-shell leakage Gram. It is not a continuum/rigging-map or physical "
            "projector certificate; cutoff, shell-depth and operator-domain convergence "
            "remain open."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--depth", type=int, default=1)
    ap.add_argument("--jmax2", type=int, default=5)
    ap.add_argument("--mode", choices=("no_prune", "legacy"), default="no_prune")
    ap.add_argument("--zero-rtol", type=float, default=1e-10)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.depth, args.jmax2, args.mode, args.zero_rtol)
    text = json.dumps(out, indent=2, ensure_ascii=False)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed_structural_checks"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
