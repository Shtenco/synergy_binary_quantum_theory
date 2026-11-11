#!/usr/bin/env python3
"""Post-HDA master-constraint spectrum from the actual Peter-Weyl H_v columns.

This script does NOT test HDA. It consumes the already-existing regulated K5
Peter-Weyl Hamiltonian constraints as input and asks the next question: what is
the common physical zero/low-energy sector selected by

    M_G = sum_{v,w} H_v^dag G[v,w] H_w

on a declared finite carrier domain?

The default domain is a six-dimensional low-logical-excitation carrier inside
the all-j=1/2 Gauss/intertwiner sector: the all-K=0 state plus the five states
with one node recoupling label K=2. --domain full uses all 32 K=0/2 states.

A zero mode of the finite Galerkin master operator is a genuine combination in
the declared domain annihilated by every included H_v (up to tolerance). If no
zero mode exists, the script reports that result rather than tuning HDA or
renaming a small eigenvalue as zero.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_peter_weyl_safe_hda_column as PW


def inner_sparse(a: dict, b: dict) -> complex:
    if len(a) > len(b):
        return np.conjugate(inner_sparse(b, a))
    return sum(np.conjugate(x) * b.get(k, 0j) for k, x in a.items())


def carrier_basis(domain: str):
    spins = (1,) * len(PW.EDGES)
    if domain == "full":
        return PW.basis_full_jhalf()
    if domain != "low":
        raise ValueError(f"unknown domain {domain!r}")
    out = [(spins, (0, 0, 0, 0, 0))]
    for v in PW.VERT:
        ks = [0] * len(PW.VERT)
        ks[v] = 2
        out.append((spins, tuple(ks)))
    return out


def positive_metrics(n: int, seed: int):
    rng = np.random.default_rng(seed)
    I = np.eye(n, dtype=complex)
    D = np.diag(np.arange(1, n + 1, dtype=float)).astype(complex)
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    R = X.conj().T @ X + 0.5 * I
    R = 0.5 * (R + R.conj().T)
    return {"identity": I, "diagonal": D, "random_positive": R}


def constraint_columns(basis, nodes, jmax2: int, mode: str = "legacy"):
    if mode not in ("legacy", "no_prune"):
        raise ValueError(f"unsupported pruning mode: {mode}")
    columns = {}
    supports = {}
    for a, key in enumerate(basis):
        source = {key: 1.0 + 0.0j}
        for v in nodes:
            raw = PW.apply_H_cached_state(source, v, jmax2, mode=mode)
            col = PW.prune_state(raw, 1e-10, mode=mode)
            columns[(v, a)] = col
            supports[f"v{v}_b{a}"] = len(col)
    return columns, supports


def cross_gram(columns, nodes, dim: int):
    # B[v,w,a,b] = <H_v a | H_w b>.
    B = np.zeros((len(nodes), len(nodes), dim, dim), dtype=complex)
    for iv, v in enumerate(nodes):
        for iw, w in enumerate(nodes):
            for a in range(dim):
                av = columns[(v, a)]
                for b in range(dim):
                    B[iv, iw, a, b] = inner_sparse(av, columns[(w, b)])
    return B


def master_from_gram(B: np.ndarray, G: np.ndarray) -> np.ndarray:
    # Return the RAW assembly. Never symmetrize before measuring its defect.
    return np.einsum("vw,vwab->ab", G, B, optimize=True)


def hermitian_diagnostics(M: np.ndarray):
    defect = float(np.linalg.norm(M - M.conj().T, ord="fro"))
    scale = max(float(np.linalg.norm(M, ord="fro")), np.finfo(float).tiny)
    return defect, defect / scale


def classify_spectrum(evals: np.ndarray, zero_threshold: float, psd_tolerance: float):
    if not np.all(np.isfinite(evals)):
        raise ValueError("master spectrum contains non-finite eigenvalues")
    if zero_threshold < 0 or psd_tolerance < 0:
        raise ValueError("spectral tolerances must be nonnegative")
    minimum = float(np.min(evals)) if len(evals) else 0.0
    if minimum < -psd_tolerance:
        raise ValueError(
            f"master operator is not PSD within tolerance: min eigenvalue={minimum}, "
            f"allowed={psd_tolerance}"
        )
    # Negative roundoff is not silently promoted to a zero mode.
    return np.abs(evals) <= zero_threshold


def projector_from_eigh(evals: np.ndarray, U: np.ndarray, threshold: float):
    mask = np.abs(evals) <= threshold
    if not np.any(mask):
        return np.zeros((U.shape[0], U.shape[0]), dtype=complex), mask
    Q = U[:, mask]
    return Q @ Q.conj().T, mask


def run(domain="low", jmax2=5, seed=20260818, zero_rtol=1e-10,
        mode="legacy", hermitian_rtol=1e-10, psd_rtol=1e-10):
    if mode not in ("legacy", "no_prune"):
        raise ValueError(f"unsupported pruning mode: {mode}")
    if min(zero_rtol, hermitian_rtol, psd_rtol) < 0:
        raise ValueError("all tolerances must be nonnegative")
    nodes = tuple(PW.VERT)
    basis = carrier_basis(domain)
    columns, supports = constraint_columns(basis, nodes, jmax2, mode)
    B = cross_gram(columns, nodes, len(basis))

    metrics = positive_metrics(len(nodes), seed)
    rows = []
    projectors = []
    all_ok = True
    for name, G in metrics.items():
        M_raw = master_from_gram(B, G)
        herm_abs, herm_rel = hermitian_diagnostics(M_raw)
        # Spectral decomposition uses the Hermitian part only after preserving
        # and checking the raw assembly defect above.
        M = 0.5 * (M_raw + M_raw.conj().T)
        evals, U = np.linalg.eigh(M)
        scale = max(1.0, float(np.max(np.abs(evals))))
        threshold = zero_rtol * scale
        psd_tolerance = psd_rtol * scale
        mask = classify_spectrum(evals, threshold, psd_tolerance)
        P0, _ = projector_from_eigh(evals, U, threshold)
        rank0 = int(np.count_nonzero(mask))
        positive = evals[~mask]
        gap = float(np.min(positive)) if len(positive) else 0.0
        min_eval = float(np.min(evals))
        min_G = float(np.min(np.linalg.eigvalsh(G)))
        # Every selected zero vector must be annihilated by each actual H_v.
        annihilation = 0.0
        if rank0:
            Q0 = U[:, mask]
            for q in range(Q0.shape[1]):
                coeff = Q0[:, q]
                for v in nodes:
                    out = {}
                    for a, amp in enumerate(coeff):
                        if abs(amp) < 1e-14:
                            continue
                        for key, val in columns[(v, a)].items():
                            out[key] = out.get(key, 0j) + amp * val
                    annihilation = max(annihilation, math.sqrt(sum(abs(x) ** 2 for x in out.values())))
        good = (herm_rel <= hermitian_rtol and min_G > 0.0
                and min_eval >= -psd_tolerance and annihilation < 2e-7)
        all_ok &= good
        projectors.append((name, P0, rank0))
        rows.append({
            "metric": name,
            "G_min_eigenvalue": min_G,
            "master_eigenvalues": [float(x) for x in evals],
            "master_min_eigenvalue": min_eval,
            "zero_threshold": threshold,
            "zero_sector_rank": rank0,
            "first_positive_master_gap": gap,
            "zero_sector_constraint_annihilation_max_norm": annihilation,
            "raw_master_hermiticity_defect_frobenius": herm_abs,
            "raw_master_hermiticity_relative_defect": herm_rel,
            "psd_tolerance": psd_tolerance,
            "spectrum_classification": "abs(lambda)<=zero_threshold; PSD checked separately",
            "passed_operator_checks": bool(good),
        })

    pair_errors = []
    for i in range(len(projectors)):
        for j in range(i + 1, len(projectors)):
            ni, Pi, ri = projectors[i]
            nj, Pj, rj = projectors[j]
            err = float(np.linalg.norm(Pi - Pj)) if ri == rj else float("inf")
            pair_errors.append({"a": ni, "b": nj, "rank_a": ri, "rank_b": rj, "projector_difference": err})
            all_ok &= ri == rj and err < 2e-6

    zero_ranks = [r["zero_sector_rank"] for r in rows]
    found = bool(zero_ranks and min(zero_ranks) > 0)
    return {
        "status": "post-HDA Peter-Weyl master-constraint physical-sector pilot",
        "passed_internal_checks": bool(all_ok),
        "scientific_result_is_not_forced_to_pass": True,
        "domain": domain,
        "pruning_mode": mode,
        "tolerances": {"zero_rtol": zero_rtol, "hermitian_rtol": hermitian_rtol, "psd_rtol": psd_rtol},
        "domain_dimension": len(basis),
        "domain_definition": "all-K=0 plus five single K=2 node excitations" if domain == "low" else "all 32 all-j=1/2 K=0/2 intertwiner states",
        "included_constraint_nodes": list(nodes),
        "Jmax": jmax2 / 2,
        "basis_keys": [{"spins": list(k[0]), "Ks": list(k[1])} for k in basis],
        "constraint_column_supports": supports,
        "positive_constraint_metric_results": rows,
        "zero_projector_metric_independence": pair_errors,
        "finite_physical_zero_sector_found": found,
        "interpretation_if_false": "No exact common H_v kernel exists in this declared carrier at the chosen cutoff/tolerance. Do not retune HDA; enlarge/refine the physical carrier or use a preregistered spectral-window/rigging sequence.",
        "interpretation_if_true": "A finite common constraint kernel exists in the declared carrier. Its projector can be used for source-deformed metric matrix elements, but physical omega still requires a history/relational construction and refinement control.",
        "claim_boundary": "This is a finite Galerkin physical-sector diagnostic using the actual regulated Peter-Weyl H_v operators. It is not the continuum rigging map, not an external-time propagator, and not a blind physical prediction.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", choices=("low", "full"), default="low")
    ap.add_argument("--jmax2", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument("--zero-rtol", type=float, default=1e-10)
    ap.add_argument("--hermitian-rtol", type=float, default=1e-10)
    ap.add_argument("--psd-rtol", type=float, default=1e-10)
    ap.add_argument("--mode", choices=("legacy", "no_prune"), default="legacy")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.domain, args.jmax2, args.seed, args.zero_rtol,
              args.mode, args.hermitian_rtol, args.psd_rtol)
    if args.output is None:
        suffix = "NO_PRUNE" if args.mode == "no_prune" else "LEGACY"
        args.output = Path(f"verification_results/PETER_WEYL_MASTER_PROJECTOR_{suffix}.json")
    text = json.dumps(out, indent=2, ensure_ascii=False)
    print(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text + "\n", encoding="utf-8")
    return 0 if out["passed_internal_checks"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
