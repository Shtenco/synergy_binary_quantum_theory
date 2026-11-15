#!/usr/bin/env python3
"""Build the actual finite K5 BQG master operator from all 60 local constraints.

This is not the synthetic random-matrix selftest.  It assembles the repository's
own graph-changing Thiemann-type local constraint maps on the 32D fully-active
K5 domain, with their outputs represented in the 140D Gauss Hilbert space.
The resulting positive master operator is an exact finite-regulator object on
that declared domain.  It is not the full Lorentzian/continuum physical theory.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import k5_thiemann_constraint_gate as K5C


def stable_norm(x) -> float:
    a = np.asarray(x, dtype=np.complex128).reshape(-1)
    if not np.all(np.isfinite(a.real)) or not np.all(np.isfinite(a.imag)):
        raise ValueError("non-finite vector/matrix entry")
    return float(math.hypot(*(abs(complex(z)) for z in a)))


def assemble_constraint_family():
    meta, states, reverse, Up, Ud, powers, full_order, V5, _target, trans = K5C.build_runtime()
    local = K5C.make_actions(Up, Ud, powers, trans)
    specs = K5C.local_specs()
    if len(meta) != 140 or len(full_order) != 32 or len(specs) != 60:
        raise AssertionError(
            f"unexpected K5 dimensions: Gauss={len(meta)}, active={len(full_order)}, constraints={len(specs)}"
        )

    # C[a] maps the 32D active domain into the complete 140D Gauss codomain.
    C = np.zeros((len(specs), 140, len(full_order)), dtype=np.complex128)
    raw_leakage = []
    for a, spec in enumerate(specs):
        max_leak = 0.0
        for col, basis_index in enumerate(full_order):
            raw = local(states[basis_index], *spec)
            raw_n = K5C.raw_norm(raw)
            projected = K5C.project_state(raw, reverse)
            projected_n = stable_norm(projected)
            if raw_n > 0.0:
                # Orthogonal projection cannot increase norm; tolerate only
                # roundoff-sized overshoot and fail on a genuine violation.
                ratio = projected_n / raw_n
                if ratio > 1.0 + 1e-9:
                    raise AssertionError(
                        f"Gauss projection increased norm for constraint {a}, column {col}: {ratio}"
                    )
                leak = max(0.0, 1.0 - min(1.0, ratio) ** 2)
                max_leak = max(max_leak, leak)
            elif projected_n != 0.0:
                raise AssertionError("zero raw constraint output has nonzero Gauss projection")
            C[a, :, col] = projected
        raw_leakage.append(max_leak)

    target = np.asarray(V5, dtype=np.complex128).reshape(-1)
    target_norm = stable_norm(target)
    if target_norm == 0.0:
        raise AssertionError("independent K5 target vector is zero")
    target = target / target_norm
    raw_target = {}
    for coefficient, basis_index in zip(target, full_order):
        for raw_index, amplitude in states[basis_index].items():
            value = raw_target.get(raw_index, 0j) + coefficient * amplitude
            if value != 0j:
                raw_target[raw_index] = value
    raw_target_residuals = [
        K5C.raw_norm(local(raw_target, *spec)) for spec in specs
    ]
    return C, specs, target, {
        "gauss_hilbert_dimension": len(meta),
        "active_domain_dimension": len(full_order),
        "constraint_count": len(specs),
        "max_gauss_projection_loss_fraction": max(raw_leakage, default=0.0),
        "mean_constraint_max_projection_loss_fraction": float(np.mean(raw_leakage)),
        "constraint_projection_loss_by_operator": raw_leakage,
        "independent_target_max_raw_constraint_residual": float(max(raw_target_residuals, default=0.0)),
        "independent_target_mean_raw_constraint_residual": float(np.mean(raw_target_residuals)),
    }


def master_from_constraints(C: np.ndarray, G: np.ndarray) -> np.ndarray:
    """M_ij = sum_ab,o conj(C[a,o,i]) G[a,b] C[b,o,j]."""
    if C.ndim != 3:
        raise ValueError("constraint tensor must have shape (m, codomain, domain)")
    m, _, n = C.shape
    if G.shape != (m, m):
        raise ValueError("constraint metric dimension mismatch")
    if not np.allclose(G, G.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("constraint metric must be Hermitian")
    ge = np.linalg.eigvalsh(G)
    if float(np.min(ge)) <= 0.0:
        raise ValueError("constraint metric must be strictly positive definite")
    GC = np.einsum("ab,boi->aoi", G, C, optimize=True)
    M = np.einsum("aoi,aoj->ij", C.conj(), GC, optimize=True)
    return 0.5 * (M + M.conj().T)


def physical_inner_product(P: np.ndarray, phi: np.ndarray, psi: np.ndarray) -> complex:
    """Finite rigging/projector inner product <eta(phi),eta(psi)>=<phi|P|psi>."""
    return complex(np.vdot(phi, P @ psi))


def physical_gram_check(P: np.ndarray, rng: np.random.Generator) -> dict:
    n = P.shape[0]
    probes = rng.normal(size=(n, 6)) + 1j * rng.normal(size=(n, 6))
    projected = P @ probes
    raw_gram = probes.conj().T @ projected
    hermiticity = stable_norm(raw_gram - raw_gram.conj().T)
    gram = 0.5 * (raw_gram + raw_gram.conj().T)
    values = np.linalg.eigvalsh(gram)
    scale = max(float(np.max(np.abs(values))), 1.0)
    rank = int(np.sum(values > 1e-9 * scale))
    return {
        "gram_dimension": int(gram.shape[0]),
        "gram_min_eigenvalue": float(np.min(values)),
        "gram_rank_after_kinematical_null_quotient": rank,
        "gram_hermiticity_residual": hermiticity,
        "positive_semidefinite": bool(float(np.min(values)) >= -1e-10 * scale),
    }


def spectral_projector(M: np.ndarray, *, relative_zero_tol: float = 1e-10):
    if not math.isfinite(relative_zero_tol) or relative_zero_tol <= 0.0:
        raise ValueError("relative_zero_tol must be finite and positive")
    evals, evecs = np.linalg.eigh(M)
    scale = max(float(np.max(np.abs(evals))), 1.0)
    if float(np.min(evals)) < -1e-10 * scale:
        raise AssertionError(f"master operator is not positive: min eigenvalue={np.min(evals)}")
    zero = np.abs(evals) <= relative_zero_tol * scale
    if not np.any(~zero):
        raise AssertionError("master operator has no positive spectral sector; gap undefined")
    P = evecs[:, zero] @ evecs[:, zero].conj().T
    gap = float(np.min(evals[~zero]))
    return evals, P, gap, int(np.sum(zero)), scale


def run(seed: int = 1729, relative_zero_tol: float = 1e-10) -> dict:
    C, specs, target, dimensions = assemble_constraint_family()
    m = C.shape[0]
    rng = np.random.default_rng(seed)

    # Include identity and dense, independently generated positive metrics.
    metrics = [np.eye(m, dtype=np.complex128)]
    for index in range(3):
        X = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        metrics.append(X.conj().T @ X / m + (0.25 + 0.1 * index) * np.eye(m))

    rows = []
    projectors = []
    target_residuals = []
    all_passed = True
    for index, G in enumerate(metrics):
        M = master_from_constraints(C, G)
        evals, P, gap, zero_rank, scale = spectral_projector(
            M, relative_zero_tol=relative_zero_tol
        )
        projector_idempotence = stable_norm(P @ P - P)
        projector_hermiticity = stable_norm(P - P.conj().T)
        constraint_residual = max(stable_norm(C[a] @ P) for a in range(m))
        target_residual = max(stable_norm(C[a] @ target) for a in range(m))
        target_projector_error = stable_norm(P @ target - target)
        physical_gram = physical_gram_check(P, rng)
        min_eig = float(np.min(evals))
        positivity_margin = min_eig / scale

        # Heat-kernel approximation on the actual master spectrum.
        heat_rows = []
        Ueval, Uvec = np.linalg.eigh(M)
        for multiplier in (0.5, 1.0, 2.0):
            T = multiplier / gap
            heat = (Uvec * np.exp(-T * Ueval)) @ Uvec.conj().T
            observed = float(np.linalg.norm(heat - P, 2))
            predicted = math.exp(-multiplier)
            heat_rows.append({
                "T_times_gap": multiplier,
                "operator_norm_heat_to_projector": observed,
                "spectral_gap_prediction": predicted,
                "absolute_error": abs(observed - predicted),
            })

        passed = (
            zero_rank == 1
            and projector_idempotence < 1e-8
            and projector_hermiticity < 1e-8
            and constraint_residual < 1e-8
            and target_residual < 1e-8
            and target_projector_error < 1e-7
            and physical_gram["positive_semidefinite"]
            and physical_gram["gram_rank_after_kinematical_null_quotient"] == zero_rank
            and physical_gram["gram_hermiticity_residual"] < 1e-8
            and dimensions["max_gauss_projection_loss_fraction"] < 1e-10
            and dimensions["independent_target_max_raw_constraint_residual"] < 1e-10
            and min_eig >= -1e-10 * scale
            and gap > 0.0
            and max(row["absolute_error"] for row in heat_rows) < 1e-7
        )
        all_passed = all_passed and passed
        projectors.append(P)
        target_residuals.append(target_residual)
        rows.append({
            "metric_index": index,
            "metric_min_eigenvalue": float(np.min(np.linalg.eigvalsh(G))),
            "master_min_eigenvalue": min_eig,
            "master_max_eigenvalue": float(np.max(evals)),
            "master_positive_gap": gap,
            "zero_sector_rank": zero_rank,
            "master_positivity_margin": positivity_margin,
            "projector_idempotence_residual": projector_idempotence,
            "projector_hermiticity_residual": projector_hermiticity,
            "all_constraint_projector_residual": constraint_residual,
            "independent_target_constraint_residual": target_residual,
            "target_to_projector_residual": target_projector_error,
            "finite_physical_inner_product": {
                "definition": "<eta(phi),eta(psi)>_phys = <phi|P_phys|psi>",
                "status": "positive semidefinite on the kinematical domain; positive definite after quotient by ker(P_phys)",
                **physical_gram,
            },
            "heat_kernel_checks": heat_rows,
            "passed": bool(passed),
        })

    pairwise = max(
        (stable_norm(projectors[i] - projectors[j])
         for i in range(len(projectors)) for j in range(i + 1, len(projectors))),
        default=0.0,
    )
    all_passed = all_passed and pairwise < 1e-7
    return {
        "schema": "bqg-k5-actual-constraint-master-projector-v1",
        "status": "PASS_K5_ACTUAL_FINITE_MASTER_PROJECTOR" if all_passed else "FAIL_K5_ACTUAL_FINITE_MASTER_PROJECTOR",
        "passed": bool(all_passed),
        "seed": seed,
        "constraint_definition": "T(v;a,b|c)=Tr[(U(vabv)-U(vbav)) U(v,c) [U(c,v),P_{k_v=4}]]",
        "constraint_family": [
            {"index": i, "spec": list(spec)} for i, spec in enumerate(specs)
        ],
        "dimensions": dimensions,
        "domain": "32D fully-active j=1/2 K5 sector",
        "codomain": "140D Gauss-invariant K5 Hilbert",
        "master_definition": "M_G = sum_ab C_a^dagger G_ab C_b, with all 60 actual local graph-changing constraint maps C_a",
        "theorem_scope": "Exact finite positive-master construction for the declared K5 domain and actual implemented Euclidean Thiemann-type local constraints; not the full Lorentzian Hamiltonian, arbitrary graph, continuum rigging map, or physical Born measure.",
        "max_pairwise_zero_projector_difference_across_positive_metrics": pairwise,
        "metric_results": rows,
        "all_metrics_same_zero_projector": bool(pairwise < 1e-7),
        "all_constraints_kill_projector": all(row["all_constraint_projector_residual"] < 1e-8 for row in rows),
        "finite_positive_physical_inner_product_passed": all(
            row["finite_physical_inner_product"]["positive_semidefinite"]
            and row["finite_physical_inner_product"]["gram_rank_after_kinematical_null_quotient"] == row["zero_sector_rank"]
            for row in rows
        ),
        "independent_target_residual_max": max(target_residuals, default=0.0),
        "independent_target_raw_residual_max": dimensions["independent_target_max_raw_constraint_residual"],
        "gauss_projection_loss_gate_passed": dimensions["max_gauss_projection_loss_fraction"] < 1e-10,
        "passed": bool(all_passed),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--relative-zero-tol", type=float, default=1e-10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.seed, args.relative_zero_tol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
