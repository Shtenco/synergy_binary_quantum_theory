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
    """Sparse complex inner product with compensated real/imaginary sums."""
    if len(a) > len(b):
        return np.conjugate(inner_sparse(b, a))
    products = (np.conjugate(x) * b.get(k, 0j) for k, x in a.items())
    values = list(products)
    return complex(
        math.fsum(z.real for z in values),
        math.fsum(z.imag for z in values),
    )


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
    if not np.all(np.isfinite(M)):
        raise ValueError("raw master matrix contains non-finite entries")
    defect = stable_complex_norm((M - M.conj().T).ravel())
    scale = stable_complex_norm(M.ravel())
    relative = defect / scale if scale > 0.0 else (0.0 if defect == 0.0 else float("inf"))
    return defect, relative


def spectral_thresholds(evals: np.ndarray, zero_rtol: float, psd_rtol: float):
    """Return scale-relative thresholds; no unit-sized floor is introduced."""
    if (not np.isfinite(zero_rtol) or not np.isfinite(psd_rtol)
            or zero_rtol < 0.0 or psd_rtol < 0.0):
        raise ValueError("spectral relative tolerances must be finite and nonnegative")
    if not np.all(np.isfinite(evals)):
        raise ValueError("master spectrum contains non-finite eigenvalues")
    scale = float(np.max(np.abs(evals))) if len(evals) else 0.0
    return scale, zero_rtol * scale, psd_rtol * scale


def classify_spectrum(evals: np.ndarray, zero_threshold: float, psd_tolerance: float):
    if not np.all(np.isfinite(evals)):
        raise ValueError("master spectrum contains non-finite eigenvalues")
    if (not np.isfinite(zero_threshold) or not np.isfinite(psd_tolerance)
            or zero_threshold < 0 or psd_tolerance < 0):
        raise ValueError("spectral tolerances must be finite and nonnegative")
    minimum = float(np.min(evals)) if len(evals) else 0.0
    if minimum < -psd_tolerance:
        raise ValueError(
            f"master operator is not PSD within tolerance: min eigenvalue={minimum}, "
            f"allowed={psd_tolerance}"
        )
    # The zero window is symmetric; PSD tolerance is checked independently
    # above, so a negative residual outside the zero window is never a zero mode.
    return np.abs(evals) <= zero_threshold


def projector_from_eigh(evals: np.ndarray, U: np.ndarray, threshold: float):
    mask = np.abs(evals) <= threshold
    if not np.any(mask):
        return np.zeros((U.shape[0], U.shape[0]), dtype=complex), mask
    Q = U[:, mask]
    return Q @ Q.conj().T, mask


def stable_complex_norm(values) -> float:
    """Scaled sum-of-squares norm for complex values; avoids naive squaring."""
    scale = 0.0
    ssq = 1.0
    for value in values:
        z = complex(value)
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise FloatingPointError("non-finite value in stable complex norm")
        for component in (abs(z.real), abs(z.imag)):
            if component == 0.0:
                continue
            if scale < component:
                ssq = 1.0 + ssq * (scale / component) ** 2
                scale = component
            else:
                ssq += (component / scale) ** 2
    return 0.0 if scale == 0.0 else scale * math.sqrt(ssq)


def constraint_annihilation_metrics(columns, coeff: np.ndarray, nodes):
    """Return absolute and Frobenius-relative residuals for every constraint."""
    coeff = np.asarray(coeff, dtype=np.complex128)
    if coeff.ndim != 1 or not np.all(np.isfinite(coeff)):
        raise ValueError("constraint coefficients must be a finite vector")
    psi_norm = stable_complex_norm(coeff)
    if psi_norm == 0.0:
        raise ValueError("zero input vector cannot certify constraint annihilation")

    per_node = {}
    for v in nodes:
        terms = {}
        for a, amp in enumerate(coeff):
            column = columns[(v, a)]
            if amp == 0j:
                continue
            for key, value in column.items():
                terms.setdefault(key, []).append(amp * value)
        residual_values = [
            complex(math.fsum(z.real for z in values),
                    math.fsum(z.imag for z in values))
            for values in terms.values()
        ]
        residual_norm = stable_complex_norm(residual_values)
        operator_frobenius = stable_complex_norm(
            value
            for a in range(len(coeff))
            for value in columns[(v, a)].values()
        )
        if operator_frobenius == 0.0:
            relative = 0.0 if residual_norm == 0.0 else float("inf")
        else:
            # Divide sequentially to avoid overflowing the product
            # ||H_v||_F * ||psi||_2.
            relative = (residual_norm / operator_frobenius) / psi_norm
        per_node[str(v)] = {
            "absolute_residual": residual_norm,
            "operator_frobenius_norm": operator_frobenius,
            "relative_residual": relative,
        }

    return {
        "max_absolute_residual": max(
            (item["absolute_residual"] for item in per_node.values()), default=0.0
        ),
        "max_relative_residual": max(
            (item["relative_residual"] for item in per_node.values()), default=0.0
        ),
        "input_vector_norm": psi_norm,
        "per_node": per_node,
    }


def constraint_annihilation_norm(columns, coeff: np.ndarray, nodes) -> float:
    """Backward-compatible absolute residual wrapper with stable norming."""
    coeff = np.asarray(coeff, dtype=np.complex128)
    if coeff.ndim == 1 and (coeff.size == 0 or np.all(coeff == 0j)):
        return 0.0
    return constraint_annihilation_metrics(columns, coeff, nodes)["max_absolute_residual"]


def volume_raw_hermiticity_diagnostics(
    registry=None, rtol=1e-10, expected_unique_calls=None, max_spin2=None
):
    """Fail-closed audit of raw volume defects and registry call coverage."""
    if not np.isfinite(rtol) or rtol < 0.0:
        raise ValueError("volume hermiticity tolerance must be finite and nonnegative")
    if expected_unique_calls is not None and (
        not isinstance(expected_unique_calls, int) or expected_unique_calls < 0
    ):
        raise ValueError("expected_unique_calls must be a nonnegative integer")
    if max_spin2 is not None and (
        not isinstance(max_spin2, int) or max_spin2 < 0
    ):
        raise ValueError("max_spin2 must be a nonnegative integer")
    if registry is None:
        registry = getattr(PW, "VOLUME_RAW_DEFECTS", {})
    records = []
    failures = []
    for key, raw in sorted(registry.items(), key=lambda item: repr(item[0])):
        key_valid = (
            isinstance(key, tuple) and len(key) == 3
            and all(
                isinstance(x, (int, np.integer)) and x >= 0
                and (max_spin2 is None or x <= max_spin2)
                for x in key
            )
        )
        try:
            defect = float(raw.get("relative_antihermitian_defect", float("nan")))
            raw_norm = float(raw.get("raw_norm", float("nan")))
            delta_frob = float(raw.get("symmetrization_delta_frobenius_norm", float("nan")))
            delta_2_bound = float(raw.get("symmetrization_delta_spectral_norm_upper_bound", float("nan")))
            relative_bound = float(raw.get("symmetrization_relative_frobenius_bound", float("nan")))
        except (AttributeError, TypeError, ValueError, OverflowError):
            defect = raw_norm = delta_frob = delta_2_bound = relative_bound = float("nan")
        finite = all(math.isfinite(x) for x in (
            defect, raw_norm, delta_frob, delta_2_bound, relative_bound
        ))
        record = {
            "spin_triple": [int(x) for x in key] if key_valid else None,
            "invalid_key_repr": None if key_valid else repr(key),
            "relative_antihermitian_defect": defect if math.isfinite(defect) else None,
            "raw_norm": raw_norm if math.isfinite(raw_norm) else None,
            "symmetrization_delta_frobenius_norm": delta_frob if math.isfinite(delta_frob) else None,
            "symmetrization_delta_spectral_norm_upper_bound": delta_2_bound if math.isfinite(delta_2_bound) else None,
            "symmetrization_relative_frobenius_bound": relative_bound if math.isfinite(relative_bound) else None,
        }
        records.append(record)
        expected_relative = (
            delta_frob / raw_norm if raw_norm > 0.0
            else (0.0 if delta_frob == 0.0 else float("inf"))
        )
        consistency_tol = 100.0 * np.finfo(float).eps * max(
            1.0, abs(expected_relative) if math.isfinite(expected_relative) else 1.0
        )
        defect_consistent = (
            math.isfinite(expected_relative)
            and abs(defect - 2.0 * relative_bound) <= consistency_tol
        )
        bounds_consistent = (
            finite and delta_frob >= 0.0 and delta_2_bound >= delta_frob
            and relative_bound >= 0.0
            and relative_bound <= 1.0 + consistency_tol
            and abs(relative_bound - expected_relative) <= consistency_tol
            and defect_consistent
        )
        if (not key_valid or not finite or defect < 0.0 or raw_norm < 0.0
                or delta_frob < 0.0 or delta_2_bound < 0.0
                or relative_bound < 0.0 or defect > rtol
                or not bounds_consistent):
            failures.append(record)
    coverage_complete = (
        expected_unique_calls is None or len(records) == expected_unique_calls
    )
    passed = bool(records) and not failures and coverage_complete
    finite_defects = [
        record["relative_antihermitian_defect"]
        for record in records
        if record["relative_antihermitian_defect"] is not None
    ]
    return {
        "passed": passed,
        "rtol": float(rtol),
        "triples_checked": len(records),
        "expected_unique_calls": expected_unique_calls,
        "spin_cutoff_doubled": max_spin2,
        "coverage_complete": coverage_complete,
        "coverage_shortfall": (
            None if expected_unique_calls is None
            else expected_unique_calls - len(records)
        ),
        "max_relative_defect": max(finite_defects) if finite_defects else None,
        "failures": failures,
        "records": records,
        "fail_closed_if_registry_empty": True,
        "fail_closed_if_coverage_incomplete": True,
    }


def run(domain="low", jmax2=5, seed=20260818, zero_rtol=1e-10,
        mode="legacy", hermitian_rtol=1e-10, psd_rtol=1e-10,
        annihilation_rtol=2e-7, volume_hermitian_rtol=1e-10):
    if mode not in ("legacy", "no_prune"):
        raise ValueError(f"unsupported pruning mode: {mode}")
    tolerances = (zero_rtol, hermitian_rtol, psd_rtol,
                  annihilation_rtol, volume_hermitian_rtol)
    if not all(np.isfinite(value) and value >= 0 for value in tolerances):
        raise ValueError("all tolerances must be finite and nonnegative")
    nodes = tuple(PW.VERT)
    basis = carrier_basis(domain)
    # Isolate this run's volume ledger. Every unique cache miss must execute
    # the producer body and create exactly one raw-defect record.
    PW.volume123_matrix.cache_clear()
    PW.VOLUME_RAW_DEFECTS.clear()
    columns, supports = constraint_columns(basis, nodes, jmax2, mode)
    volume_cache_misses = PW.volume123_matrix.cache_info().misses
    volume_diagnostics = volume_raw_hermiticity_diagnostics(
        rtol=volume_hermitian_rtol,
        expected_unique_calls=volume_cache_misses,
        max_spin2=jmax2,
    )
    volume_diagnostics["volume_cache_misses"] = volume_cache_misses
    volume_ok = volume_diagnostics["passed"]
    B = cross_gram(columns, nodes, len(basis))

    metrics = positive_metrics(len(nodes), seed)
    rows = []
    projectors = []
    all_ok = volume_ok
    for name, G in metrics.items():
        M_raw = master_from_gram(B, G)
        herm_abs, herm_rel = hermitian_diagnostics(M_raw)
        # Spectral decomposition uses the Hermitian part only after preserving
        # and checking the raw assembly defect above.
        M = 0.5 * (M_raw + M_raw.conj().T)
        evals, U = np.linalg.eigh(M)
        scale, threshold, psd_tolerance = spectral_thresholds(
            evals, zero_rtol, psd_rtol
        )
        mask = classify_spectrum(evals, threshold, psd_tolerance)
        P0, _ = projector_from_eigh(evals, U, threshold)
        rank0 = int(np.count_nonzero(mask))
        positive = evals[evals > threshold]
        negative_nonzero = evals[evals < -threshold]
        gap = float(np.min(positive)) if len(positive) else 0.0
        min_eval = float(np.min(evals))
        min_G = float(np.min(np.linalg.eigvalsh(G)))
        # Every selected zero vector must be annihilated by each actual H_v.
        annihilation = 0.0
        annihilation_relative = 0.0
        annihilation_by_node = {}
        if rank0:
            Q0 = U[:, mask]
            for q in range(Q0.shape[1]):
                residuals = constraint_annihilation_metrics(columns, Q0[:, q], nodes)
                annihilation = max(annihilation, residuals["max_absolute_residual"])
                annihilation_relative = max(
                    annihilation_relative, residuals["max_relative_residual"]
                )
                for vertex, item in residuals["per_node"].items():
                    previous = annihilation_by_node.get(vertex, {
                        "absolute_residual": 0.0, "relative_residual": 0.0
                    })
                    annihilation_by_node[vertex] = {
                        "absolute_residual": max(
                            previous["absolute_residual"], item["absolute_residual"]
                        ),
                        "relative_residual": max(
                            previous["relative_residual"], item["relative_residual"]
                        ),
                    }
        good = (volume_ok and herm_rel <= hermitian_rtol and min_G > 0.0
                and min_eval >= -psd_tolerance
                and annihilation_relative <= annihilation_rtol)
        all_ok &= good
        projectors.append((name, P0, rank0))
        rows.append({
            "metric": name,
            "G_min_eigenvalue": min_G,
            "master_eigenvalues": [float(x) for x in evals],
            "master_spectral_scale": scale,
            "master_min_eigenvalue": min_eval,
            "zero_threshold": threshold,
            "zero_sector_rank": rank0,
            "first_positive_master_gap": gap,
            "negative_nonzero_eigenvalues": [float(x) for x in negative_nonzero],
            "negative_nonzero_count": int(len(negative_nonzero)),
            "zero_sector_constraint_annihilation_max_norm": annihilation,
            "zero_sector_constraint_annihilation_max_relative_residual": annihilation_relative,
            "zero_sector_constraint_annihilation_by_node": annihilation_by_node,
            "annihilation_relative_tolerance": annihilation_rtol,
            "volume_raw_hermiticity_passed": volume_ok,
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
            err = float(np.linalg.norm(Pi - Pj)) if ri == rj else None
            pair_errors.append({
                "a": ni, "b": nj, "rank_a": ri, "rank_b": rj,
                "projector_difference": err,
                "rank_mismatch": bool(ri != rj),
            })
            all_ok &= (ri == rj and err is not None and err < 2e-6)

    zero_ranks = [r["zero_sector_rank"] for r in rows]
    found = bool(zero_ranks and min(zero_ranks) > 0)
    return {
        "status": "post-HDA Peter-Weyl master-constraint physical-sector pilot",
        "passed_internal_checks": bool(all_ok),
        "scientific_result_is_not_forced_to_pass": True,
        "domain": domain,
        "pruning_mode": mode,
        "tolerances": {
            "zero_rtol": zero_rtol,
            "hermitian_rtol": hermitian_rtol,
            "psd_rtol": psd_rtol,
            "annihilation_relative_rtol": annihilation_rtol,
            "volume_hermitian_rtol": volume_hermitian_rtol,
        },
        "volume_raw_hermiticity": volume_diagnostics,
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
    ap.add_argument("--annihilation-rtol", type=float, default=2e-7)
    ap.add_argument("--volume-hermitian-rtol", type=float, default=1e-10)
    ap.add_argument("--mode", choices=("legacy", "no_prune"), default="legacy")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.domain, args.jmax2, args.seed, args.zero_rtol,
              args.mode, args.hermitian_rtol, args.psd_rtol,
              args.annihilation_rtol, args.volume_hermitian_rtol)
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
