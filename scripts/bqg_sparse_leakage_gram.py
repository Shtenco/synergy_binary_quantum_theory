#!/usr/bin/env python3
"""Coherent sparse leakage-Gram assembly for finite-shell BQG operators.

Input contract: columns[j] is the COMPLETE, coherently aggregated sparse output
map H e_j -> {canonical_basis_key: complex_amplitude}. No pruning is performed
here. Any upstream pruning must be accounted for separately and prevents an
'exact' label unless its residual is certified.

This module is a reusable algebraic kernel called by the K5 graph-shell
wrapper. Its output is exact only relative to complete, coherently aggregated
columns supplied by the upstream producer.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Mapping, Sequence
import math
import numpy as np


@dataclass(frozen=True)
class LeakageGramResult:
    shell_action: np.ndarray
    leakage_gram: np.ndarray
    pullback_gram: np.ndarray
    shell_keys: tuple[Hashable, ...]
    diagnostics: dict


def assemble_leakage_gram(
    shell_keys: Sequence[Hashable],
    output_columns: Sequence[Mapping[Hashable, complex]],
) -> LeakageGramResult:
    """Build C, G_leak=L†L, and C†C+G_leak without dropping amplitudes.

    Repeated contributions to a common outside key must already have been
    coherently summed in each input mapping. Contributions from different
    columns are then paired by the common canonical output key, retaining all
    off-diagonal interference terms.
    """
    keys = tuple(shell_keys)
    n = len(keys)
    if len(set(keys)) != n:
        raise ValueError("shell_keys must be unique canonical basis keys")
    if len(output_columns) != n:
        raise ValueError("one complete output mapping is required per shell column")
    shell_index = {key: i for i, key in enumerate(keys)}
    C = np.zeros((n, n), dtype=np.complex128)
    outside_rows: dict[Hashable, dict[int, complex]] = {}
    amplitude_magnitudes = []
    for j, out in enumerate(output_columns):
        for key, raw_amp in out.items():
            amp = complex(raw_amp)
            if not (np.isfinite(amp.real) and np.isfinite(amp.imag)):
                raise ValueError(f"non-finite amplitude in column {j}, key {key!r}")
            if amp == 0:
                continue
            amp_abs = float(math.hypot(amp.real, amp.imag))
            if not math.isfinite(amp_abs):
                raise FloatingPointError("amplitude magnitude exceeds float64 range")
            amplitude_magnitudes.append(amp_abs)
            i = shell_index.get(key)
            if i is not None:
                C[i, j] += amp
            else:
                row = outside_rows.setdefault(key, {})
                row[j] = row.get(j, 0j) + amp

    if not (np.all(np.isfinite(C.real)) and np.all(np.isfinite(C.imag))):
        raise FloatingPointError("shell action overflowed during coherent assembly")

    G_leak = np.zeros((n, n), dtype=np.complex128)
    for row in outside_rows.values():
        cols = sorted(row)
        vals = np.asarray([row[j] for j in cols], dtype=np.complex128)
        G_leak[np.ix_(cols, cols)] += np.outer(vals.conj(), vals)
    if not (np.all(np.isfinite(G_leak.real)) and np.all(np.isfinite(G_leak.imag))):
        raise FloatingPointError(
            "coherent leakage Gram is not representable in float64"
        )

    shell_gram = C.conj().T @ C
    if not (np.all(np.isfinite(shell_gram.real)) and np.all(np.isfinite(shell_gram.imag))):
        raise FloatingPointError("shell Gram is not representable in float64")
    pullback = shell_gram + G_leak
    if not (np.all(np.isfinite(pullback.real)) and np.all(np.isfinite(pullback.imag))):
        raise FloatingPointError("full pullback Gram is not representable in float64")
    leak_herm = float(np.linalg.norm(G_leak - G_leak.conj().T))
    pull_herm = float(np.linalg.norm(pullback - pullback.conj().T))
    leak_ev = np.linalg.eigvalsh(0.5 * G_leak + 0.5 * G_leak.conj().T) if n else np.array([])
    scale = float(np.linalg.norm(pullback, ord=2)) if n else 0.0
    if not all(np.isfinite(x) for x in (leak_herm, pull_herm, scale)):
        raise FloatingPointError("leakage/pullback diagnostics exceed float64 range")
    input_l1 = math.fsum(amplitude_magnitudes)
    squared_magnitudes = [x * x for x in amplitude_magnitudes]
    input_l2_sq = (
        math.fsum(squared_magnitudes)
        if all(math.isfinite(x) for x in squared_magnitudes)
        else None
    )
    if not math.isfinite(input_l1):
        input_l1 = None
    if input_l2_sq is not None and not math.isfinite(input_l2_sq):
        input_l2_sq = None
    diagnostics = {
        "shell_dimension": n,
        "distinct_outside_keys": len(outside_rows),
        "outside_nonzero_column_entries": sum(len(row) for row in outside_rows.values()),
        "coherent_offdiagonal_leakage_norm": float(
            np.linalg.norm(G_leak - np.diag(np.diag(G_leak)))
        ),
        "leakage_hermiticity_defect": leak_herm,
        "pullback_hermiticity_defect": pull_herm,
        "leakage_min_eigenvalue": float(leak_ev[0]) if len(leak_ev) else 0.0,
        "pullback_spectral_norm": scale,
        "input_amplitude_l1_sum": float(input_l1),
        "input_amplitude_l2_sum_squares": float(input_l2_sq),
        "pruning_performed_here": False,
        "exactness_condition": (
            "exact relative to supplied sparse columns only if each column is the "
            "complete unpruned coherent H e_j expansion in the correct orthonormal "
            "basis (or the appropriate Gram-weighted basis is used)"
        ),
    }
    return LeakageGramResult(C, G_leak, pullback, keys, diagnostics)


def classify_psd_spectrum(
    matrix: np.ndarray,
    rel_tol: float = 1e-10,
    eps_factor: float = 64.0,
) -> dict:
    """Scale-invariant PSD/nullspace diagnostic; no unit absolute scale floor.

    The zero threshold is relative to the matrix's own spectral scale. An
    identically zero matrix is treated as entirely null. Significant negative
    eigenvalues are reported as a PSD failure, not silently classified as zero.
    """
    for name, value in (("rel_tol", rel_tol), ("eps_factor", eps_factor)):
        if (isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, float, np.integer, np.floating))
                or not np.isfinite(value) or value < 0.0):
            raise ValueError(f"{name} must be finite and nonnegative")
    A = np.asarray(matrix, dtype=np.complex128)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("matrix must be square")
    if not (np.all(np.isfinite(A.real)) and np.all(np.isfinite(A.imag))):
        raise ValueError("matrix contains non-finite entries")
    herm_def = float(np.linalg.norm(A - A.conj().T))
    # Halve before addition to avoid overflow for finite near-maximal entries.
    H = 0.5 * A + 0.5 * A.conj().T
    if not (np.all(np.isfinite(H.real)) and np.all(np.isfinite(H.imag))):
        raise FloatingPointError("Hermitian representative is not representable")
    ev, U = np.linalg.eigh(H)
    if not (np.all(np.isfinite(ev)) and np.all(np.isfinite(U.real))
            and np.all(np.isfinite(U.imag))):
        raise FloatingPointError("PSD eigensolver returned non-finite values")
    scale = float(np.max(np.abs(ev))) if len(ev) else 0.0
    roundoff_rel = float(eps_factor) * np.finfo(float).eps * max(1, A.shape[0])
    threshold = max(float(rel_tol), roundoff_rel) * scale
    if not np.isfinite(threshold):
        raise FloatingPointError("relative PSD threshold exceeds float64 range")
    if scale == 0.0:
        zero = np.ones(len(ev), dtype=bool)
        significant_negative = False
    else:
        zero = np.abs(ev) <= threshold
        significant_negative = bool(np.min(ev) < -threshold)
    zero_projector = U[:, zero] @ U[:, zero].conj().T
    if not (np.all(np.isfinite(zero_projector.real))
            and np.all(np.isfinite(zero_projector.imag))):
        raise FloatingPointError("zero projector contains non-finite entries")
    return {
        "eigenvalues": ev,
        "eigenvectors": U,
        "spectral_scale": scale,
        "relative_zero_threshold": threshold,
        "nullity": int(np.sum(zero)),
        "positive_rank": int(np.sum(ev > threshold)),
        "significant_negative_eigenvalue": significant_negative,
        "hermiticity_defect": herm_def,
        "psd_pass": bool(
            not significant_negative
            and herm_def <= 1e-10 * max(scale, np.finfo(float).tiny)
        ),
        "zero_projector": zero_projector,
    }


def _self_test() -> None:
    # Two shell columns reach the same outside key. The off-diagonal imaginary
    # entry is coherent interference and must survive.
    keys = ("a", "b")
    cols = (
        {"a": 1.0, "outside-x": 1.0},
        {"b": 1.0, "outside-x": 1.0j},
    )
    got = assemble_leakage_gram(keys, cols)
    expected_leak = np.array([[1.0, 1.0j], [-1.0j, 1.0]], dtype=complex)
    assert np.allclose(got.leakage_gram, expected_leak, atol=1e-14)
    assert abs(got.leakage_gram[0, 1]) > 0.9

    full = np.array([[1, 0], [0, 1], [1, 1j]], dtype=complex)
    assert np.allclose(got.pullback_gram, full.conj().T @ full, atol=1e-14)

    base = np.diag([0.0, 2.0, 5.0]).astype(complex)
    ref = classify_psd_spectrum(base)
    for factor in (1e-24, 1e-12, 1e12, 1e24):
        cur = classify_psd_spectrum(factor * base)
        assert cur["nullity"] == ref["nullity"] == 1
        assert cur["positive_rank"] == ref["positive_rank"] == 2
        assert cur["psd_pass"]
        assert np.allclose(cur["zero_projector"], ref["zero_projector"], atol=1e-12)

    bad = classify_psd_spectrum(np.diag([-1.0, 2.0]))
    assert bad["significant_negative_eigenvalue"] and not bad["psd_pass"]


if __name__ == "__main__":
    _self_test()
    print("PASS: coherent leakage Gram, full pullback identity, scale invariance, PSD rejection")
