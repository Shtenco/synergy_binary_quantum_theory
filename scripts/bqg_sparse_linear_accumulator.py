#!/usr/bin/env python3
"""Strict sparse complex accumulation for finite-dimensional operator prototypes.

This module is deliberately NOT wired into the V09/V10/V11 production path yet.
It provides a dependency-free reference implementation for L1/L3 tests:
accumulate all contributions first, then optionally freeze one matrix-level
sparsification mask. Never threshold a state amplitude during accumulation.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Hashable


def _finite_complex(value: complex, *, name: str) -> complex:
    z = complex(value)
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError(f"{name} must be finite")
    return z


class SparseLinearAccumulator:
    """Collect sparse complex contributions without state-amplitude pruning."""

    def __init__(self) -> None:
        self._terms: dict[Hashable, list[complex]] = {}

    def add(self, state: Mapping[Hashable, complex], scale: complex = 1.0) -> None:
        factor = _finite_complex(scale, name="scale")
        if factor == 0:
            return
        for key, value in state.items():
            z = _finite_complex(value, name="state amplitude") * factor
            if not (math.isfinite(z.real) and math.isfinite(z.imag)):
                raise ValueError("scaled state amplitude must be finite")
            if z != 0:
                self._terms.setdefault(key, []).append(z)

    def finish(self) -> dict[Hashable, complex]:
        out: dict[Hashable, complex] = {}
        for key, terms in self._terms.items():
            try:
                z = complex(
                    math.fsum(term.real for term in terms),
                    math.fsum(term.imag for term in terms),
                )
            except OverflowError as exc:
                raise ValueError("accumulated amplitude overflowed") from exc
            if not (math.isfinite(z.real) and math.isfinite(z.imag)):
                raise ValueError("accumulated amplitude must be finite")
            if z != 0:
                out[key] = z
        return out


def sparsify_matrix_columns(
    columns: Sequence[Mapping[Hashable, complex]], tolerance: float
) -> list[dict[Hashable, complex]]:
    """Apply one fixed elementwise mask to an already assembled matrix."""
    tol = float(tolerance)
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("tolerance must be finite and nonnegative")
    result: list[dict[Hashable, complex]] = []
    for column in columns:
        kept: dict[Hashable, complex] = {}
        for key, value in column.items():
            z = _finite_complex(value, name="matrix element")
            if abs(z) > tol:
                kept[key] = z
        result.append(kept)
    return result


def _stable_l2_norm(values: Sequence[float]) -> float:
    """Compute an L2 norm without squaring large finite values directly."""
    norm = 0.0
    for value in values:
        norm = math.hypot(norm, float(value))
    return norm


def sparsify_matrix_columns_with_error_bound(
    columns: Sequence[Mapping[Hashable, complex]], tolerance: float
) -> tuple[list[dict[Hashable, complex]], dict[str, object]]:
    """Sparsify an assembled matrix and certify an absolute spectral-error bound.

    For the discarded matrix E, the report uses ||E||_2 <= ||E||_F.
    This is a rigorous absolute operator-norm bound for the supplied finite
    matrix. It is not a relative spectral-norm bound and says nothing about
    truncation performed before this matrix was assembled.
    """
    tol = float(tolerance)
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("tolerance must be finite and nonnegative")

    result: list[dict[Hashable, complex]] = []
    all_abs: list[float] = []
    kept_abs: list[float] = []
    discarded_abs: list[float] = []
    input_entry_count = 0
    zero_entry_count = 0
    per_column: list[dict[str, float | int]] = []

    for column in columns:
        kept: dict[Hashable, complex] = {}
        column_all: list[float] = []
        column_kept: list[float] = []
        column_discarded: list[float] = []
        column_zero_entries = 0
        for key, value in column.items():
            z = _finite_complex(value, name="matrix element")
            magnitude = abs(z)
            input_entry_count += 1
            if magnitude == 0.0:
                zero_entry_count += 1
                column_zero_entries += 1
                continue
            column_all.append(magnitude)
            all_abs.append(magnitude)
            if magnitude > tol:
                kept[key] = z
                column_kept.append(magnitude)
                kept_abs.append(magnitude)
            else:
                column_discarded.append(magnitude)
                discarded_abs.append(magnitude)
        result.append(kept)
        per_column.append({
            "column_index": len(per_column),
            "input_entry_count": len(column),
            "input_nonzero_count": len(column_all),
            "zero_entry_count": column_zero_entries,
            "retained_nonzero_count": len(column_kept),
            "discarded_nonzero_count": len(column_discarded),
            "discarded_frobenius_norm": _stable_l2_norm(column_discarded),
        })

    matrix_fro = _stable_l2_norm(all_abs)
    retained_fro = _stable_l2_norm(kept_abs)
    discarded_fro = _stable_l2_norm(discarded_abs)
    report: dict[str, object] = {
        "tolerance": tol,
        "column_count": len(columns),
        "input_entry_count": input_entry_count,
        "input_nonzero_count": len(all_abs),
        "zero_entry_count": zero_entry_count,
        "retained_nonzero_count": len(kept_abs),
        "discarded_nonzero_count": len(discarded_abs),
        "matrix_frobenius_norm": matrix_fro,
        "retained_frobenius_norm": retained_fro,
        "discarded_frobenius_norm": discarded_fro,
        "spectral_norm_error_upper_bound": discarded_fro,
        "relative_frobenius_error_bound": (
            discarded_fro / matrix_fro if matrix_fro else 0.0
        ),
        "bound": "||A-A_tau||_2 <= ||A-A_tau||_F",
        "scope": (
            "post-assembly elementwise matrix mask only; excludes all "
            "pre-assembly pruning, spin-wall truncation, and support caps"
        ),
        "per_column": per_column,
    }
    return result, report


def apply_sparse_matrix_columns(
    columns: Sequence[Mapping[Hashable, complex]],
    coefficients: Sequence[complex],
) -> dict[Hashable, complex]:
    """Apply fixed sparse matrix columns linearly to a coefficient vector."""
    if len(columns) != len(coefficients):
        raise ValueError("column count must match coefficient count")
    accumulator = SparseLinearAccumulator()
    for column, coefficient in zip(columns, coefficients):
        accumulator.add(column, coefficient)
    return accumulator.finish()
