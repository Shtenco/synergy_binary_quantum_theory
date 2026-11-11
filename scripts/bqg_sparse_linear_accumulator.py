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
