#!/usr/bin/env python3
"""Threshold-free sparse matrix assembly contract for finite W5 operators.

The input basis is fixed before matrix construction. Every returned coefficient
is retained (including tiny nonzero values); only exact complex zero is omitted
from sparse storage. This module deliberately does not discover a basis by
amplitude thresholding and does not claim to construct the BQG W5 basis.
"""
from __future__ import annotations
from math import sqrt
from typing import Callable, Hashable, Iterable, Mapping

Key = Hashable
State = Mapping[Key, complex]
ColumnOperator = Callable[[Key], State]


def assemble_columns(
    basis: Iterable[Key], column_operator: ColumnOperator
) -> tuple[tuple[Key, ...], tuple[dict[int, complex], ...], tuple[tuple[int, Key, complex], ...]]:
    """Assemble a fixed-basis matrix; return columns and exact support leakage.

    No tolerance is accepted by design. A nonzero output outside the supplied
    basis is reported as leakage and never silently discarded.
    """
    keys = tuple(basis)
    if len(set(keys)) != len(keys):
        raise ValueError("basis contains duplicate keys")
    index = {key: i for i, key in enumerate(keys)}
    columns: list[dict[int, complex]] = []
    leakage: list[tuple[int, Key, complex]] = []
    for col, key in enumerate(keys):
        raw = column_operator(key)
        out: dict[int, complex] = {}
        for out_key, value in raw.items():
            z = complex(value)
            if z == 0j:
                continue
            row = index.get(out_key)
            if row is None:
                leakage.append((col, out_key, z))
            else:
                out[row] = out.get(row, 0j) + z
        columns.append({r: z for r, z in out.items() if z != 0j})
    return keys, tuple(columns), tuple(leakage)


def apply_columns(columns: tuple[dict[int, complex], ...], x: Mapping[int, complex]) -> dict[int, complex]:
    """Apply a column-oriented sparse matrix without magnitude pruning."""
    out: dict[int, complex] = {}
    for col, xval in x.items():
        if xval == 0j:
            continue
        if col < 0 or col >= len(columns):
            raise IndexError(f"input basis index {col} is out of range")
        for row, a in columns[col].items():
            out[row] = out.get(row, 0j) + a * xval
    return {r: z for r, z in out.items() if z != 0j}


def adjoint_columns(columns: tuple[dict[int, complex], ...], dimension: int) -> tuple[dict[int, complex], ...]:
    """Return columns of the exact algebraic conjugate transpose."""
    rows: list[dict[int, complex]] = [dict() for _ in range(dimension)]
    for col, column in enumerate(columns):
        for row, value in column.items():
            if row < 0 or row >= dimension:
                raise IndexError(f"output basis index {row} is out of range")
            rows[row][col] = complex(value).conjugate()
    return tuple(rows)


def sparse_difference_norm(a: Mapping[int, complex], b: Mapping[int, complex]) -> float:
    return sqrt(sum(abs(a.get(k, 0j) - b.get(k, 0j)) ** 2 for k in set(a) | set(b)))


def sparse_norm(a: Mapping[int, complex]) -> float:
    return sqrt(sum(abs(z) ** 2 for z in a.values()))


def linearity_residuals(columns: tuple[dict[int, complex], ...], x: Mapping[int, complex], y: Mapping[int, complex], lam: complex) -> dict[str, float]:
    ax = apply_columns(columns, x)
    ay = apply_columns(columns, y)
    alx = apply_columns(columns, {k: lam * v for k, v in x.items()})
    axy = apply_columns(columns, {k: x.get(k, 0j) + y.get(k, 0j) for k in set(x) | set(y)})
    lam_ax = {k: lam * v for k, v in ax.items()}
    ax_plus_ay = {k: ax.get(k, 0j) + ay.get(k, 0j) for k in set(ax) | set(ay)}
    return {
        "homogeneity_absolute_residual": sparse_difference_norm(alx, lam_ax),
        "additivity_absolute_residual": sparse_difference_norm(axy, ax_plus_ay),
    }


def adjoint_residual(direct_adjoint: tuple[dict[int, complex], ...], forward: tuple[dict[int, complex], ...], dimension: int) -> dict[str, float]:
    expected = adjoint_columns(forward, dimension)
    err2 = 0.0
    scale2 = 0.0
    for actual_col, expected_col in zip(direct_adjoint, expected):
        keys = set(actual_col) | set(expected_col)
        err2 += sum(abs(actual_col.get(k, 0j) - expected_col.get(k, 0j)) ** 2 for k in keys)
        scale2 += sum(abs(expected_col.get(k, 0j)) ** 2 for k in keys)
    absolute = sqrt(err2)
    scale = sqrt(scale2)
    return {"absolute_frobenius_residual": absolute, "reference_frobenius_norm": scale,
            "relative_frobenius_residual": absolute / max(scale, 1e-300)}
