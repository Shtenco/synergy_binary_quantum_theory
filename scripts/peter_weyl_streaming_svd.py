#!/usr/bin/env python3
"""Streaming TSQR + SVD for tall, narrow complex constraint maps.

The input is an iterable of row blocks of the original operator C. No Gram
matrix is formed. Each block and merge is reduced with Householder QR; the
final right singular vectors are obtained from an SVD of the compact R.

This module certifies properties of the supplied floating-point row blocks
only. It does not estimate the upstream operator-construction error.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class StreamingSVDResult:
    singular_values: np.ndarray
    right_singular_vectors: np.ndarray
    row_count: int
    block_count: int
    input_scale: float
    max_abs_component: float
    rank: int
    rank_rtol: float
    scale_separation_loss: bool = False
    method: str = "binary-tree Householder TSQR followed by SVD(R)"
    exact_rank_certified: bool = False


def _checked_block(block: np.ndarray, n_columns: int) -> np.ndarray:
    arr = np.asarray(block)
    if arr.ndim != 2:
        raise ValueError("each row block must be a two-dimensional array")
    if arr.shape[1] != n_columns:
        raise ValueError(f"row block has {arr.shape[1]} columns; expected {n_columns}")
    if not np.issubdtype(arr.dtype, np.number):
        raise TypeError("row blocks must contain numeric values")
    arr = np.asarray(arr, dtype=np.complex128 if np.iscomplexobj(arr) else np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError("row block contains non-finite entries")
    return arr


def _qr_r(matrix: np.ndarray) -> np.ndarray:
    # NumPy's LAPACK reduced QR uses Householder reflectors.
    _, r = np.linalg.qr(matrix, mode="reduced")
    return r


def streaming_svd(row_blocks: Iterable[np.ndarray], n_columns: int,
                  rank_rtol: float = 1e-12) -> StreamingSVDResult:
    """Compute the SVD of a tall/narrow row stream without materializing C.

    rank_rtol is a relative singular-value threshold:
    sigma_i > rank_rtol * sigma_max. A dynamic common scale reduces overflow
    risk while retaining a common scale across all blocks.
    """
    if not isinstance(n_columns, int) or n_columns <= 0:
        raise ValueError("n_columns must be a positive integer")
    if not np.isfinite(rank_rtol) or rank_rtol < 0:
        raise ValueError("rank_rtol must be finite and nonnegative")

    # Binary-counter TSQR: at most one compact R per tree level.
    levels: list[np.ndarray | None] = []
    scale = 0.0
    row_count = 0
    block_count = 0
    max_abs = 0.0
    complex_input = False
    scale_separation_loss = False

    for raw in row_blocks:
        block = _checked_block(raw, n_columns)
        if block.shape[0] == 0:
            continue
        block_count += 1
        row_count += int(block.shape[0])
        complex_input = complex_input or np.iscomplexobj(block)
        # Use the largest real/imaginary component as scale. Unlike abs(z),
        # this cannot overflow merely because both complex components are large.
        local_scale = float(max(np.max(np.abs(block.real)),
                                np.max(np.abs(block.imag)) if np.iscomplexobj(block) else 0.0))
        max_abs = max(max_abs, local_scale)
        if local_scale == 0.0:
            continue

        if local_scale > scale:
            factor = 0.0 if scale == 0.0 else scale / local_scale
            for i, old in enumerate(levels):
                if old is not None:
                    scaled_old = old * factor
                    if factor == 0.0 or np.any((old != 0.0) & (scaled_old == 0.0)):
                        scale_separation_loss = True
                    levels[i] = scaled_old
            scale = local_scale

        normalized = block / scale
        if np.any((block != 0.0) & (normalized == 0.0)):
            scale_separation_loss = True
        r = _qr_r(normalized)
        level = 0
        while True:
            if level == len(levels):
                levels.append(r)
                break
            if levels[level] is None:
                levels[level] = r
                break
            r = _qr_r(np.vstack((levels[level], r)))
            levels[level] = None
            level += 1

    if row_count == 0:
        raise ValueError("row stream contains no rows")

    if scale == 0.0:
        singular = np.zeros(n_columns, dtype=np.float64)
        vectors = np.eye(n_columns, dtype=np.complex128 if complex_input else np.float64)
        rank = 0
    else:
        compact_r = None
        # Deterministic merge of remaining binary-tree nodes.
        for r in reversed(levels):
            if r is not None:
                compact_r = r if compact_r is None else _qr_r(np.vstack((compact_r, r)))
        assert compact_r is not None
        # Full Vh is required when rows < columns to expose the full right nullspace.
        _, s_scaled, vh = np.linalg.svd(compact_r, full_matrices=True)
        singular = np.zeros(n_columns, dtype=np.float64)
        singular[:len(s_scaled)] = s_scaled
        scaled_singular = singular.copy()
        singular *= scale
        if np.any((scaled_singular != 0.0) & (singular == 0.0)):
            # The normalized SVD retained a nonzero value, but conversion back
            # to the input's physical scale erased it. Any downstream absolute
            # spectral comparison (notably the Weyl delta bound) is inconclusive.
            scale_separation_loss = True
        vectors = vh.conj().T
        sigma_max = float(s_scaled[0]) if len(s_scaled) else 0.0
        rank = (int(np.count_nonzero(s_scaled > rank_rtol * sigma_max))
                if sigma_max > 0.0 else 0)

    if not np.all(np.isfinite(singular)) or not np.all(np.isfinite(vectors)):
        raise FloatingPointError("SVD result is non-finite; input scale exceeds representable range")

    return StreamingSVDResult(
        singular_values=singular,
        right_singular_vectors=vectors,
        row_count=row_count,
        block_count=block_count,
        input_scale=scale,
        max_abs_component=max_abs,
        rank=rank,
        rank_rtol=float(rank_rtol),
        scale_separation_loss=scale_separation_loss,
    )


def classify_singular_values(singular_values: np.ndarray, rank_rtol: float) -> dict:
    """Scale-invariant rank classification with explicit singular-value units."""
    s = np.asarray(singular_values, dtype=np.float64)
    if s.ndim != 1 or not np.all(np.isfinite(s)) or np.any(s < 0.0):
        raise ValueError("singular_values must be a finite nonnegative vector")
    if not np.isfinite(rank_rtol) or rank_rtol < 0.0:
        raise ValueError("rank_rtol must be finite and nonnegative")
    scale = float(np.max(s)) if s.size else 0.0
    normalized = np.zeros_like(s) if scale == 0.0 else s / scale
    rank = 0 if scale == 0.0 else int(np.count_nonzero(normalized > rank_rtol))
    return {
        "rank": rank,
        "rank_rtol_on_singular_values": float(rank_rtol),
        "normalized_singular_values": normalized.tolist(),
        "normalized_threshold": float(rank_rtol),
        "spectral_scale": scale,
        "classification": "NUMERICAL_ONLY_NO_UPSTREAM_ERROR_BOUND",
    }
