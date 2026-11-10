"""Shape and finiteness guards for finite master-rank certificates."""
from __future__ import annotations

import numpy as np


def validate_full_column_matrix(matrix, expected_columns: int, context: str):
    """Fail closed unless a stacked master image has exactly the claimed columns."""
    if type(expected_columns) is not int or expected_columns < 1:
        raise ValueError(f"{context}: invalid expected column count {expected_columns!r}")
    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValueError(f"{context}: master image must be 2D, got shape {array.shape}")
    if array.shape[1] != expected_columns:
        raise ValueError(
            f"{context}: master image has {array.shape[1]} columns, "
            f"but the orbit declares {expected_columns}"
        )
    if not np.isfinite(array).all():
        raise ValueError(f"{context}: master image contains NaN or infinity")
    return array
