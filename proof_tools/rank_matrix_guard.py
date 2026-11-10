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



def validate_index_value(value, context: str) -> int:
    """Reject bools, numeric strings and negative values as artifact indices."""
    if type(value) is not int or value < 0:
        raise ValueError(f"{context}: index must be a nonnegative JSON integer, got {value!r}")
    return value


def validate_positive_integer(value, context: str) -> int:
    """Require a strictly positive JSON integer for dimensions/counts."""
    checked = validate_index_value(value, context)
    if checked < 1:
        raise ValueError(f"{context}: value must be a positive integer, got {value!r}")
    return checked


def parse_indexed_map_filename(stem: str, context: str) -> int:
    """Parse canonical b<index> map stems without coercive int() acceptance."""
    if not isinstance(stem, str) or not stem.startswith("b"):
        raise ValueError(f"{context}: map filename must start with 'b', got {stem!r}")
    digits = stem[1:]
    if (
        not digits
        or not digits.isascii()
        or not digits.isdigit()
        or (len(digits) > 1 and digits.startswith("0"))
    ):
        raise ValueError(f"{context}: map filename index is not canonical: {stem!r}")
    return int(digits)


def validate_shard_coverage(shard_ids, expected_shards: int, context: str):
    """Require one unique shard ID for every integer in [0, expected_shards)."""
    if type(expected_shards) is not int or expected_shards < 1:
        raise ValueError(f"{context}: invalid expected shard count {expected_shards!r}")
    ids = list(shard_ids)
    if len(ids) != expected_shards:
        raise ValueError(
            f"{context}: expected {expected_shards} shard IDs, got {len(ids)}"
        )
    if any(type(shard_id) is not int for shard_id in ids):
        raise ValueError(f"{context}: shard IDs must be integers, got {ids!r}")
    expected = list(range(expected_shards))
    if sorted(ids) != expected:
        raise ValueError(
            f"{context}: shard IDs must cover exactly 0..{expected_shards - 1} "
            f"with no duplicates; got {sorted(ids)!r}"
        )
    return tuple(sorted(ids))


def validate_index_coverage(record_ids, map_ids, expected_count: int, context: str):
    """Require complete zero-based record and map indices with no gaps."""
    validate_shard_coverage(record_ids, expected_count, f"{context} records")
    validate_shard_coverage(map_ids, expected_count, f"{context} maps")
