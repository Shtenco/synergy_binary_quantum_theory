import numpy as np
import pytest

from scripts.q2_pl_generation3_chunked_export import (
    _sha256_memmap,
    stable_child_id,
)


def test_stable_child_ids_are_parent_major_and_permutation_minor():
    assert stable_child_id(0, 0) == 0
    assert stable_child_id(0, 23) == 23
    assert stable_child_id(1, 0) == 24
    assert stable_child_id(9215, 23) == 221183
    with pytest.raises(ValueError):
        stable_child_id(-1, 0)
    with pytest.raises(ValueError):
        stable_child_id(0, 24)


def test_memmap_hash_is_independent_of_read_chunk_size(tmp_path):
    path = tmp_path / "array.npy"
    array = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.uint32, shape=(1003, 4)
    )
    array[:] = np.arange(4012, dtype=np.uint32).reshape(1003, 4)
    array.flush()
    assert _sha256_memmap(array, 17) == _sha256_memmap(array, 256)


def test_chunk_size_must_be_positive(tmp_path):
    from scripts.q2_pl_generation3_chunked_export import export
    with pytest.raises(ValueError, match="chunk_size"):
        export(tmp_path, chunk_size=0)


def test_external_merge_detects_duplicates_across_chunk_boundaries(tmp_path):
    from scripts.q2_pl_generation3_chunked_export import _count_unique_edges_external

    edges = np.array(
        [[0, 1], [2, 3], [1, 4], [0, 1], [3, 5], [2, 3]],
        dtype=np.uint32,
    )
    unique, duplicates = _count_unique_edges_external(
        edges, chunk_rows=2, scratch_dir=tmp_path
    )
    assert unique == 4
    assert duplicates == 2


def test_streaming_union_find_counts_components_without_sparse_graph():
    from scripts.q2_pl_generation3_chunked_export import _streaming_component_count

    edges = np.array([[0, 1], [1, 2], [3, 4]], dtype=np.uint32)
    assert _streaming_component_count(5, edges, chunk_rows=2) == 2
