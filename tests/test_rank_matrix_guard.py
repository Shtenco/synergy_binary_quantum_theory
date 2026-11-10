import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proof_tools"))
from rank_matrix_guard import (validate_full_column_matrix, validate_shard_coverage, validate_index_coverage, validate_index_value, parse_indexed_map_filename)


class RankMatrixGuardTests(unittest.TestCase):
    def test_accepts_exact_finite_column_dimension(self):
        matrix = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        checked = validate_full_column_matrix(matrix, 2, "test orbit")
        np.testing.assert_array_equal(checked, matrix)

    def test_rejects_extra_columns_even_when_row_rank_equals_expected(self):
        # This matrix has rank 2, but three columns: rank == m alone is not
        # enough to certify full column rank of a claimed two-column block.
        matrix = np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
        self.assertEqual(np.linalg.matrix_rank(matrix), 2)
        with self.assertRaisesRegex(ValueError, "has 3 columns"):
            validate_full_column_matrix(matrix, 2, "test orbit")

    def test_rejects_missing_columns(self):
        with self.assertRaisesRegex(ValueError, "has 1 columns"):
            validate_full_column_matrix(np.ones((4, 1)), 2, "test orbit")

    def test_rejects_non_matrix_input(self):
        with self.assertRaisesRegex(ValueError, "must be 2D"):
            validate_full_column_matrix(np.ones(3), 3, "test orbit")

    def test_rejects_nonfinite_values(self):
        with self.assertRaisesRegex(ValueError, "NaN or infinity"):
            validate_full_column_matrix(np.array([[1.0, np.nan]]), 2, "test orbit")

    def test_accepts_complete_shard_coverage(self):
        self.assertEqual(validate_shard_coverage([2, 0, 1], 3, "test"), (0, 1, 2))

    def test_rejects_duplicate_shard_ids(self):
        with self.assertRaisesRegex(ValueError, "no duplicates"):
            validate_shard_coverage([0, 1, 1], 3, "test")

    def test_rejects_missing_shard_id(self):
        with self.assertRaisesRegex(ValueError, "expected 2 shard IDs"):
            validate_shard_coverage([0], 2, "test")

    def test_rejects_noninteger_shard_ids(self):
        with self.assertRaisesRegex(ValueError, "must be integers"):
            validate_shard_coverage([0, "1"], 2, "test")

    def test_rejects_boolean_shard_ids(self):
        with self.assertRaisesRegex(ValueError, "must be integers"):
            validate_shard_coverage([False, 1], 2, "test")

    def test_accepts_complete_record_and_map_indices(self):
        self.assertIsNone(validate_index_coverage([0, 1, 2], [2, 0, 1], 3, "test"))

    def test_rejects_record_index_gap(self):
        with self.assertRaisesRegex(ValueError, "records"):
            validate_index_coverage([0, 2, 3], [0, 1, 2], 3, "test")

    def test_rejects_map_index_gap(self):
        with self.assertRaisesRegex(ValueError, "maps"):
            validate_index_coverage([0, 1, 2], [0, 2, 3], 3, "test")

    def test_rejects_boolean_index(self):
        with self.assertRaisesRegex(ValueError, "integers"):
            validate_index_coverage([False, 1], [0, 1], 2, "test")

    def test_rejects_boolean_or_nonpositive_expected_dimension(self):
        for value in (True, 0, -1, 2.0):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "invalid expected column count"):
                    validate_full_column_matrix(np.eye(2), value, "test orbit")

    def test_validate_index_value_accepts_only_nonnegative_json_integers(self):
        self.assertEqual(validate_index_value(0, "record"), 0)
        for value in (True, False, "1", 1.0, -1, None):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "nonnegative JSON integer"):
                    validate_index_value(value, "record")

    def test_parse_map_filename_requires_canonical_decimal_index(self):
        self.assertEqual(parse_indexed_map_filename("b0", "map"), 0)
        self.assertEqual(parse_indexed_map_filename("b123", "map"), 123)
        for value in ("b", "b01", "b+1", "b-1", "b1x", "b１２"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "not canonical"):
                    parse_indexed_map_filename(value, "map")
        with self.assertRaisesRegex(ValueError, "must start with"):
            parse_indexed_map_filename("x1", "map")


if __name__ == "__main__":
    unittest.main()
