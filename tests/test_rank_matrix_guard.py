import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proof_tools"))
from rank_matrix_guard import validate_full_column_matrix


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

    def test_rejects_boolean_or_nonpositive_expected_dimension(self):
        for value in (True, 0, -1, 2.0):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "invalid expected column count"):
                    validate_full_column_matrix(np.eye(2), value, "test orbit")


if __name__ == "__main__":
    unittest.main()
