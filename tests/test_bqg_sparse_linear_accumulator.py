import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_sparse_linear_accumulator import (
    apply_sparse_matrix_columns,
    sparsify_matrix_columns,
    sparsify_matrix_columns_with_error_bound,
)


class TestSparseLinearAccumulator(unittest.TestCase):
    def test_fixed_mask_reports_rigorous_absolute_spectral_error_bound(self):
        columns = [
            {"r0": 3.0 + 4.0j, "r1": 0.3 + 0.4j},
            {"r0": 0.0, "r1": 0.6},
        ]
        kept, report = sparsify_matrix_columns_with_error_bound(columns, 0.5)
        self.assertEqual(kept, [{"r0": 3.0 + 4.0j}, {"r1": 0.6}])
        expected_discarded = math.hypot(0.5, 0.0)
        self.assertAlmostEqual(report["discarded_frobenius_norm"], expected_discarded)
        self.assertEqual(report["spectral_norm_error_upper_bound"], expected_discarded)
        self.assertEqual(report["discarded_nonzero_count"], 1)
        self.assertEqual(report["zero_entry_count"], 1)
        self.assertEqual(report["retained_nonzero_count"], 2)
        self.assertEqual(report["bound"], "||A-A_tau||_2 <= ||A-A_tau||_F")
        self.assertEqual(report["scope"], (
            "post-assembly elementwise matrix mask only; excludes all "
            "pre-assembly pruning, spin-wall truncation, and support caps"
        ))

    def test_bound_controls_error_for_arbitrary_test_vector(self):
        columns = [
            {"a": 0.2 + 0.1j, "b": 2.0},
            {"a": 0.3, "c": 0.4j},
        ]
        kept, report = sparsify_matrix_columns_with_error_bound(columns, 0.5)
        vector = [1.5 - 0.5j, -2.0 + 0.25j]
        full_out = apply_sparse_matrix_columns(columns, vector)
        sparse_out = apply_sparse_matrix_columns(kept, vector)
        rows = set(full_out) | set(sparse_out)
        error_norm = math.sqrt(sum(abs(full_out.get(k, 0j) - sparse_out.get(k, 0j)) ** 2 for k in rows))
        vector_norm = math.sqrt(sum(abs(x) ** 2 for x in vector))
        self.assertLessEqual(
            error_norm,
            report["spectral_norm_error_upper_bound"] * vector_norm + 1e-14,
        )

    def test_legacy_api_matches_certified_helper(self):
        columns = [{"x": 1e-4, "y": 1.0}, {"z": 2e-4}]
        expected = sparsify_matrix_columns(columns, 1e-3)
        actual, report = sparsify_matrix_columns_with_error_bound(columns, 1e-3)
        self.assertEqual(actual, expected)
        self.assertEqual(report["discarded_nonzero_count"], 2)

    def test_empty_and_zero_matrices_have_zero_bound(self):
        kept, report = sparsify_matrix_columns_with_error_bound([{}, {"x": 0j}], 0.0)
        self.assertEqual(kept, [{}, {}])
        self.assertEqual(report["spectral_norm_error_upper_bound"], 0.0)
        self.assertEqual(report["relative_frobenius_error_bound"], 0.0)

    def test_unrepresentable_frobenius_bound_is_rejected(self):
        columns = [{"a": 1e308}, {"b": 1e308}]
        with self.assertRaisesRegex(ValueError, "no finite error certificate"):
            sparsify_matrix_columns_with_error_bound(columns, 0.0)

    def test_invalid_tolerance_and_nonfinite_entries_are_rejected(self):
        with self.assertRaises(ValueError):
            sparsify_matrix_columns_with_error_bound([{"x": 1.0}], -1.0)
        with self.assertRaises(ValueError):
            sparsify_matrix_columns_with_error_bound([{"x": complex(float("inf"), 0)}], 0.1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
