import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bqg_v10_w5_amp_tol_shadow_compare import (
    closure_matrix_available, compare_vectors, unpack_columns,
)


class TestV10ShadowVectorComparison(unittest.TestCase):
    def test_compares_complex_sparse_vectors_and_support(self):
        reference = {
            "input-A": {"out-1": 1 + 2j, "out-2": 0.5 - 0.25j},
        }
        candidate = {
            "input-A": {"out-1": 1 + 2j, "out-3": 0.5 - 0.25j},
        }
        result = compare_vectors(reference, candidate)
        self.assertEqual(result["support_symmetric_difference_total"], 2)
        self.assertAlmostEqual(result["max_coefficient_absolute_residual"], abs(0.5 - 0.25j))
        self.assertGreater(result["absolute_l2_residual"], 0.0)

    def test_full_shadow_matrix_is_comparable_but_not_certified(self):
        shadow = {
            "status": "INCOMPLETE_NO_PRUNING_COVERAGE",
            "pruning_mode": "shadow_no_prune",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
            "forward_B_columns_sparse": [],
        }
        self.assertTrue(closure_matrix_available(shadow, shadow=True))
        self.assertFalse(closure_matrix_available(shadow, shadow=False))

    def test_shadow_matrix_is_not_comparable_when_closure_is_capped(self):
        shadow = {
            "status": "INCOMPLETE_CLOSURE_LIMIT_REACHED",
            "pruning_mode": "shadow_no_prune",
            "algebraic_status": "PASS_ADJOINT_CONSISTENCY",
        }
        self.assertFalse(closure_matrix_available(shadow, shadow=True))

    def test_decodes_full_complex_columns(self):
        result = {
            "basis_keys_repr": ["input", "output"],
            "forward_B_columns_sparse": [
                [[1, 0.25, -0.5]],
                [[0, 1.0, 2.0]],
            ],
        }
        decoded = unpack_columns(result)
        self.assertEqual(decoded["input"]["output"], 0.25 - 0.5j)
        self.assertEqual(decoded["output"]["input"], 1.0 + 2.0j)


if __name__ == "__main__":
    unittest.main()
