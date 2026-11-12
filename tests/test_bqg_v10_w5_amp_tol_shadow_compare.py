import unittest

from bqg_v10_w5_amp_tol_shadow_compare import compare_vectors, unpack_columns


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
