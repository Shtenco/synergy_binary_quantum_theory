#!/usr/bin/env python3
"""Synthetic dense-oracle tests for legacy/no-prune constraint stability."""
from __future__ import annotations

import json
import unittest
import numpy as np

import peter_weyl_constraint_stability as stability


class ConstraintStabilityTests(unittest.TestCase):
    def test_sparse_gram_matches_independent_dense_assembly(self):
        columns = {
            0: {"x": 1.0 + 2.0j, "y": 0.5},
            1: {"x": -0.25j, "z": 2.0},
        }
        G = stability.gram_from_columns(columns, 2)
        keys = sorted({k for col in columns.values() for k in col})
        dense = np.array([[columns[a].get(k, 0j) for a in range(2)] for k in keys])
        np.testing.assert_allclose(G, dense.conj().T @ dense, rtol=1e-13, atol=1e-13)
        self.assertLess(np.linalg.norm(G - G.conj().T), 1e-13)

    def test_delta_operator_norm_matches_dense_svd(self):
        legacy = {
            0: {(0, "x"): 1.0 + 0j, (1, "y"): 2.0 + 0j},
            1: {(0, "y"): 1.0 + 1.0j},
        }
        no_prune = {
            0: {(0, "x"): 1.1 + 0j, (1, "y"): 2.0 + 0j},
            1: {(0, "y"): 1.0 + 1.0j, (1, "z"): 0.2j},
        }
        delta = stability.subtract_stacked(legacy, no_prune, 2)
        gram = stability.gram_from_columns(delta, 2)
        observed = stability.spectrum_of_gram(gram, 1e-12)["largest_singular_value"]
        keys = sorted({k for col in delta.values() for k in col})
        dense = np.array([[delta[a].get(k, 0j) for a in range(2)] for k in keys])
        expected = float(np.linalg.svd(dense, compute_uv=False)[0])
        self.assertAlmostEqual(observed, expected, places=12)

    def test_streaming_row_blocks_match_dense_operator_rows(self):
        columns = {
            0: {(0, "x"): 1.0 + 2.0j, (0, "y"): 0.5},
            1: {(0, "x"): -0.25j, (1, "z"): 2.0},
        }
        blocks = list(stability.iter_operator_row_blocks(columns, 2, block_rows=2))
        observed = np.vstack(blocks)
        keys = sorted(
            {key for column in columns.values() for key in column}, key=repr
        )
        expected = np.array([
            [columns[a].get(key, 0j) for a in range(2)] for key in keys
        ])
        np.testing.assert_allclose(observed, expected, rtol=0.0, atol=0.0)

    def test_direct_svd_certificate_is_strict_json_serializable(self):
        columns = {
            0: {(0, "x"): 1.0 + 1.0j},
            1: {(0, "x"): 0.5 - 0.25j},
        }
        result = stability.direct_svd_from_columns(columns, 2, block_rows=1)
        encoded = json.dumps(result, allow_nan=False)
        self.assertTrue(encoded)
        self.assertEqual(
            len(result["right_singular_vectors_real"]),
            len(result["right_singular_vectors_imag"]),
        )
        self.assertEqual(result["exact_rank_certified"], False)

    def test_empty_sparse_support_is_explicit_zero_operator(self):
        blocks = list(stability.iter_operator_row_blocks({0: {}, 1: {}}, 2))
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0].shape, (1, 2))
        np.testing.assert_array_equal(blocks[0], np.zeros((1, 2), dtype=complex))
        result = stability.direct_svd_from_columns({0: {}, 1: {}}, 2)
        self.assertEqual(result["rank"], 0)
        np.testing.assert_array_equal(result["singular_values"], [0.0, 0.0])

    def test_subtraction_preserves_tiny_nonzero_difference(self):
        left = {0: {(0, "x"): 1.0 + 0j}}
        right = {0: {(0, "x"): 1.0 + 1e-15j}}
        delta = stability.subtract_stacked(left, right, 1)
        self.assertEqual(delta[0][(0, "x")], 1e-15j)

    def test_spectrum_rejects_materially_non_psd_gram(self):
        with self.assertRaisesRegex(ValueError, "not PSD"):
            stability.spectrum_of_gram(np.diag([-1.0, 2.0]), 1e-10)

    def test_gram_rank_classification_is_scale_invariant(self):
        base = np.diag([0.0, 1e-12]).astype(complex)
        ranks = [
            stability.spectrum_of_gram(base * scale, 1e-10)["estimated_rank"]
            for scale in (1e-12, 1e-6, 1.0, 1e6, 1e12)
        ]
        self.assertEqual(ranks, [1, 1, 1, 1, 1])

    def test_zero_gram_has_zero_scale_and_zero_rank(self):
        result = stability.spectrum_of_gram(np.zeros((3, 3), dtype=complex), 1e-10)
        self.assertEqual(result["gram_eigenvalue_scale"], 0.0)
        self.assertEqual(result["estimated_rank"], 0)
        self.assertEqual(result["rank_threshold_on_gram_eigenvalues"], 0.0)

    def test_subnormal_gram_rank_uses_normalized_threshold(self):
        gram = np.diag([0.0, 1e-320]).astype(complex)
        result = stability.spectrum_of_gram(gram, 1e-10)
        self.assertEqual(result["gram_eigenvalue_scale"], 1e-320)
        self.assertEqual(result["rank_threshold_on_gram_eigenvalues"], 0.0)
        self.assertEqual(result["normalized_rank_threshold_on_gram_eigenvalues"], 1e-10)
        self.assertEqual(result["estimated_rank"], 1)

    def test_column_digest_is_order_independent(self):
        a = {0: {(0, "x"): 1.0 + 2.0j, (1, "y"): 3.0}}
        b = {0: {(1, "y"): 3.0, (0, "x"): 1.0 + 2.0j}}
        self.assertEqual(stability.column_digest(a), stability.column_digest(b))


if __name__ == "__main__":
    unittest.main(verbosity=2)
