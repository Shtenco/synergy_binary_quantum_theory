#!/usr/bin/env python3
"""Synthetic dense-oracle tests for legacy/no-prune constraint stability."""
from __future__ import annotations

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

    def test_subtraction_preserves_tiny_nonzero_difference(self):
        left = {0: {(0, "x"): 1.0 + 0j}}
        right = {0: {(0, "x"): 1.0 + 1e-15j}}
        delta = stability.subtract_stacked(left, right, 1)
        self.assertEqual(delta[0][(0, "x")], 1e-15j)

    def test_spectrum_rejects_materially_non_psd_gram(self):
        with self.assertRaisesRegex(ValueError, "not PSD"):
            stability.spectrum_of_gram(np.diag([-1.0, 2.0]), 1e-10)

    def test_column_digest_is_order_independent(self):
        a = {0: {(0, "x"): 1.0 + 2.0j, (1, "y"): 3.0}}
        b = {0: {(1, "y"): 3.0, (0, "x"): 1.0 + 2.0j}}
        self.assertEqual(stability.column_digest(a), stability.column_digest(b))


if __name__ == "__main__":
    unittest.main(verbosity=2)
