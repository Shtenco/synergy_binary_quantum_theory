#!/usr/bin/env python3
"""Independent dense-oracle tests for streaming Peter-Weyl TSQR + SVD."""
from __future__ import annotations

import unittest
import numpy as np

import peter_weyl_streaming_svd as tsqr


class StreamingSVDTests(unittest.TestCase):
    def test_matches_dense_complex_svd(self):
        rng = np.random.default_rng(71)
        a = rng.normal(size=(113, 9)) + 1j * rng.normal(size=(113, 9))
        blocks = [a[i:i + 13] for i in range(0, len(a), 13)]
        got = tsqr.streaming_svd(blocks, 9)
        expected = np.linalg.svd(a, compute_uv=False)
        np.testing.assert_allclose(got.singular_values, expected, rtol=2e-13, atol=2e-13)
        self.assertEqual(got.rank, 9)
        self.assertEqual(got.row_count, 113)

    def test_preserves_controlled_small_singular_values(self):
        rng = np.random.default_rng(72)
        q1, _ = np.linalg.qr(
            rng.normal(size=(160, 8)) + 1j * rng.normal(size=(160, 8))
        )
        q2, _ = np.linalg.qr(
            rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        )
        target = np.array([1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14])
        a = q1 @ np.diag(target) @ q2.conj().T
        got = tsqr.streaming_svd(
            [a[i:i + 17] for i in range(0, 160, 17)], 8
        )
        np.testing.assert_allclose(got.singular_values, target, rtol=2e-3, atol=1e-16)

    def test_rank_deficient_and_complete_right_nullspace(self):
        a = np.array([[1., 2., 3., 4.], [2., 4., 6., 8.]])
        got = tsqr.streaming_svd([a[:1], a[1:]], 4, rank_rtol=1e-12)
        self.assertEqual(got.right_singular_vectors.shape, (4, 4))
        self.assertEqual(got.rank, 1)
        np.testing.assert_allclose(
            a @ got.right_singular_vectors[:, 1:], 0.0, atol=1e-12
        )

    def test_scale_invariant_classifier(self):
        base = np.array([0.0, 1e-12, 1.0])
        ranks = [
            tsqr.classify_singular_values(base * c, 1e-10)["rank"]
            for c in (1e-12, 1e-6, 1.0, 1e6, 1e12)
        ]
        self.assertEqual(ranks, [1, 1, 1, 1, 1])

    def test_row_permutation_and_block_size_invariance(self):
        rng = np.random.default_rng(73)
        a = rng.normal(size=(97, 7))
        p = rng.permutation(len(a))
        x = tsqr.streaming_svd(
            [a[i:i + 11] for i in range(0, len(a), 11)], 7
        )
        y = tsqr.streaming_svd(
            [a[p][i:i + 19] for i in range(0, len(a), 19)], 7
        )
        np.testing.assert_allclose(
            x.singular_values, y.singular_values, rtol=2e-13, atol=2e-13
        )

    def test_streaming_svd_is_scale_invariant_over_wide_range(self):
        rng = np.random.default_rng(74)
        a = rng.normal(size=(41, 6)) + 1j * rng.normal(size=(41, 6))
        base = tsqr.streaming_svd([a[:9], a[9:]], 6)
        for factor in (1e-150, 1e150):
            scaled = tsqr.streaming_svd(
                [factor * a[:9], factor * a[9:]], 6
            )
            np.testing.assert_allclose(
                scaled.singular_values / factor,
                base.singular_values,
                rtol=5e-13,
                atol=5e-13,
            )
            self.assertEqual(scaled.rank, base.rank)

    def test_reports_underflow_from_extreme_block_scale_separation(self):
        blocks = [
            np.array([[1e308, 0.0]]),
            np.array([[0.0, 1e-308]]),
        ]
        result = tsqr.streaming_svd(blocks, 2)
        self.assertTrue(result.scale_separation_loss)

    def test_warns_on_nonzero_sub_ulp_component_even_without_underflow(self):
        matrix = np.array([[1.0, 0.0], [0.0, 1e-20]])
        result = tsqr.streaming_svd([matrix], 2)
        self.assertFalse(np.any(matrix == 0.0) is False)
        self.assertTrue(result.scale_separation_loss)
        self.assertEqual(result.rank, 2)

    def test_complex_component_scaling_avoids_abs_overflow(self):
        value = 1.2e308 + 1.2e308j
        result = tsqr.streaming_svd([np.array([[value]])], 1)
        expected = np.hypot(1.2e308, 1.2e308)
        self.assertAlmostEqual(result.singular_values[0] / expected, 1.0, places=14)
        self.assertFalse(result.scale_separation_loss)

    def test_all_zero_matrix(self):
        got = tsqr.streaming_svd(
            [np.zeros((5, 3)), np.zeros((2, 3))], 3
        )
        np.testing.assert_array_equal(got.singular_values, np.zeros(3))
        self.assertEqual(got.rank, 0)

    def test_rejects_nonfinite_and_bad_shapes(self):
        with self.assertRaisesRegex(ValueError, "non-finite"):
            tsqr.streaming_svd([np.array([[np.nan, 0.]])], 2)
        with self.assertRaisesRegex(ValueError, "two-dimensional"):
            tsqr.streaming_svd([np.array([1., 2.])], 2)

    def test_rejects_empty_stream(self):
        with self.assertRaisesRegex(ValueError, "no rows"):
            tsqr.streaming_svd([], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
