#!/usr/bin/env python3
"""Synthetic dense-oracle tests for legacy/no-prune constraint stability."""
from __future__ import annotations

import json
import unittest
from unittest.mock import patch
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

    def test_subnormal_gram_rank_survives_threshold_underflow(self):
        gram = np.diag([0.0, 1.0]).astype(complex) * 1e-320
        result = stability.spectrum_of_gram(gram, 1e-10)
        self.assertEqual(result["rank_threshold_on_gram_eigenvalues"], 0.0)
        self.assertEqual(result["normalized_gram_eigenvalues"], [0.0, 1.0])
        self.assertEqual(result["estimated_rank"], 1)
        self.assertTrue(result["rank_threshold_resolves_roundoff"])

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


    def test_kernel_projector_distance_matches_analytic_principal_angle(self):
        theta = 0.173
        legacy = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        no_prune = stability.direct_svd_from_columns(
            {0: {(0, "r"): np.sin(theta)},
             1: {(0, "r"): -np.cos(theta)}}, 2
        )
        result = stability.compare_kernel_projectors(legacy, no_prune, 2)
        self.assertEqual(result["status"], "SAME_DIMENSION_SUBSPACE_DISTANCE")
        self.assertTrue(result["same_nullity"])
        self.assertAlmostEqual(
            result["max_principal_angle_sine"], abs(np.sin(theta)), places=12
        )
        self.assertLess(result["legacy_projector_idempotence_defect_frobenius"], 1e-12)
        self.assertLess(result["no_prune_projector_idempotence_defect_frobenius"], 1e-12)

    def test_kernel_projector_identical_subspaces_have_zero_distance(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        result = stability.compare_kernel_projectors(direct, direct, 2)
        self.assertEqual(result["status"], "SAME_DIMENSION_SUBSPACE_DISTANCE")
        self.assertAlmostEqual(result["projector_difference_spectral_norm"], 0.0, places=14)

    def test_kernel_nullity_mismatch_is_not_reported_as_principal_angle(self):
        legacy = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        no_prune = stability.direct_svd_from_columns(
            {0: {(0, "r0"): 1.0, (1, "r1"): 0.0},
             1: {(0, "r0"): 0.0, (1, "r1"): 1.0}}, 2
        )
        result = stability.compare_kernel_projectors(legacy, no_prune, 2)
        self.assertEqual(result["status"], "NULLITY_MISMATCH")
        self.assertFalse(result["same_nullity"])
        self.assertIsNone(result["max_principal_angle_sine"])
        self.assertEqual(result["legacy_nullity"], 1)
        self.assertEqual(result["no_prune_nullity"], 0)

    def test_kernel_comparison_blocks_scale_separation_warning(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        warned = dict(direct)
        warned["scale_separation_loss"] = True
        result = stability.compare_kernel_projectors(direct, warned, 2)
        self.assertEqual(result["status"], "INCONCLUSIVE_SCALE_SEPARATION")
        self.assertIsNone(result["max_principal_angle_sine"])


    def test_kernel_projector_rejects_nonunitary_full_right_singular_vectors(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        malformed = dict(direct)
        malformed["right_singular_vectors_real"] = np.zeros((2, 2)).tolist()
        with self.assertRaisesRegex(ValueError, "not orthonormal"):
            stability.kernel_projector_from_direct_svd(malformed, 2)

    def test_kernel_projector_rejects_rank_inconsistent_with_spectrum(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        malformed = dict(direct)
        malformed["rank"] = 0
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            stability.kernel_projector_from_direct_svd(malformed, 2)

    def test_kernel_projector_rejects_boolean_rank_and_dimension(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        malformed = dict(direct)
        malformed["rank"] = True
        with self.assertRaisesRegex(ValueError, "rank"):
            stability.kernel_projector_from_direct_svd(malformed, 2)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            stability.kernel_projector_from_direct_svd(direct, True)

    def test_cutoff_sweep_compares_adjacent_kernel_projectors(self):
        left = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        theta = 0.12
        right = stability.direct_svd_from_columns(
            {0: {(0, "r"): np.sin(theta)},
             1: {(0, "r"): -np.cos(theta)}}, 2
        )
        def result(direct, digest):
            maps = {
                0: {(0, "a"): 1.0},
                1: {(0, "b"): 2.0},
            } if digest == "cut1" else {
                0: {(0, "a"): 1.5},
                1: {(0, "b"): 2.0, (1, "c"): 0.25},
            }
            return {
                "domain_dimension": 2,
                "legacy": {
                    "direct_svd": direct,
                    "stacked_column_sha256": digest + "-L",
                },
                "no_prune": {
                    "direct_svd": direct,
                    "stacked_column_sha256": digest + "-N",
                },
                "_operator_maps": {"legacy": maps, "no_prune": maps},
            }
        with patch.object(
            stability, "run",
            side_effect=[result(left, "cut1"), result(right, "cut2")]
        ) as mocked:
            sweep = stability.cutoff_convergence_sweep(
                "low", [1, 3], block_rows=16
            )
        self.assertEqual(mocked.call_count, 2)
        self.assertTrue(all(call.kwargs["retain_operator_maps"] for call in mocked.call_args_list))
        self.assertEqual(sweep["doubled_spin_cutoffs"], [1, 3])
        self.assertEqual(len(sweep["points"]), 2)
        self.assertEqual(len(sweep["adjacent_cutoff_comparisons"]), 1)
        comparison = sweep["adjacent_cutoff_comparisons"][0]
        self.assertAlmostEqual(
            comparison["no_prune"]["operator_difference_norm_2"], 0.5, places=12
        )
        self.assertEqual(
            comparison["no_prune_kernel_subspace"]["status"],
            "SAME_DIMENSION_SUBSPACE_DISTANCE"
        )
        self.assertAlmostEqual(
            comparison["no_prune_kernel_subspace"]["max_principal_angle_sine"],
            abs(np.sin(theta)), places=12
        )

    def test_adjacent_cutoff_operator_difference_uses_common_row_union(self):
        previous = {
            0: {(0, "a"): 1.0},
            1: {(0, "b"): 2.0},
        }
        current = {
            0: {(0, "a"): 1.5},
            1: {(0, "b"): 2.0, (1, "new_state"): 0.25},
        }
        result = stability.adjacent_cutoff_operator_difference(
            previous, current, 2, singular_rtol=1e-12, block_rows=2
        )
        self.assertEqual(
            result["classification"], "FINITE_FLOAT64_COMMON_ROW_EMBEDDING"
        )
        self.assertAlmostEqual(result["operator_difference_norm_2"], 0.5, places=12)
        self.assertEqual(result["row_count_union"], 3)
        self.assertFalse(result["scale_separation_warning"])
        self.assertFalse(result["exact_operator_norm_certified"])

    def test_adjacent_cutoff_operator_difference_blocks_scale_loss(self):
        previous = {0: {(0, "a"): 1.0}, 1: {(0, "b"): 1.0}}
        current = {0: {(0, "a"): 1.0 + 1e-200}, 1: {(0, "b"): 1.0}}
        with patch.object(
            stability, "direct_svd_from_columns",
            return_value={"scale_separation_loss": True, "rank": 2,
                          "singular_values": [0.0, 0.0], "row_count": 2,
                          "block_count": 1}
        ):
            result = stability.adjacent_cutoff_operator_difference(
                previous, current, 2
            )
        self.assertEqual(result["classification"], "INCONCLUSIVE_SCALE_SEPARATION")
        self.assertIsNone(result["operator_difference_norm_2"])
        self.assertIsNone(result["difference_rank_at_declared_tolerance"])

    def test_cutoff_sweep_rejects_nonincreasing_or_invalid_cutoffs(self):
        with patch.object(stability, "run") as mocked:
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                stability.cutoff_convergence_sweep("low", [1, 1])
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                stability.cutoff_convergence_sweep("low", [3, 1])
            with self.assertRaisesRegex(ValueError, "integer >= 1"):
                stability.cutoff_convergence_sweep("low", [1, True])
            mocked.assert_not_called()

    def test_run_rejects_invalid_regulator_and_tolerance_arguments_before_work(self):
        with self.assertRaisesRegex(ValueError, "jmax2"):
            stability.run(jmax2=0)
        with self.assertRaisesRegex(ValueError, "jmax2"):
            stability.run(jmax2=True)
        with self.assertRaisesRegex(ValueError, "jmax2"):
            stability.run(jmax2=1.0)
        with self.assertRaisesRegex(ValueError, "rank_rtol"):
            stability.run(jmax2=1, rank_rtol=True)
        with self.assertRaisesRegex(ValueError, "singular_rtol"):
            stability.run(jmax2=1, singular_rtol=True)
        with self.assertRaisesRegex(ValueError, "block_rows"):
            stability.run(jmax2=1, block_rows=0)

    def test_scale_separation_blocks_kernel_and_threshold_rank_claims(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        warned = dict(direct)
        warned["scale_separation_loss"] = True
        warned["singular_values"] = [0.0, 0.0]
        kernel = stability.compare_kernel_projectors(direct, warned, 2)
        self.assertEqual(kernel["status"], "INCONCLUSIVE_SCALE_SEPARATION")
        self.assertIsNone(kernel["no_prune_rank"])
        standalone = stability.singular_threshold_sweep(
            warned, 2, [1e-12, 1e-8]
        )
        self.assertEqual(standalone["status"], "INCONCLUSIVE_SCALE_SEPARATION")
        self.assertEqual(standalone["threshold_results"], [])
        sweep = stability.compare_threshold_sweeps(
            direct, warned, 2, [1e-12, 1e-8]
        )
        self.assertEqual(sweep["status"], "INCONCLUSIVE_SCALE_SEPARATION")
        self.assertEqual(sweep["comparisons"], [])

    def test_singular_threshold_sweep_detects_rank_plateaus_and_transitions(self):
        columns = {
            0: {(0, "r0"): 1.0},
            1: {(1, "r1"): 1e-6},
            2: {(2, "r2"): 1e-12},
        }
        direct = stability.direct_svd_from_columns(columns, 3)
        sweep = stability.singular_threshold_sweep(
            direct, 3, [1e-13, 1e-9, 1e-5]
        )
        self.assertEqual(
            [row["rank"] for row in sweep["threshold_results"]],
            [3, 2, 1],
        )
        self.assertEqual(
            [row["nullity"] for row in sweep["threshold_results"]],
            [0, 1, 2],
        )
        self.assertEqual(
            [row["rank"] for row in sweep["rank_plateaus_on_declared_grid"]],
            [3, 2, 1],
        )

    def test_singular_threshold_sweep_rejects_unsorted_or_invalid_grid(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 1.0}, 1: {(1, "s"): 1e-5}}, 2
        )
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            stability.singular_threshold_sweep(direct, 2, [1e-5, 1e-8])
        with self.assertRaisesRegex(ValueError, "finite"):
            stability.singular_threshold_sweep(direct, 2, [float("nan")])
        with self.assertRaisesRegex(ValueError, "boolean"):
            stability.singular_threshold_sweep(direct, 2, [True])
        with self.assertRaisesRegex(ValueError, "positive integer"):
            stability.singular_threshold_sweep(direct, True, [1e-8])

    def test_threshold_sweep_rejects_nonorthonormal_right_singular_vectors(self):
        direct = stability.direct_svd_from_columns(
            {0: {(0, "r"): 1.0}, 1: {(1, "s"): 1e-5}}, 2
        )
        direct["right_singular_vectors_real"] = np.zeros((2, 2)).tolist()
        with self.assertRaisesRegex(ValueError, "not orthonormal"):
            stability.singular_threshold_sweep(direct, 2, [1e-8])

    def test_threshold_sweep_comparison_reports_projector_distance(self):
        left = stability.direct_svd_from_columns(
            {0: {(0, "r"): 0.0}, 1: {(0, "r"): 1.0}}, 2
        )
        right = stability.direct_svd_from_columns(
            {0: {(0, "r"): np.sin(0.1)},
             1: {(0, "r"): -np.cos(0.1)}}, 2
        )
        result = stability.compare_threshold_sweeps(
            left, right, 2, [1e-12, 1e-8]
        )
        self.assertEqual(len(result["comparisons"]), 2)
        self.assertTrue(all(row["same_nullity"] for row in result["comparisons"]))
        self.assertAlmostEqual(
            result["comparisons"][0]["max_principal_angle_sine"],
            abs(np.sin(0.1)), places=12
        )

if __name__ == "__main__":
    unittest.main(verbosity=2)
