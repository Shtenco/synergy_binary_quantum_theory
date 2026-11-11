#!/usr/bin/env python3
"""Regression gates for the finite Peter-Weyl master-projector audit."""
from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np
import peter_weyl_master_projector_pilot as pilot


class MasterProjectorIntegrityTests(unittest.TestCase):
    def test_master_assembly_returns_raw_nonhermitian_matrix(self):
        B = np.zeros((1, 1, 2, 2), dtype=complex)
        B[0, 0] = np.array([[1.0, 2.0j], [0.0, 1.0]])
        G = np.ones((1, 1), dtype=complex)
        raw = pilot.master_from_gram(B, G)
        self.assertGreater(np.linalg.norm(raw - raw.conj().T), 0.0)
        absolute, relative = pilot.hermitian_diagnostics(raw)
        self.assertGreater(absolute, 0.0)
        self.assertGreater(relative, 0.0)

    def test_hermitian_diagnostics_detects_defect_before_symmetrization(self):
        M = np.array([[1.0, 1.0j], [0.0, 2.0]], dtype=complex)
        absolute, relative = pilot.hermitian_diagnostics(M)
        self.assertAlmostEqual(absolute, np.sqrt(2.0))
        self.assertGreater(relative, 0.0)

    def test_negative_roundoff_is_not_promoted_to_zero_mode(self):
        evals = np.array([-5e-10, 1.0])
        mask = pilot.classify_spectrum(
            evals, zero_threshold=1e-10, psd_tolerance=1e-9
        )
        np.testing.assert_array_equal(mask, [False, False])

    def test_nonpositive_spectrum_beyond_psd_tolerance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not PSD"):
            pilot.classify_spectrum(
                np.array([-5e-10, 1.0]),
                zero_threshold=1e-10,
                psd_tolerance=1e-11,
            )

    def test_true_small_positive_and_negative_roundoff_are_separated(self):
        evals = np.array([-5e-12, 5e-12, 1.0])
        mask = pilot.classify_spectrum(
            evals, zero_threshold=1e-10, psd_tolerance=1e-10
        )
        np.testing.assert_array_equal(mask, [True, True, False])

    def test_projector_uses_absolute_eigenvalue_threshold(self):
        evals = np.array([-5e-10, 5e-11, 1.0])
        U = np.eye(3, dtype=complex)
        P, mask = pilot.projector_from_eigh(evals, U, 1e-10)
        np.testing.assert_array_equal(mask, [False, True, False])
        np.testing.assert_allclose(P, np.diag([0.0, 1.0, 0.0]))

    def test_constraint_columns_explicitly_forwards_no_prune_policy(self):
        basis = [("basis-state",)]
        observed = []
        def fake_apply(state, vertex, wall, mode="legacy"):
            observed.append(("apply", mode, vertex, wall))
            return {"out": 1e-14}
        def fake_prune(state, tol, mode="legacy"):
            observed.append(("prune", mode, tol))
            return dict(state)
        with patch.object(pilot.PW, "apply_H_cached_state", side_effect=fake_apply), \
             patch.object(pilot.PW, "prune_state", side_effect=fake_prune):
            columns, supports = pilot.constraint_columns(
                basis, (0,), 5, mode="no_prune"
            )
        self.assertEqual(columns[(0, 0)], {"out": 1e-14})
        self.assertEqual(supports["v0_b0"], 1)
        self.assertEqual(observed, [
            ("apply", "no_prune", 0, 5),
            ("prune", "no_prune", 1e-10),
        ])

    def test_constraint_columns_rejects_unknown_mode(self):
        with self.assertRaisesRegex(ValueError, "unsupported pruning mode"):
            pilot.constraint_columns([], (), 5, mode="magic")

    def test_sparse_inner_product_uses_compensated_complex_sums(self):
        a = {"large_pos": 1e16 + 1e16j, "small": 1.0 + 2.0j,
             "large_neg": -1e16 - 1e16j}
        b = {key: 1.0 + 0.0j for key in a}
        self.assertEqual(pilot.inner_sparse(a, b), 1.0 - 2.0j)

    def test_cross_gram_master_matches_independent_dense_column_assembly(self):
        nodes = (0, 1)
        dim = 2
        columns = {
            (0, 0): {"x": 1.0 + 2.0j, "y": 0.5},
            (0, 1): {"x": -0.25j, "z": 2.0},
            (1, 0): {"y": 1.5j, "z": 0.25},
            (1, 1): {"x": 1.0, "z": -1.0j},
        }
        B = pilot.cross_gram(columns, nodes, dim)
        G = np.array([[2.0, 0.3 + 0.4j], [0.3 - 0.4j, 1.5]], dtype=complex)
        gram_master = pilot.master_from_gram(B, G)

        # Independent direct construction: align every sparse output key into
        # one dense row space, then form sum_vw G_vw H_v^dag H_w.
        output_keys = sorted(
            {key for mapping in columns.values() for key in mapping}
        )
        dense = {}
        for v in nodes:
            H = np.zeros((len(output_keys), dim), dtype=complex)
            for a in range(dim):
                for row, key in enumerate(output_keys):
                    H[row, a] = columns[(v, a)].get(key, 0j)
            dense[v] = H
        direct = np.zeros((dim, dim), dtype=complex)
        for v in nodes:
            for w in nodes:
                direct += G[v, w] * (dense[v].conj().T @ dense[w])
        np.testing.assert_allclose(gram_master, direct, rtol=1e-13, atol=1e-13)

    def test_spectral_classifier_rejects_nonfinite_tolerances(self):
        with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
            pilot.classify_spectrum(
                np.array([0.0, 1.0]), zero_threshold=float("nan"),
                psd_tolerance=1e-9
            )

    def test_hermitian_relative_defect_is_scale_invariant(self):
        base = np.array([[1.0, 1.0j], [0.0, 2.0]], dtype=complex)
        relative = [
            pilot.hermitian_diagnostics(base * scale)[1]
            for scale in (1e-150, 1e-12, 1.0, 1e12, 1e150)
        ]
        np.testing.assert_allclose(relative, np.full(5, relative[2]), rtol=2e-15)

    def test_hermitian_diagnostics_rejects_nonfinite_matrix(self):
        with self.assertRaisesRegex(ValueError, "non-finite"):
            pilot.hermitian_diagnostics(np.array([[complex(float("nan"), 0.0)]]))

    def test_annihilation_keeps_tiny_nonzero_coefficients(self):
        columns = {(0, 0): {"x": 1e15 + 0j}}
        norm = pilot.constraint_annihilation_norm(
            columns, np.array([1e-15 + 0j]), nodes=(0,)
        )
        self.assertAlmostEqual(norm, 1.0, places=14)

    def test_annihilation_skips_only_exactly_zero_coefficients(self):
        columns = {
            (0, 0): {"x": 1e15 + 0j},
            (0, 1): {"x": -1e15 + 0j},
        }
        # The tiny first contribution must survive and cancel only through
        # the explicitly supplied second contribution, not a magic cutoff.
        norm = pilot.constraint_annihilation_norm(
            columns, np.array([1e-15 + 0j, 0j]), nodes=(0,)
        )
        self.assertAlmostEqual(norm, 1.0, places=14)

    def test_spectral_thresholds_are_scale_invariant(self):
        base = np.array([0.0, 1e-12])
        masks = []
        normalized_thresholds = []
        for scale in (1e-12, 1e-6, 1.0, 1e6, 1e12):
            values = base * scale
            spectral_scale, zero, psd = pilot.spectral_thresholds(values, 1e-10, 1e-10)
            masks.append(pilot.classify_spectrum(values, zero, psd).tolist())
            normalized_thresholds.append((zero / spectral_scale, psd / spectral_scale))
        self.assertEqual(masks, [[True, False]] * 5)
        np.testing.assert_allclose(normalized_thresholds, [[1e-10, 1e-10]] * 5)

    def test_relative_annihilation_residual_is_scale_invariant(self):
        observed = []
        for scale in (1e-100, 1e-12, 1.0, 1e12, 1e100):
            columns = {(0, 0): {"x": scale + 0j}}
            result = pilot.constraint_annihilation_metrics(
                columns, np.array([1.0 + 0j]), nodes=(0,)
            )
            observed.append(result["max_relative_residual"])
        np.testing.assert_allclose(observed, np.ones(5), rtol=2e-15, atol=0.0)

    def test_exact_common_kernel_has_zero_relative_annihilation(self):
        columns = {
            (0, 0): {"x": 1.0 + 0j},
            (0, 1): {"x": 1.0 + 0j},
        }
        coeff = np.array([1.0, -1.0], dtype=complex) / np.sqrt(2.0)
        result = pilot.constraint_annihilation_metrics(columns, coeff, nodes=(0,))
        self.assertEqual(result["max_absolute_residual"], 0.0)
        self.assertEqual(result["max_relative_residual"], 0.0)

    def test_stable_complex_norm_avoids_naive_square_overflow(self):
        observed = pilot.stable_complex_norm([1.2e308 + 1.2e308j])
        expected = np.hypot(1.2e308, 1.2e308)
        self.assertAlmostEqual(observed / expected, 1.0, places=14)

    def test_volume_hermiticity_registry_fails_closed(self):
        self.assertFalse(
            pilot.volume_raw_hermiticity_diagnostics({}, rtol=1e-10)["passed"]
        )
        good = {(1, 1, 1): {
            "relative_antihermitian_defect": 1e-14, "raw_norm": 2.0,
            "symmetrization_delta_frobenius_norm": 1e-14,
            "symmetrization_delta_spectral_norm_upper_bound": 1e-14,
            "symmetrization_relative_frobenius_bound": 5e-15,
            "hermitian_spectrum_resolution": {
                "spectral_scale": 1.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 1e-13,
                "resolved_numerical_nullity": 0,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 1.0,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.5,
            },
        }}
        self.assertTrue(
            pilot.volume_raw_hermiticity_diagnostics(good, rtol=1e-10)["passed"]
        )
        bad = {(1, 1, 1): {
            "relative_antihermitian_defect": 1e-4, "raw_norm": 2.0,
            "symmetrization_delta_frobenius_norm": 1e-4,
            "symmetrization_delta_spectral_norm_upper_bound": 1e-4,
            "symmetrization_relative_frobenius_bound": 5e-5,
            "hermitian_spectrum_resolution": {
                "spectral_scale": 1.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 1e-13,
                "resolved_numerical_nullity": 0,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 1.0,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.5,
            },
        }}
        self.assertFalse(
            pilot.volume_raw_hermiticity_diagnostics(bad, rtol=1e-10)["passed"]
        )
        nonfinite = {(1, 1, 1): {
            "relative_antihermitian_defect": float("nan"), "raw_norm": 2.0,
            "symmetrization_delta_frobenius_norm": 0.0,
            "symmetrization_delta_spectral_norm_upper_bound": 0.0,
            "symmetrization_relative_frobenius_bound": 0.0,
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(nonfinite, rtol=1e-10)
        self.assertFalse(result["passed"])
        self.assertIsNone(result["records"][0]["relative_antihermitian_defect"])
        self.assertIsNone(result["max_relative_defect"])

    def test_volume_registry_detects_incomplete_unique_call_coverage(self):
        good = {(1, 1, 1): {
            "relative_antihermitian_defect": 0.0, "raw_norm": 1.0,
            "symmetrization_delta_frobenius_norm": 0.0,
            "symmetrization_delta_spectral_norm_upper_bound": 0.0,
            "symmetrization_relative_frobenius_bound": 0.0,
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(
            good, rtol=1e-10, expected_unique_calls=2
        )
        self.assertFalse(result["passed"])
        self.assertFalse(result["coverage_complete"])
        self.assertEqual(result["coverage_shortfall"], 1)

    def test_volume_registry_checks_spectral_nullity_against_spin_dimension(self):
        base = {
            "relative_antihermitian_defect": 0.0,
            "raw_norm": 1.0,
            "symmetrization_delta_frobenius_norm": 0.0,
            "symmetrization_delta_spectral_norm_upper_bound": 0.0,
            "symmetrization_relative_frobenius_bound": 0.0,
        }
        invalid = {(1, 1, 1): {
            **base,
            "hermitian_spectrum_resolution": {
                "spectral_scale": 1.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 1e-13,
                "resolved_numerical_nullity": 8,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 1.0,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.5,
            },
        }}
        self.assertFalse(
            pilot.volume_raw_hermiticity_diagnostics(invalid)["passed"]
        )

        zero = {(0, 0, 0): {
            **base,
            "raw_norm": 0.0,
            "hermitian_spectrum_resolution": {
                "spectral_scale": 0.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 0.0,
                "resolved_numerical_nullity": 1,
                "min_resolved_abs_eigenvalue": None,
                "normalized_gap_to_zero": None,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": None,
            },
        }}
        self.assertTrue(
            pilot.volume_raw_hermiticity_diagnostics(zero)["passed"]
        )

    def test_volume_registry_rejects_inconsistent_spectral_sensitivity(self):
        registry = {(1, 1, 1): {
            "relative_antihermitian_defect": 0.0,
            "raw_norm": 2.0,
            "symmetrization_delta_frobenius_norm": 0.0,
            "symmetrization_delta_spectral_norm_upper_bound": 0.0,
            "symmetrization_relative_frobenius_bound": 0.0,
            "hermitian_spectrum_resolution": {
                "spectral_scale": 4.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 4e-13,
                "resolved_numerical_nullity": 0,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 0.25,
                # Incorrect: exact scalar derivative should be 0.5.
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.25,
            },
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(registry)
        self.assertFalse(result["passed"])

    def test_volume_registry_rejects_inconsistent_symmetrization_bound(self):
        registry = {(1, 1, 1): {
            "relative_antihermitian_defect": 1e-4,
            "raw_norm": 2.0,
            "symmetrization_delta_frobenius_norm": 1e-4,
            "symmetrization_delta_spectral_norm_upper_bound": 1e-4,
            "symmetrization_relative_frobenius_bound": 1e-8,
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(registry)
        self.assertFalse(result["passed"])

    def test_volume_registry_accepts_exact_unique_call_coverage(self):
        good = {
            (1, 1, 1): {
                "relative_antihermitian_defect": 0.0, "raw_norm": 1.0,
                "symmetrization_delta_frobenius_norm": 0.0,
                "symmetrization_delta_spectral_norm_upper_bound": 0.0,
                "symmetrization_relative_frobenius_bound": 0.0,
                "hermitian_spectrum_resolution": {
                "spectral_scale": 1.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 1e-13,
                "resolved_numerical_nullity": 0,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 1.0,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.5,
            },
            },
            (1, 1, 2): {
                "relative_antihermitian_defect": 1e-15, "raw_norm": 2.0,
                "symmetrization_delta_frobenius_norm": 1e-15,
                "symmetrization_delta_spectral_norm_upper_bound": 1e-15,
                "symmetrization_relative_frobenius_bound": 5e-16,
                "hermitian_spectrum_resolution": {
                "spectral_scale": 1.0,
                "eigenvalue_resolution_rtol": 1e-13,
                "eigenvalue_resolution_threshold": 1e-13,
                "resolved_numerical_nullity": 0,
                "min_resolved_abs_eigenvalue": 1.0,
                "normalized_gap_to_zero": 1.0,
                "max_scalar_sqrt_abs_derivative_on_resolved_modes": 0.5,
            },
            },
        }
        result = pilot.volume_raw_hermiticity_diagnostics(
            good, rtol=1e-10, expected_unique_calls=2
        )
        self.assertTrue(result["passed"])
        self.assertTrue(result["coverage_complete"])
        self.assertEqual(result["max_symmetrization_delta_frobenius_norm"], 1e-15)
        self.assertEqual(result["max_symmetrization_relative_frobenius_bound"], 5e-16)
        self.assertEqual(result["max_hermitian_spectral_scale"], 1.0)
        self.assertEqual(result["max_resolved_scalar_sqrt_abs_sensitivity"], 0.5)

    def test_volume_registry_rejects_boolean_spin_indices_and_limits(self):
        malformed = {(True, 0, 0): {
            "relative_antihermitian_defect": 0.0,
            "raw_norm": 1.0,
        }}
        self.assertFalse(
            pilot.volume_raw_hermiticity_diagnostics(malformed)["passed"]
        )
        with self.assertRaisesRegex(ValueError, "max_spin2"):
            pilot.volume_raw_hermiticity_diagnostics({}, max_spin2=True)
        with self.assertRaisesRegex(ValueError, "expected_unique_calls"):
            pilot.volume_raw_hermiticity_diagnostics({}, expected_unique_calls=True)

    def test_volume_registry_rejects_malformed_spin_key(self):
        malformed = {"not-a-spin-triple": {
            "relative_antihermitian_defect": 0.0, "raw_norm": 1.0
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(malformed)
        self.assertFalse(result["passed"])
        self.assertIsNone(result["records"][0]["spin_triple"])

    def test_volume_registry_rejects_malformed_record_without_crashing(self):
        malformed = {
            (1, 1, 1): None,
            (1, 1, 2): {"relative_antihermitian_defect": "bad", "raw_norm": 1.0},
        }
        result = pilot.volume_raw_hermiticity_diagnostics(malformed)
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["failures"]), 2)
        import json
        json.dumps(result, allow_nan=False)

    def test_volume_registry_rejects_spin_triple_above_declared_cutoff(self):
        registry = {(1, 1, 6): {
            "relative_antihermitian_defect": 0.0, "raw_norm": 1.0,
            "symmetrization_delta_frobenius_norm": 0.0,
            "symmetrization_delta_spectral_norm_upper_bound": 0.0,
            "symmetrization_relative_frobenius_bound": 0.0,
        }}
        result = pilot.volume_raw_hermiticity_diagnostics(
            registry, max_spin2=5
        )
        self.assertFalse(result["passed"])
        self.assertIsNone(result["records"][0]["spin_triple"])

    def test_absolute_annihilation_wrapper_keeps_zero_vector_compatibility(self):
        columns = {(0, 0): {"x": 1.0 + 0j}}
        self.assertEqual(
            pilot.constraint_annihilation_norm(
                columns, np.array([0.0 + 0j]), nodes=(0,)
            ),
            0.0,
        )

    def test_run_rejects_negative_tolerances(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            pilot.run(zero_rtol=-1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
