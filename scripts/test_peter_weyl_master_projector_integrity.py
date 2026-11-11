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

    def test_run_rejects_negative_tolerances(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            pilot.run(zero_rtol=-1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
