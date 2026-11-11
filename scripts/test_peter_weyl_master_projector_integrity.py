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

    def test_run_rejects_negative_tolerances(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            pilot.run(zero_rtol=-1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
