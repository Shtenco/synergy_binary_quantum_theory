#!/usr/bin/env python3
"""Synthetic dense-oracle tests for the K5 graph-shell leakage master."""
from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

import bqg_k5_graph_sector_master_shell as shell


class K5GraphShellLeakageTests(unittest.TestCase):
    def test_apply_H_to_key_requests_no_prune_without_second_prune(self):
        key = shell.K5.basis_full_jhalf()[0]
        expected = {"out": 1.0 + 2.0j}
        with patch.object(
            shell.K5, "apply_H_cached_state", return_value=expected
        ) as apply, patch.object(
            shell.K5, "prune_state", side_effect=AssertionError("second prune")
        ):
            got = shell.apply_H_to_key(key, 2, 5, "no_prune")
        self.assertEqual(got, expected)
        apply.assert_called_once_with(
            {key: 1.0 + 0.0j}, 2, 5, mode="no_prune"
        )

    def test_compressed_constraints_retains_coherent_offdiagonal_leakage(self):
        basis = ("a", "b")
        columns = [
            {"a": 1.0, "outside": 1.0},
            {"b": 1.0, "outside": 1.0j},
        ]
        with patch.object(shell.K5, "VERT", (0,)), patch.object(
            shell, "apply_H_to_key", side_effect=columns
        ):
            constraints, leakage_grams, ledger = shell.compressed_constraints(
                basis, 5, "no_prune"
            )
        np.testing.assert_allclose(constraints[0], np.eye(2), atol=0, rtol=0)
        expected = np.array([[1.0, 1.0j], [-1.0j, 1.0]], dtype=complex)
        np.testing.assert_allclose(leakage_grams[0], expected, atol=1e-14)
        self.assertGreater(abs(leakage_grams[0][0, 1]), 0.9)
        self.assertGreater(ledger[0]["coherent_offdiagonal_leakage_norm"], 1.0)

    def test_full_master_equals_projected_plus_coherent_leakage(self):
        C = np.eye(2, dtype=complex)
        G = np.array([[1.0, 1.0j], [-1.0j, 1.0]], dtype=complex)
        got = shell.master_matrices([C], [G])
        expected_projected = C.conj().T @ C
        expected_full = expected_projected + G
        np.testing.assert_allclose(got["projected"], expected_projected, atol=0, rtol=0)
        np.testing.assert_allclose(got["leakage"], G, atol=0, rtol=0)
        np.testing.assert_allclose(got["full"], expected_full, atol=0, rtol=0)
        self.assertGreater(np.linalg.norm(got["full"] - got["projected"]), 0.0)

    def test_raw_hermiticity_is_measured_before_hermitian_cleanup(self):
        C = np.eye(2, dtype=complex)
        deliberately_nonhermitian = np.array(
            [[1.0, 1.0j], [0.0, 1.0]], dtype=complex
        )
        got = shell.master_matrices([C], [deliberately_nonhermitian])
        defect = got["diagnostics"]["full_master_raw_hermiticity"][
            "relative_frobenius_defect"
        ]
        self.assertGreater(defect, 0.1)
        np.testing.assert_allclose(
            got["full"], got["full"].conj().T, atol=0, rtol=0
        )

    def test_relative_zero_projector_is_invariant_under_positive_rescaling(self):
        base = np.diag([0.0, 2.0, 5.0]).astype(complex)
        reference = shell.projector_from_master(base, reltol=1e-10)
        self.assertEqual(reference["zero_rank"], 1)
        for factor in (1e-24, 1e-12, 1e12, 1e24):
            current = shell.projector_from_master(factor * base, reltol=1e-10)
            self.assertEqual(current["zero_rank"], 1)
            self.assertTrue(current["psd_pass"])
            np.testing.assert_allclose(
                current["zero_projector"],
                reference["zero_projector"],
                atol=1e-12,
                rtol=0,
            )

    def test_identically_zero_master_has_full_nullity_without_unit_floor(self):
        result = shell.projector_from_master(np.zeros((3, 3), dtype=complex))
        self.assertEqual(result["zero_rank"], 3)
        self.assertIsNone(result["positive_gap"])
        self.assertEqual(result["zero_threshold"], 0.0)

    def test_significant_negative_master_eigenvalue_fails_psd(self):
        result = shell.projector_from_master(
            np.diag([-1.0, 2.0]).astype(complex)
        )
        self.assertFalse(result["psd_pass"])
        self.assertTrue(result["significant_negative_eigenvalue"])

    def test_leakage_gram_fails_closed_when_true_result_overflows_float64(self):
        with self.assertRaisesRegex(FloatingPointError, "float64 range"):
            shell.assemble_leakage_gram(
                ("inside",),
                ({"outside": 1e200 + 0j},),
            )

    def test_psd_classifier_handles_near_maximal_finite_hermitian_matrix(self):
        result = shell.classify_psd_spectrum(
            np.array([[1e308 + 0j]], dtype=complex)
        )
        self.assertTrue(result["psd_pass"])
        self.assertEqual(result["positive_rank"], 1)

    def test_psd_classifier_rejects_invalid_tolerance_controls(self):
        for rel_tol, eps_factor in ((float("nan"), 64.0), (1e-10, -1.0), (True, 64.0)):
            with self.subTest(rel_tol=rel_tol, eps_factor=eps_factor):
                with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                    shell.classify_psd_spectrum(
                        np.eye(2, dtype=complex),
                        rel_tol=rel_tol,
                        eps_factor=eps_factor,
                    )

    def test_generated_shell_checks_all_used_local_basis_sectors(self):
        basis = shell.K5.basis_full_jhalf()
        result = shell.validate_shell_basis_orthonormality(basis)
        self.assertTrue(result["passed"])
        self.assertGreater(result["local_sectors_checked"], 0)
        self.assertLess(result["maximum_relative_gram_defect"], 2e-12)

    def test_depth_two_requires_complete_pullback_spin_wall(self):
        self.assertEqual(shell.required_spin_wall2(2), 7)
        with patch.object(shell, "generate_shell") as generate:
            with self.assertRaisesRegex(ValueError, "complete depth-2 pullback wall"):
                shell.run(depth=2, jmax2=5)
        generate.assert_not_called()


if __name__ == "__main__":
    unittest.main()
