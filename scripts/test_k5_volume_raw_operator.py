#!/usr/bin/env python3
"""Small-sector tests of the production K5 raw volume triple product.

These tests call the actual spin-matrix and volume routines from the K5
producer. They inspect the unsymmetrized Q before the historical Hermitian
completion, then test the derived sqrt(|Q|) matrix. They do not test the full
K5 HDA column or establish a continuum/physical volume theorem.
"""
from __future__ import annotations

import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


class K5RawVolumeOperatorTests(unittest.TestCase):
    CASES = ((1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2))

    def test_raw_triple_product_is_hermitian_before_symmetrization(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                q = k5.volume123_raw_operator(*spins)
                scale = max(np.linalg.norm(q), np.finfo(float).tiny)
                defect = np.linalg.norm(q - q.conj().T) / scale
                self.assertLessEqual(defect, 5e-14)

    def test_volume_matrix_is_hermitian_positive_semidefinite(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                volume = k5.volume123_matrix(*spins)
                np.testing.assert_allclose(
                    volume, volume.conj().T, rtol=2e-13, atol=2e-13
                )
                eigenvalues = np.linalg.eigvalsh(volume)
                self.assertGreaterEqual(float(eigenvalues.min()), -2e-12)

    def test_raw_defect_ledger_is_populated_from_unsymmetrized_operator(self):
        for spins in self.CASES:
            with self.subTest(spins=spins):
                k5.volume123_matrix(*spins)
                record = k5.VOLUME_RAW_DEFECTS[spins]
                q = k5.volume123_raw_operator(*spins)
                expected = np.linalg.norm(q - q.conj().T) / max(
                    np.linalg.norm(q), np.finfo(float).tiny
                )
                self.assertAlmostEqual(
                    record["relative_antihermitian_defect"], expected, places=14
                )
                self.assertAlmostEqual(record["raw_norm"], np.linalg.norm(q), places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
