#!/usr/bin/env python3
"""Small-sector tests of the production K5 raw volume triple product.

These tests call the actual spin-matrix and volume routines from the K5
producer. They inspect the unsymmetrized Q before the historical Hermitian
completion, then test the derived sqrt(|Q|) matrix. They do not test the full
K5 HDA column or establish a continuum/physical volume theorem.
"""
from __future__ import annotations

import itertools
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


    def test_raw_volume_is_covariant_under_all_leg_permutations(self):
        spins = (1, 2, 2)
        dims = tuple(s + 1 for s in spins)
        q = k5.volume123_raw_operator(*spins)
        old_indices = list(itertools.product(*(range(d) for d in dims)))
        old_flat = {index: n for n, index in enumerate(old_indices)}

        for permutation in itertools.permutations(range(3)):
            with self.subTest(permutation=permutation):
                new_spins = tuple(spins[i] for i in permutation)
                new_dims = tuple(dims[i] for i in permutation)
                new_indices = list(
                    itertools.product(*(range(d) for d in new_dims))
                )
                new_flat = {index: n for n, index in enumerate(new_indices)}
                unitary = np.zeros((len(new_indices), len(old_indices)), dtype=complex)
                for old_index in old_indices:
                    new_index = tuple(old_index[i] for i in permutation)
                    unitary[new_flat[new_index], old_flat[old_index]] = 1.0

                inversions = sum(
                    permutation[i] > permutation[j]
                    for i in range(3) for j in range(i + 1, 3)
                )
                parity = -1.0 if inversions % 2 else 1.0
                q_permuted = k5.volume123_raw_operator(*new_spins)
                np.testing.assert_allclose(
                    q_permuted,
                    parity * unitary @ q @ unitary.conj().T,
                    rtol=3e-13,
                    atol=3e-13,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
