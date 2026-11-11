#!/usr/bin/env python3
"""Finite-sector adjoint checks for the production K5 link-hit matrices.

This checks the local Clebsch-Gordan transition matrices and their orientation
reversal. It does not establish that the complete primitive/T/H operators are
mutual adjoints; those require full state-space matrix-element comparisons.
"""
from __future__ import annotations

import unittest

import numpy as np

import k5_peter_weyl_safe_hda_column as k5


class K5HitAdjointTests(unittest.TestCase):
    def test_reversed_hit_is_adjoint_for_every_binary_channel(self):
        # Doubled spins 1 -> 2 correspond to j=1/2 -> j=1.
        s_in, s_out = 1, 2
        for i in range(2):
            for j in range(2):
                with self.subTest(i=i, j=j):
                    left, right, norm = k5.hit_mats(
                        s_in, s_out, 0, 1, i, j
                    )
                    rev_left, rev_right, rev_norm = k5.hit_mats(
                        s_out, s_in, 1, 0, j, i
                    )
                    np.testing.assert_allclose(
                        rev_left, left.conj().T, rtol=2e-14, atol=2e-14
                    )
                    np.testing.assert_allclose(
                        rev_right, right.conj().T, rtol=2e-14, atol=2e-14
                    )
                    self.assertAlmostEqual(norm, rev_norm, places=14)

    def test_reversing_a_hit_path_twice_restores_original(self):
        path = (0, 2, 4, 1, 3)
        for iout in range(2):
            for iin in range(2):
                reversed_once = k5.adjoint_path(path, iout, iin)
                reversed_twice = k5.adjoint_path(*reversed_once)
                self.assertEqual(reversed_twice, (path, iout, iin))

    def test_adjoint_sequence_reversal_is_an_involution(self):
        seq = (
            ("V", 0),
            ("P", (0, 1, 2, 0), 1, 0),
            ("P", (2, 3, 4), 0, 1),
        )
        self.assertEqual(k5.adjoint_sequence(k5.adjoint_sequence(seq)), seq)


if __name__ == "__main__":
    unittest.main(verbosity=2)
