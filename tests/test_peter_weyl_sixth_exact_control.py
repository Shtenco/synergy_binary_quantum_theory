import sys
import unittest
from pathlib import Path

import numpy as np
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from peter_weyl_cutoff_sixth_audit import (
    exact_quaternion_haar_control,
    sparse_transitions_without_hidden_threshold,
)
import peter_weyl_cutoff_staircase as PW


class ExactQuaternionHaarControlTests(unittest.TestCase):
    def test_exact_control_reproduces_existing_five_values(self):
        expected = [6, 120, 4650, 254604, 16807392]
        self.assertEqual(
            [exact_quaternion_haar_control(n) for n in range(1, 6)],
            expected,
        )

    def test_transition_conversion_does_not_hide_small_nonzero_entries(self):
        matrix = np.array([[1.0, 0.0], [1e-16, 1.0]], dtype=complex)
        transitions, diag = sparse_transitions_without_hidden_threshold(matrix)
        self.assertEqual(len(transitions[0]), 2)
        self.assertEqual(transitions[0][1][1], 1e-16)
        self.assertEqual(diag["transition_entries_at_or_below_1e-14"], 1)

    def test_link_builder_preserves_tiny_nonzero_cg_coefficients(self):
        with patch.object(PW, "CG_exact", return_value=PW.sp.Rational(1, 10**16)):
            _, _, _, _, U, _ = PW.build_link(0.5)
        self.assertTrue(any(
            value != 0.0 and abs(value) < 1e-14
            for matrix in U.values()
            for value in matrix.ravel()
        ))

    def test_exact_sixth_moment(self):
        self.assertEqual(exact_quaternion_haar_control(6), 1250089632)


if __name__ == "__main__":
    unittest.main()
