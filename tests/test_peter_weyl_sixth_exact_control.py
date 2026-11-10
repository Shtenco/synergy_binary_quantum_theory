import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from peter_weyl_cutoff_sixth_audit import (
    exact_quaternion_haar_control,
    sparse_transitions_without_hidden_threshold,
)


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

    def test_exact_sixth_moment(self):
        self.assertEqual(exact_quaternion_haar_control(6), 1250089632)


if __name__ == "__main__":
    unittest.main()
