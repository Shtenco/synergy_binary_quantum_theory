import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from peter_weyl_cutoff_sixth_audit import exact_quaternion_haar_control


class ExactQuaternionHaarControlTests(unittest.TestCase):
    def test_exact_control_reproduces_existing_five_values(self):
        expected = [6, 120, 4650, 254604, 16807392]
        self.assertEqual(
            [exact_quaternion_haar_control(n) for n in range(1, 6)],
            expected,
        )

    def test_exact_sixth_moment(self):
        self.assertEqual(exact_quaternion_haar_control(6), 1250089632)


if __name__ == "__main__":
    unittest.main()
