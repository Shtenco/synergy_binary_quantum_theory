import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from bqg_v09_w5_lorentzian_raw_column import (
    w5_input_key, covariant_scalar_to_gauss, ordered_epsilon_triples
)

class W5LorentzianAdapterTests(unittest.TestCase):
    def test_exact_w5_basis_keys(self):
        for K in (0, 2, 4, 6, 8):
            spins, labels = w5_input_key(K)
            self.assertEqual(spins, (4,) * 10)
            self.assertEqual(labels, (K, 0, 0, 0, 0))

    def test_rejects_non_w5_K(self):
        with self.assertRaises(ValueError):
            w5_input_key(10)

    def test_covariant_scalar_roundtrip_key(self):
        key = ((4,) * 10, (-1, 0, 2, 4, 6), 0, 0, 8, 8)
        out = covariant_scalar_to_gauss(key, 0)
        self.assertEqual(out, ((4,) * 10, (8, 0, 2, 4, 6)))

    def test_rejects_non_scalar(self):
        key = ((4,) * 10, (-1, 0, 0, 0, 0), 2, 0, 0, 2)
        with self.assertRaises(ValueError):
            covariant_scalar_to_gauss(key, 0)

    def test_epsilon_orbit_has_24_terms(self):
        neighbors = (1, 2, 3, 4)
        terms = ordered_epsilon_triples(neighbors)
        self.assertEqual(len(terms), 24)
        self.assertEqual(sum(sign for sign, _ in terms), 0)

if __name__ == "__main__":
    unittest.main()
