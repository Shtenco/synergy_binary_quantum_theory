import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import k5_peter_weyl_safe_hda_column as PW
sys.path.insert(0, str(ROOT / "scripts" / "archive"))
import bqg_v09_w5_lorentzian_raw_column_v09_frozen as V09
import bqg_v10_w5_lorentzian_adjoint_closure as V10


class TestW5V09V10PairedRegression(unittest.TestCase):
    def test_all_five_seed_keys_are_identical(self):
        for K in (0, 2, 4, 6, 8):
            self.assertEqual(V09.w5_input_key(K), V10.w5_input_key(K))

    def test_all_five_nodes_have_identical_ordered_24_term_orbits(self):
        for node in range(5):
            neighbors = PW.NEIG[node]
            old = V09.ordered_epsilon_triples(neighbors)
            new = V10.epsilon_triples(neighbors)
            self.assertEqual(len(old), 24)
            self.assertEqual(len({p for _, p in old}), 24)
            self.assertEqual(old, new)
            self.assertEqual(sum(sign for sign, _ in old), 0)

    def test_v09_serializer_loss_is_explicit_and_complex_values_roundtrip(self):
        state = {
            ("large",): 1.25 - 0.75j,
            ("medium",): 2.0e-12 + 3.0e-12j,
            ("boundary",): 1.0e-13 + 0j,
            ("small",): 0.5e-13 - 0.5e-13j,
        }
        rows = V09.serialize_sparse(state)
        restored = {row[0]: complex(row[1], row[2]) for row in rows}
        self.assertEqual(len(rows), 2)
        self.assertIn(repr(("large",)), restored)
        self.assertIn(repr(("medium",)), restored)
        self.assertNotIn(repr(("boundary",)), restored)
        self.assertNotIn(repr(("small",)), restored)
        self.assertEqual(restored[repr(("large",))], state[("large",)])
        self.assertEqual(restored[repr(("medium",))], state[("medium",)])

    def test_v10_sparse_adjoint_residual_uses_full_complex_values(self):
        actual = [{("a",): 1.0 + 2.0j, ("b",): -3.0 + 0.5j}]
        expected = [{("a",): 1.0 + 2.0j, ("b",): -3.0 + 0.5j}]
        error, scale = V10.sparse_matrix_residual(actual, expected)
        self.assertEqual(error, 0.0)
        self.assertGreater(scale, 0.0)


if __name__ == "__main__":
    unittest.main()
