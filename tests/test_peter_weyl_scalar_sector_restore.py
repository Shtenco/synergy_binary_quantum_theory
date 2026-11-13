import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import peter_weyl_covariant_K_composition_cached_gate as gate


def key(j2=0, m2=0, node1_j2=0, node1_m2=0):
    spins = (4,) * 10
    labels = (
        ("C", j2, m2, 0, 0),
        ("C", node1_j2, node1_m2, 0, 0),
        ("G", 0), ("G", 0), ("G", 0),
    )
    return spins, labels


class TestScalarSectorNormContract(unittest.TestCase):
    def test_roundoff_residue_below_existing_squared_weight_gate_is_removed(self):
        state = {
            key(): 2.577039668 + 0j,
            key(j2=2, m2=0): 1.7136e-10 + 0j,
        }
        diagnostics = []
        restored = gate.restore_scalar_charge_M(
            state, (0,), {0: (0, 0)}, diagnostics=diagnostics
        )
        self.assertEqual(len(restored), 1)
        self.assertEqual(next(iter(restored.values())), 2.577039668 + 0j)
        self.assertLess(diagnostics[0]["0"]["wrong_JM_weight_fraction"], 1e-18)

    def test_material_wrong_J_sector_fails_closed(self):
        state = {key(): 1.0 + 0j, key(j2=2, m2=0): 1e-4 + 0j}
        with self.assertRaisesRegex(RuntimeError, "wrong_weight_fraction"):
            gate.restore_scalar_charge_M(state, (0,), {0: (0, 0)})

    def test_wrong_M_in_same_J_sector_fails_closed(self):
        state = {key(): 1.0 + 0j, key(j2=0, m2=2): 1e-4 + 0j}
        with self.assertRaisesRegex(RuntimeError, "wrong_weight_fraction"):
            gate.restore_scalar_charge_M(state, (0,), {0: (0, 0)})

    def test_each_charged_endpoint_is_checked_independently(self):
        state = {
            key(): 1.0 + 0j,
            key(node1_j2=2, node1_m2=0): 1e-4 + 0j,
        }
        with self.assertRaisesRegex(RuntimeError, "node 1"):
            gate.restore_scalar_charge_M(state, (0, 1), {0: (0, 0), 1: (0, 0)})

    def test_small_residual_on_second_endpoint_is_reported_and_removed(self):
        state = {
            key(): 1.0 + 0j,
            key(node1_j2=2, node1_m2=0): 1e-10 + 0j,
        }
        diagnostics = []
        restored = gate.restore_scalar_charge_M(
            state, (0, 1), {0: (0, 0), 1: (0, 0)}, diagnostics=diagnostics
        )
        self.assertEqual(len(restored), 1)
        self.assertIn("1", diagnostics[0])
        self.assertGreater(diagnostics[0]["1"]["wrong_JM_weight"], 0.0)

    def test_empty_output_is_valid_and_diagnostic(self):
        diagnostics = []
        restored = gate.restore_scalar_charge_M(
            {}, (0,), {0: (0, 0)}, diagnostics=diagnostics
        )
        self.assertEqual(restored, {})
        self.assertEqual(diagnostics[0]["0"]["wrong_JM_weight_fraction"], 0.0)

    def test_invalid_tolerance_is_rejected(self):
        for value in (-1.0, float("nan"), float("inf")):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    gate.restore_scalar_charge_M(
                        {key(): 1.0 + 0j}, (0,), {0: (0, 0)},
                        wrong_weight_fraction_tolerance=value,
                    )


if __name__ == "__main__":
    unittest.main(verbosity=2)
