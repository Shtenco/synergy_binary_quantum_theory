import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import peter_weyl_w5_direct_all_M_su2_gate as gate


class TestDirectAllMSU2FailClosed(unittest.TestCase):
    def test_empty_commutator_inputs_are_rejected(self):
        with self.assertRaisesRegex(AssertionError, "empty state"):
            gate.relative_error({}, {})

    def test_zero_norm_commutator_inputs_are_rejected(self):
        zero = {"state": 0j}
        with self.assertRaisesRegex(AssertionError, "zero or invalid norm"):
            gate.relative_error(zero, zero)

    def test_nonzero_output_is_required(self):
        with self.assertRaisesRegex(AssertionError, "empty sparse state"):
            gate.require_nonzero_output({}, label="synthetic H_E")
        with self.assertRaisesRegex(AssertionError, "nonzero norm"):
            gate.require_nonzero_output({"state": 0j}, label="synthetic H_E")

    def test_tiny_but_nonzero_operator_output_is_not_misclassified_as_empty(self):
        state = {"state": 1.0e-200 + 0j}
        self.assertEqual(gate.require_nonzero_output(state, label="tiny H_E"), 1.0e-200)
        rel, absolute, scale = gate.relative_error(state, state)
        self.assertEqual(rel, 0.0)
        self.assertEqual(absolute, 0.0)
        self.assertEqual(scale, 1.0e-200)

    def test_nonfinite_output_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-finite amplitude"):
            gate.require_nonzero_output({"state": complex(float("inf"), 0)}, label="H_E")

    def test_wrong_sector_fraction_uses_full_output_weight(self):
        state = {
            ((4,) * 10, (("C", 2, 0, 0, 2),)): 1.0 + 0j,
            ((4,) * 10, (("C", 4, 0, 0, 2),)): 1.0 + 0j,
        }
        fraction = gate.wrong_sector_fraction(
            state, node=0, expected_j2=2, expected_m2=0
        )
        self.assertEqual(fraction, 0.5)

    def test_full_gate_rejects_empty_production_operator(self):
        with patch.object(gate.SINE, "complete_HE_sine", return_value=({}, 0.0, 0.0)):
            with self.assertRaisesRegex(AssertionError, "H_E\(M2=-2\).*empty"):
                gate.run(Jmax2=10, pruning_mode="production")


if __name__ == "__main__":
    unittest.main(verbosity=2)
