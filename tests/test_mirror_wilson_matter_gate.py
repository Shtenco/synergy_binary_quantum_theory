from __future__ import annotations
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mirror_wilson_matter_gate as gate


class MirrorWilsonMatterGateTests(unittest.TestCase):
    def test_full_finite_gate_passes(self):
        result = gate.run()
        self.assertTrue(result["passed"])
        self.assertLess(result["mirror"]["mirror_rest_mass_difference"], gate.TOL)
        self.assertLess(result["mirror"]["mirror_rest_source_sum"], gate.TOL)
        self.assertLess(result["mirror"]["max_mirror_spectrum_difference"], gate.TOL)
        self.assertLess(result["mirror"]["max_moving_source_error"], gate.TOL)
        self.assertTrue(result["wilson_massless_corner_control"]["only_physical_corner_zero"])

    def test_mirror_partners_have_positive_equal_rest_mass_and_opposite_source(self):
        result = gate.run()["mirror"]
        plus, minus = result["branches"]["1"], result["branches"]["-1"]
        self.assertGreater(plus["rest_mass"], 0)
        self.assertGreater(minus["rest_mass"], 0)
        self.assertAlmostEqual(plus["rest_mass"], minus["rest_mass"], delta=gate.TOL)
        self.assertAlmostEqual(plus["rest_source"], -minus["rest_source"], delta=gate.TOL)


if __name__ == "__main__":
    unittest.main()
