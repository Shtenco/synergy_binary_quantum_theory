import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bqg_v10_w5_lorentzian_adjoint_closure as V10


class TestV10StrictComplexAccumulation(unittest.TestCase):
    def setUp(self):
        self.key = ((4,) * 10, (-1, 0, 0, 0, 0, 0))
        self.value = 0.6e-12 + 0.6e-12j

    def cv(self, state, source_v, c, k, i, Jmax2):
        return {self.key: self.value}, 0.0

    def ck(self, state, source_v, a, i, j, Jmax2):
        return state, {}

    def test_forward_triple_accumulates_all_subthreshold_complex_paths(self):
        result, leakage, diagnostics = V10.triple_forward(
            {self.key: 1.0 + 0j}, 0, 0, 1, 2, 10, self.cv, self.ck
        )
        expected = 8 * self.value
        self.assertEqual(result[self.key], expected)
        self.assertEqual(leakage, 0.0)

    def test_adjoint_triple_accumulates_all_subthreshold_complex_paths(self):
        result, leakage, diagnostics = V10.triple_adjoint(
            {self.key: 1.0 + 0j}, 0, 0, 1, 2, 10, self.cv, self.ck
        )
        expected = 8 * self.value
        self.assertEqual(result[self.key], expected)
        self.assertEqual(leakage, 0.0)

    def test_strict_accumulator_is_used_by_production_v10_module(self):
        self.assertEqual(V10.SparseLinearAccumulator.__name__, "SparseLinearAccumulator")
        self.assertEqual(V10.SparseLinearAccumulator.__module__, "bqg_sparse_linear_accumulator")
        self.assertEqual(V10.TOL, 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
