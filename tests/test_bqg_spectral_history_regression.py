from __future__ import annotations
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_spectral_history_graph_gate as spectral


class SpectralHistoryRegressionTests(unittest.TestCase):
    def test_block_lanczos_regression_closes_exact_reachable_space(self):
        result = spectral.run_regression()
        self.assertTrue(result["passed"])
        self.assertTrue(result["closed"])
        self.assertEqual(result["reachable_dimension"], 7)
        self.assertEqual(result["reduced_dimension"], 7)
        self.assertLess(result["orthogonality_error"], 1e-9)
        self.assertLess(result["compression_error"], 1e-9)
        self.assertLess(result["zero_projector_error"], 1e-9)
        self.assertLess(result["spectral_weight_sum_rule_error"], 1e-9)
        self.assertLess(max(result["heat_kernel_errors"].values()), 1e-9)


if __name__ == "__main__":
    unittest.main()
