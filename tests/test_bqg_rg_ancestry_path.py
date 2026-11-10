from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bqg_rg_ancestry_path as ancestry


class TestBQGRGAncestryPath(unittest.TestCase):
    def test_four_step_full_matrix_prefix(self):
        files = [
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json",
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_4.json",
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_5.json",
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_6.json",
        ]
        matrices = ancestry.load(files)
        result = ancestry.solve(matrices)

        self.assertEqual(result["channel_path"], [0, 1, 1, 1, 2])
        self.assertAlmostEqual(
            result["best_cumulative_fidelity"], 0.5685020686576171, places=12
        )
        self.assertAlmostEqual(
            result["naive_lowest_cumulative_fidelity"],
            4.732255132750176e-05,
            places=15,
        )
        self.assertGreater(
            result["ancestry_to_lowest_product_ratio"], 12000.0
        )

    def test_loader_rejects_nonconsecutive_scales(self):
        files = [
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json",
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_5.json",
        ]
        with self.assertRaisesRegex(ValueError, "nonconsecutive"):
            ancestry.load(files)


if __name__ == "__main__":
    unittest.main()
