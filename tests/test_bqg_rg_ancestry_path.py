from __future__ import annotations

import json
import tempfile
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bqg_rg_ancestry_path as ancestry


class TestBQGRGAncestryPath(unittest.TestCase):
    def test_five_step_full_matrix_ancestry(self):
        files = [
            ROOT / f"docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_{s2}.json"
            for s2 in (3, 4, 5, 6, 7)
        ]
        matrices = ancestry.load(files)
        result = ancestry.solve(matrices)

        self.assertEqual(result["channel_path"], [0, 1, 1, 1, 2, 2])
        self.assertAlmostEqual(
            result["best_cumulative_fidelity"], 0.49282787562801145, places=12
        )
        self.assertAlmostEqual(
            result["naive_lowest_cumulative_fidelity"],
            4.694570240506031e-05,
            places=15,
        )
        self.assertAlmostEqual(
            result["ancestry_to_lowest_product_ratio"],
            10497.827285142275,
            places=8,
        )

        saved = json.loads(
            (ROOT / "docs/numerical_continuation/BQG_RG_ANCESTRY_PATH_5STEP.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(saved["channel_path"], result["channel_path"])
        self.assertAlmostEqual(
            saved["best_cumulative_fidelity"],
            result["best_cumulative_fidelity"],
            places=12,
        )

    def test_loader_requires_explicit_pass(self):
        source = ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json"
        item = json.loads(source.read_text(encoding="utf-8"))
        item.pop("passed", None)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "missing_pass.json"
            path.write_text(json.dumps(item), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "explicitly pass"):
                ancestry.load([path])

    def test_loader_rejects_failed_shard(self):
        source = ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json"
        item = json.loads(source.read_text(encoding="utf-8"))
        item["passed"] = False
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "failed.json"
            path.write_text(json.dumps(item), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "explicitly pass"):
                ancestry.load([path])

    def test_loader_requires_single_half_spin_step(self):
        source = ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json"
        item = json.loads(source.read_text(encoding="utf-8"))
        item["target_s2"] = item["source_s2"] + 2
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wrong_step.json"
            path.write_text(json.dumps(item), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "advance s2 by one"):
                ancestry.load([path])

    def test_loader_rejects_nonfinite_master_spectrum(self):
        source = ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json"
        item = json.loads(source.read_text(encoding="utf-8"))
        item["source_master_eigenvalues"][0] = float("nan")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nonfinite_spectrum.json"
            path.write_text(json.dumps(item), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-finite"):
                ancestry.load([path])

    def test_loader_rejects_nonconsecutive_scales(self):
        files = [
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_3.json",
            ROOT / "docs/numerical_continuation/BQG_FULL_CHANNEL_FLOW_S2_5.json",
        ]
        with self.assertRaisesRegex(ValueError, "nonconsecutive"):
            ancestry.load(files)


if __name__ == "__main__":
    unittest.main()
