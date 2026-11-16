import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_direct_all_m_shadow_equivalence import compare_directory


class TestDirectAllMShadowEquivalence(unittest.TestCase):
    def write_pair(self, directory, *, node, mode, channel, pair, perturb=0.0):
        k12, k34 = channel
        tag = "neg2_0" if pair == (-2, 0) else "0_2"
        columns = []
        for m2 in pair:
            columns.append({
                "source_node": node, "K12": k12, "K34": k34, "input_M2": m2,
                "entries": [
                    [f"out:{node}:{k12}:{k34}:{m2}", 1.0 + perturb, -0.25],
                    [f"outside:{node}:{k12}:{k34}:{m2}", 1e-14, 2e-14],
                ],
            })
        report = {
            "schema": "bqg-w5-direct-all-M-su2-commutator-v1",
            "status": "PASS_DIRECT_ALL_M_SU2_COMMUTATOR",
            "passed": True, "Jmax2": 10, "source_nodes": [node],
            "input_charged_sector": {"channels_K12_K34": [[k12, k34]]},
            "magnetic_components_tested_M2": list(pair),
            "pruning_mode": mode, "operator_columns_sparse": columns,
            "cases": [{"pruning_events": 3}, {"pruning_events": 4}],
        }
        (directory / (
            f"direct_all_M_{mode}_node{node}_channel{k12}_{k34}_mpair{tag}.json"
        )).write_text(json.dumps(report), encoding="utf-8")

    def make_complete_directory(self, directory, *, perturb=0.0):
        for node in range(5):
            for channel in ((0, 2), (2, 0)):
                for mode in ("production", "shadow_no_prune"):
                    for pair in ((-2, 0), (0, 2)):
                        self.write_pair(
                            directory, node=node, mode=mode, channel=channel,
                            pair=pair, perturb=perturb if mode == "shadow_no_prune" else 0.0,
                        )

    def test_compares_all_nodes_channels_and_all_m_columns(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.make_complete_directory(directory)
            result = compare_directory(directory)
        self.assertTrue(result["passed"], result)
        self.assertEqual(result["input_columns_compared"], 30)
        self.assertEqual(result["expected_input_columns"], 30)
        self.assertEqual(result["max_relative_l2_error"], 0.0)

    def test_detects_material_production_shadow_difference(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.make_complete_directory(directory, perturb=1e-4)
            result = compare_directory(directory)
        self.assertFalse(result["passed"])
        self.assertEqual(result["status"], "FAIL_PRODUCTION_SHADOW_EQUIVALENCE")
        self.assertGreater(result["max_relative_l2_error"], 1e-10)

    def test_missing_adjacent_m_shard_is_not_treated_as_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.make_complete_directory(directory)
            (directory / "direct_all_M_production_node0_channel0_2_mpair-2_0.json").unlink()
            with self.assertRaises(FileNotFoundError):
                compare_directory(directory)

    def test_disagreeing_overlapping_m0_columns_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.make_complete_directory(directory)
            path = directory / "direct_all_M_production_node0_channel0_2_mpair0_2.json"
            report = json.loads(path.read_text())
            report["operator_columns_sparse"][0]["entries"][0][1] += 0.1
            path.write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "overlapping M2=0"):
                compare_directory(directory)

    def test_failed_underlying_gate_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.make_complete_directory(directory)
            path = directory / "direct_all_M_production_node0_channel0_2_mpair-2_0.json"
            report = json.loads(path.read_text())
            report["passed"] = False
            report["status"] = "FAIL_DIRECT_ALL_M_SU2_COMMUTATOR"
            path.write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaises(ValueError):
                compare_directory(directory)


if __name__ == "__main__":
    unittest.main(verbosity=2)
