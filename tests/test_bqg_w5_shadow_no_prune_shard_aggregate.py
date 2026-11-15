import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_shadow_no_prune_shard_aggregate import aggregate


class TestW5ShadowShardAggregate(unittest.TestCase):
    def write_shard(self, directory, jmax2, node, *, valid=True):
        report = {
            "status": "INCOMPLETE_NO_PRUNING_REFERENCE",
            "source_nodes": [node],
            "Jmax2": jmax2,
            "comparison_performed": True,
            "parameters": {"unbounded_closure": True},
            "lower_level_pruning_bypassed_at_instrumented_sites": True,
            "all_lower_level_pruning_bypassed": False,
            "unverified_lower_level_sites_may_remain": True,
            "production_coefficient_loss_ledger": {},
            "shadow_coefficient_loss_ledger": {},
            "production_coherent_node_sum_ledger": [],
            "shadow_coherent_node_sum_ledger": [],
            "forward_matrix_comparison": {"relative_l2_residual": 0.0},
            "direct_adjoint_matrix_comparison": {"relative_l2_residual": 0.0},
        }
        if not valid:
            report["comparison_performed"] = False
        path = directory / f"W5_FULL_SHADOW_J{jmax2}_NODE{node}_RUN1.json"
        path.write_text(json.dumps(report), encoding="utf-8")

    def fill(self, directory):
        for j in (10, 12, 14):
            for node in range(5):
                self.write_shard(directory, j, node)

    def test_all_fifteen_shards_complete_but_claim_stays_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory)
        self.assertTrue(result["finite_matrix_comparison_complete"])
        self.assertFalse(result["passed"])
        self.assertEqual(result["status"], "INCOMPLETE_NO_PRUNING_COVERAGE_ALL_SHARDS_COMPLETED")
        self.assertEqual(result["valid_shard_count"], 15)

    def test_missing_shard_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_FULL_SHADOW_J14_NODE4_RUN1.json").unlink()
            result = aggregate(directory)
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertEqual(len(result["missing_shards"]), 1)

    def test_bad_comparison_cannot_be_aggregated(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            self.write_shard(directory, 10, 0, valid=False)
            result = aggregate(directory)
        self.assertFalse(result["finite_matrix_comparison_complete"])
        self.assertTrue(any("failed checks" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
