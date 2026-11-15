import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bqg_w5_cache_equivalence_shard_aggregate import aggregate


class TestW5CacheShardAggregate(unittest.TestCase):
    def write_shard(self, directory, jmax2, node, seed_k, *, passed=True):
        report = {
            "status": "PASS_CACHE_EQUIVALENCE" if passed else "INCOMPLETE_CACHE_EQUIVALENCE",
            "comparison_performed": passed,
            "cache_on": {"status": "PASS" if passed else "INCOMPLETE"},
            "cache_off": {"status": "PASS" if passed else "INCOMPLETE"},
            "comparison": {
                "status": "PASS_CACHE_EQUIVALENCE" if passed else "INCOMPLETE_CACHE_EQUIVALENCE",
                "comparison_completed": passed,
                "matrix_exactly_equal": passed,
                "equivalence_passed": passed,
                "observed_ledger_exactly_equal": passed,
                "static_coverage_complete": passed,
            },
            "parameters": {
                "nodes": [node],
                "seed_Ks": [seed_k],
                "Jmax2": jmax2,
                "unbounded_closure": True,
            },
        }
        path = directory / f"W5_CACHE_EQUIVALENCE_J{jmax2}_NODE{node}_K{seed_k}_RUN1.json"
        path.write_text(json.dumps(report), encoding="utf-8")

    def fill(self, directory, *, failed=None):
        for jmax2 in (10, 12, 14):
            for node in range(5):
                for seed_k in (0, 2, 4, 6, 8):
                    key = (jmax2, node, seed_k)
                    self.write_shard(directory, *key, passed=(key != failed))

    def test_all_75_exact_shards_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            result = aggregate(directory)
        self.assertTrue(result["passed"], result["errors"])
        self.assertEqual(result["valid_shard_count"], 75)
        self.assertEqual(result["expected_shard_count"], 75)

    def test_missing_shard_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            (directory / "W5_CACHE_EQUIVALENCE_J14_NODE4_K8_RUN1.json").unlink()
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["missing_shards"]), 1)

    def test_failed_cache_comparison_cannot_pass_aggregate(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory, failed=(10, 0, 0))
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertTrue(any("failed checks" in error for error in result["errors"]))

    def test_duplicate_shard_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            self.fill(directory)
            duplicate = directory / "W5_CACHE_EQUIVALENCE_J10_NODE0_K0_RUN2.json"
            duplicate.write_text((directory / "W5_CACHE_EQUIVALENCE_J10_NODE0_K0_RUN1.json").read_text())
            result = aggregate(directory)
        self.assertFalse(result["passed"])
        self.assertTrue(any("duplicate shard identity" in error for error in result["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
