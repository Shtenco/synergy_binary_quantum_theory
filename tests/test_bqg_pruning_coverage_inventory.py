import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_coverage_inventory as inventory


class TestV10TransitivePruningCoverageInventory(unittest.TestCase):
    def test_every_threshold_site_is_instrumented_or_explicitly_classified(self):
        result = inventory.run_inventory()
        self.assertTrue(result["expected_module_set_matches"], result["module_set_mismatch"])
        self.assertEqual(result["plain_lru_caches"], [])
        self.assertEqual(result["unclassified_sites"], [])
        self.assertTrue(result["complete"])
        self.assertGreaterEqual(result["threshold_site_count"], 20)
        self.assertEqual(result["status"], "PASS_TRANSITIVE_THRESHOLD_INVENTORY")


if __name__ == "__main__":
    unittest.main(verbosity=2)
