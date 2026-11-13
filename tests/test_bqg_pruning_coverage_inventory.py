import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_coverage_inventory as inventory


class TestV10TransitivePruningCoverageInventory(unittest.TestCase):
    def test_plain_lru_cache_is_detected_without_threshold_comparisons(self):
        import ast
        tree = ast.parse(
            "import functools\\n@functools.lru_cache(None)\\ndef cached():\\n    return 1\\n"
        )
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
        self.assertEqual(
            inventory._plain_lru_cache_decorators(function),
            ["functools.lru_cache(None)"],
        )

    def test_audit_lru_cache_is_not_reported_as_plain(self):
        import ast
        tree = ast.parse(
            "@audit_lru_cache(None)\\ndef cached():\\n    return 1\\n"
        )
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
        self.assertEqual(inventory._plain_lru_cache_decorators(function), [])

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
