import sys
import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_coverage_inventory as inventory


class TestV10TransitivePruningCoverageInventory(unittest.TestCase):
    def test_plain_lru_cache_is_detected_without_threshold_comparisons(self):
        import ast
        tree = ast.parse(
            "import functools\n@functools.lru_cache(None)\ndef cached():\n    return 1\n"
        )
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
        self.assertEqual(
            inventory._plain_lru_cache_decorators(function),
            ["functools.lru_cache(None)"],
        )

    def test_audit_lru_cache_is_not_reported_as_plain(self):
        import ast
        tree = ast.parse(
            "@audit_lru_cache(None)\ndef cached():\n    return 1\n"
        )
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef))
        self.assertEqual(inventory._plain_lru_cache_decorators(function), [])

    def test_inventory_detects_plain_cache_even_without_threshold_comparisons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "bqg_v10_w5_lorentzian_adjoint_closure.py").write_text(
                "import functools\n"
                "@functools.lru_cache(None)\n"
                "def cached_basis():\n"
                "    return (1, 2, 3)\n",
                encoding="utf-8",
            )
            with patch.object(inventory, "SCRIPTS", root), patch.object(
                inventory, "EXPECTED_MODULES",
                {"bqg_v10_w5_lorentzian_adjoint_closure"},
            ):
                result = inventory.run_inventory()
        self.assertFalse(result["complete"])
        self.assertEqual(len(result["plain_lru_caches"]), 1)
        self.assertEqual(
            result["plain_lru_caches"][0]["function"], "cached_basis"
        )

    def test_every_threshold_site_is_instrumented_or_explicitly_classified(self):
        result = inventory.run_inventory()
        self.assertTrue(result["expected_module_set_matches"], result["module_set_mismatch"])
        self.assertEqual(result["plain_lru_caches"], [])
        self.assertEqual(
            result["audited_lru_caches"],
            [{
                "module": "bqg_pruning_audit",
                "function": "cached",
                "decorator": "functools.lru_cache(maxsize=maxsize)",
                "classification": "audit_lru_cache_internal_value_and_trace_cache",
            }],
        )
        fixed_masks = [
            site for site in result["threshold_sites"]
            if site["module"] == "bqg_sparse_linear_accumulator"
            and site["function"] == "sparsify_matrix_columns"
        ]
        self.assertEqual(len(fixed_masks), 2)
        self.assertTrue(all(
            site["classification"] == "post_assembly_fixed_operator_mask_requires_error_bound"
            for site in fixed_masks
        ))
        self.assertEqual(result["unclassified_sites"], [])
        self.assertTrue(result["complete"])
        self.assertGreaterEqual(result["threshold_site_count"], 20)
        self.assertEqual(result["status"], "PASS_TRANSITIVE_THRESHOLD_INVENTORY")


if __name__ == "__main__":
    unittest.main(verbosity=2)
